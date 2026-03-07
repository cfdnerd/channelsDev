//***********************************//update mu and lambda
opt++;
const label updateScheduleIterationCurrent =
  Foam::max(continuationScheduleIterationActive, label(1));
const scalar updatePowerConstraintScheduleDenomCurrent =
  Foam::max(PowerDissRelax - scalar(0.2)*updateScheduleIterationCurrent, PowerDissMax);
const scalar updatePowerConstraintScheduleGap =
  Foam::max(updatePowerConstraintScheduleDenomCurrent - PowerDissMax, scalar(0));
const scalar updatePowerConstraintScheduleStep = 0.2;
const label updatePowerConstraintScheduleIterationsToSettle =
  updatePowerConstraintScheduleGap > SMALL
  ? static_cast<label>(updatePowerConstraintScheduleGap/updatePowerConstraintScheduleStep + 0.999999)
  : 0;
const scalar updateOptimizationIterationsRemainingRaw =
  Foam::max(runTime.endTime().value() - static_cast<scalar>(opt), scalar(0));
const label updateOptimizationIterationsRemaining =
  static_cast<label>(updateOptimizationIterationsRemainingRaw + 0.000001);
const bool laminarReferenceWarmupActive =
  opt <= optLaminarReferenceWarmupIterations;
const label updateScheduleHorizonLagIterations =
  Foam::max
  (
    updatePowerConstraintScheduleIterationsToSettle - updateOptimizationIterationsRemaining,
    label(0)
  );
const bool updatePowerConstraintScheduleSettled =
  updatePowerConstraintScheduleDenomCurrent <= PowerDissMax + SMALL;
const bool updatePowerConstraintScheduleWillSettleInRun =
  updatePowerConstraintScheduleIterationsToSettle <= updateOptimizationIterationsRemaining;
const bool updateContinuationFreezePhaseActive =
  updatePowerConstraintScheduleSettled || !updatePowerConstraintScheduleWillSettleInRun;
continuationFrozenByReject = false;
continuationFrozenByLowTrust = false;
const scalar continuationClipHoldNoProgressTolClamped =
  Foam::max(optContinuationClipHoldNoProgressTol, scalar(1e-8));
continuationFrozenByClipHold =
  static_cast<bool>(optContinuationClipHoldEnable)
  && updateContinuationFreezePhaseActive
  && (
       (
            hasPreviousClipDiagnostics
         && previousMaxFluidClipFraction >= optContinuationClipHoldFluidFraction
         && previousGuardConstraintViolation > optConstraintTol
         && previousGuardConstraintViolationDelta > -continuationClipHoldNoProgressTolClamped
       )
    || previousHealthGateFailed
    || previousThermalEnvelopeViolation > 0.0
     )
  && !continuationAdvanceOverrideOnReject;
const bool continuationFrozen =
  continuationFrozenByReject || continuationFrozenByLowTrust || continuationFrozenByClipHold;
if (continuationIteration <= 0)
{
  continuationIteration = 1;
}
else if (!continuationFrozen)
{
  const bool continuationHealthAcceptable =
       !previousHealthGateFailed
    && previousThermalEnvelopeViolation <= 0.0
    && !sensitivityQualityHardFail
    && !sensitivityQualityRequireConservativeStep;
  const bool continuationAdvanceAllowed =
    !static_cast<bool>(optContinuationEventDriven)
    || (
         true
      && continuationHealthAcceptable
       );
  if (continuationAdvanceAllowed)
  {
    continuationIteration += 1 + Foam::max(continuationAdvanceCatchUpSteps, label(0));
  }
}
continuationScheduleIterationActive =
  continuationIteration;
continuationAdvanceOverrideOnReject = false;
continuationAdvanceCatchUpSteps = 0;

const scalar recoveryMoveLimitScale =
  Foam::max(Foam::min(optRecoveryMoveLimitScale, scalar(1.0)), scalar(0.05));
const scalar recoveryContinuationRampScale =
  Foam::max(Foam::min(optRecoveryContinuationRampScale, scalar(1.0)), scalar(0.10));
const scalar globalizationTrustScale = scalar(1.0);
const scalar adaptiveBaseMoveLower =
  Foam::max(Foam::min(optAdaptiveBaseMoveMin, scalar(1.0)), scalar(0.01));
const scalar adaptiveBaseMoveUpper =
  Foam::max(Foam::min(optAdaptiveBaseMoveMax, scalar(1.0)), adaptiveBaseMoveLower);
const scalar previousClipFractionForMove =
  hasPreviousClipDiagnostics
  ? Foam::max(previousMaxClipFraction, previousMaxFluidClipFraction)
  : scalar(0.0);
const scalar clipHealthScale =
  hasPreviousClipDiagnostics
  ? Foam::max
    (
      scalar(0.15),
      scalar(1.0) - Foam::min(previousClipFractionForMove, scalar(0.85))
    )
  : scalar(1.0);
const scalar sensitivityHealthScale =
  laminarReferenceWarmupActive
  ? scalar(1.0)
  : Foam::max(Foam::min(sensitivityQualityScore, scalar(1.0)), scalar(0.0));
const scalar multiphysicsHealthScale =
  Foam::max(Foam::min(Foam::min(globalizationTrustScale, sensitivityHealthScale), clipHealthScale), scalar(0.0));
adaptiveMmaBaseMoveActive =
  adaptiveBaseMoveLower
  + (adaptiveBaseMoveUpper - adaptiveBaseMoveLower)*multiphysicsHealthScale;
adaptiveMmaAsymIncActive =
  optAdaptiveAsymIncMin
  + (optAdaptiveAsymIncMax - optAdaptiveAsymIncMin)*multiphysicsHealthScale;
adaptiveMmaAsymDecActive =
  optAdaptiveAsymDecMax
  - (optAdaptiveAsymDecMax - optAdaptiveAsymDecMin)*multiphysicsHealthScale;
const bool recoverySignal =
     maxViolationWorsenStreak >= optRecoveryWorsenStreak
  || meritRecoveryRequest;
const bool earlyRecoveryAllowed = false;
const bool recoveryIterationGate =
  opt >= optConstraintWorsenMinIterations || earlyRecoveryAllowed;
const bool recoveryEligible =
  recoverySignal
  && recoveryIterationGate
  && !laminarReferenceWarmupActive;
constraintRecoveryActive =
  static_cast<bool>(optRecoveryEnable) && recoveryEligible;
continuationRampScaleActive =
  constraintRecoveryActive ? recoveryContinuationRampScale : scalar(1.0);
delRampStepActive =
  constraintRecoveryActive ? scalar(0.05) : scalar(0.2);
optimizerMoveLimitActive =
  constraintRecoveryActive
  ? Foam::max(adaptiveMmaBaseMoveActive*recoveryMoveLimitScale*globalizationTrustScale, optMmaMoveLimitFloor)
  : Foam::max(adaptiveMmaBaseMoveActive*globalizationTrustScale, optMmaMoveLimitFloor);
if (sensitivityQualityRequireConservativeStep)
{
  optimizerMoveLimitActive =
    Foam::max
    (
      optimizerMoveLimitActive*optSensitivityLowQualityMoveScale,
      optMmaMoveLimitFloor
    );
}
if (laminarReferenceWarmupActive)
{
  optimizerMoveLimitActive =
    Foam::max
    (
      optimizerMoveLimitActive,
      Foam::max(optLaminarReferenceWarmupMoveLimit, optMmaMoveLimitFloor)
    );
}

primalThermalUnderRelaxActive =
  optThermalUnderRelaxMin
  + (optThermalUnderRelaxMax - optThermalUnderRelaxMin)*multiphysicsHealthScale;
primalTurbulenceUnderRelaxActive =
  optTurbulenceUnderRelaxMin
  + (optTurbulenceUnderRelaxMax - optTurbulenceUnderRelaxMin)*multiphysicsHealthScale;
primalPseudoTimeCoeffActive =
  Foam::max(scalar(0.0), scalar(1.0) - multiphysicsHealthScale)
  *optPrimalPseudoTimeCoeffMax;
thermalSourceRampActive =
  optThermalSourceRampMin
  + (scalar(1.0) - optThermalSourceRampMin)*multiphysicsHealthScale;
if(solid_area)
{
   forAll(cells_solid, i)
   {
      xh[cells_solid[i]]=0;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      xh[cells_fluid[i]]=1.0;
   }
}
xh.correctBoundaryConditions();

const scalar continuationAlphaBase = alphamax.value();
const scalar continuationAlphaAt63 = continuationAlphaBase*(63.0/7.0 + 1.0);
scalar scheduledAlphaMaxValue = continuationAlphaAt63;
if (continuationScheduleIterationActive <= 63)
{
  scheduledAlphaMaxValue =
    continuationAlphaBase*(continuationScheduleIterationActive/7.0 + 1.0);
}
else if
(
    continuationScheduleIterationActive > 100
 && !static_cast<bool>(optDisableLateStageContinuationAcceleration)
)
{
  const scalar continuationGrowthSteps =
    scalar(continuationScheduleIterationActive - 100);
  const scalar continuationGrowthBase =
    1.0 + 0.05*continuationRampScaleActive;
  const scalar continuationGrowthFactor =
    Foam::pow(continuationGrowthBase, continuationGrowthSteps);
  scheduledAlphaMaxValue = continuationAlphaAt63*continuationGrowthFactor;
}
alphaMax.value() = Foam::min(scheduledAlphaMaxValue, alphaMAX.value());

const scalar scheduledQu =
  continuationScheduleIterationActive >= 80
  ? 0.005 + (continuationScheduleIterationActive-80)*1e-4*continuationRampScaleActive
  : 0.005;
qu = scheduledQu;
alphaMax.value()=Foam::min(alphaMax.value(),alphaMAX.value());
qu=Foam::min(qu,0.01);
alpha=alphaMax*qu*(1-xh)/(qu+xh);
DT=(ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoCp;
