//***********************************//update mu and lambda
opt++;
const label updateScheduleIterationCurrent =
  Foam::max(opt, label(1));
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
const label updateScheduleHorizonLagIterations =
  Foam::max
  (
    updatePowerConstraintScheduleIterationsToSettle - updateOptimizationIterationsRemaining,
    label(0)
  );
(void)updateScheduleHorizonLagIterations;
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

const scalar recoveryContinuationRampScale =
  Foam::max(Foam::min(optRecoveryContinuationRampScale, scalar(1.0)), scalar(0.10));
const bool recoverySignal =
  maxViolationWorsenStreak >= optRecoveryWorsenStreak;
const bool earlyRecoveryAllowed = false;
const bool recoveryIterationGate =
  opt >= optConstraintWorsenMinIterations || earlyRecoveryAllowed;
const bool recoveryEligible =
  recoverySignal
  && recoveryIterationGate;
constraintRecoveryActive =
  static_cast<bool>(optRecoveryEnable) && recoveryEligible;
continuationRampScaleActive =
  constraintRecoveryActive ? recoveryContinuationRampScale : scalar(1.0);
delRampStepActive =
  constraintRecoveryActive ? scalar(0.05) : scalar(0.2);
optimizerMoveLimitActive = movlim;
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
