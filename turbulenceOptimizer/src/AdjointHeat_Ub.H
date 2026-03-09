const bool optimizationCycleObjectiveAdjointMode =
    customTurbulenceActive && optimizationCycleTurbulentOptimizationMode;
const bool proModeHeatAdjointMode =
    customTurbulenceActive
 && turbulenceProModeActive
 && proModeAdjointStabilityPackActive;
const bool proModeSAHeatAdjointMode =
    proModeHeatAdjointMode
 && proModeSAAdjointClosureSelected;
const bool optimizationLikeHeatAdjointMode =
    optimizationCycleObjectiveAdjointMode || proModeHeatAdjointMode;
const bool optimizationCycleHeatRollbackGuardsEnabled =
    adjointRollbackGuardsEffective;
const volScalarField& thermalAdjointSourceField =
    optimizationCycleObjectiveAdjointMode ? optimizationCycleThermalAdjointField : Tb;
const label adjointHeatCorrectors =
    optimizationCycleObjectiveAdjointMode
  ? optimizationCycleHeatAdjointCorrectors
  : (proModeHeatAdjointMode ? proModeHeatAdjointCorrectors : 40);
const scalar optimizationCycleHeatAdjointMagnitudeLimit =
    optimizationCycleObjectiveAdjointMode
  ? optimizationCycleAdjointMaxMagnitude
  : (proModeHeatAdjointMode ? proModeAdjointMaxMagnitude : adjointCycleAdjointMaxMagnitude);
const bool optimizationCycleHeatNoImproveCheckEnabled =
    optimizationCycleHeatRollbackGuardsEnabled
 && (
      optimizationCycleObjectiveAdjointMode
      ? static_cast<bool>(optimizationCycleNoImproveCheck)
      : (proModeHeatAdjointMode ? static_cast<bool>(proModeAdjointNoImproveCheck) : static_cast<bool>(adjointCycleStagnationCheck))
    );
const label optimizationCycleHeatNoImprovePatience =
    optimizationCycleObjectiveAdjointMode
  ? optimizationCycleNoImprovePatience
  : (proModeHeatAdjointMode ? Foam::max(proModeAdjointNoImprovePatience, label(2)) : adjointCycleStagnationPatience);
const label optimizationCycleHeatNoImproveMinCorrectors =
    optimizationCycleObjectiveAdjointMode
  ? optimizationCycleNoImproveMinCorrectors
  : (proModeHeatAdjointMode ? Foam::max(proModeAdjointNoImproveMinCorrectors, label(3)) : label(1));
const scalar optimizationCycleHeatNoImproveTolerance =
    optimizationCycleObjectiveAdjointMode
  ? scalar(0.02)
  : (proModeHeatAdjointMode ? scalar(0.05) : scalar(0.02));
const scalar optimizationCycleHeatGrowthFactor =
    optimizationCycleObjectiveAdjointMode
  ? scalar(1.75)
  : (proModeHeatAdjointMode ? scalar(2.00) : scalar(1.75));
const bool resetHeatAdjointState =
    optimizationCycleObjectiveAdjointMode
 ? static_cast<bool>(optimizationCycleResetAdjointState)
 : (proModeHeatAdjointMode
    ? (static_cast<bool>(proModeResetAdjointState) || proModeForceAdjointResetNextIter)
    : false);
const bool useNuEffInUbEquation =
    optimizationCycleObjectiveAdjointMode
 || (proModeHeatAdjointMode && static_cast<bool>(proModeAdjointUseNuEff));
const bool useLegacyPbRelax =
    optimizationCycleObjectiveAdjointMode
  ? static_cast<bool>(adjointCycleUseLegacyPressureRelax)
  : (proModeHeatAdjointMode ? static_cast<bool>(proModeAdjointUseLegacyPressureRelax) : static_cast<bool>(adjointCycleUseLegacyPressureRelax));
const scalar proModeWarmupProgressUb =
    proModeHeatAdjointMode
  ? Foam::min
    (
      scalar(1.0),
      Foam::max
      (
        scalar(0.0),
        scalar(opt - 1)
       /Foam::max(scalar(proModeAdjointWarmupIterations - 1), scalar(1.0))
      )
    )
  : scalar(1.0);
const scalar proModeRollbackFreeProgressUb =
    proModeHeatAdjointMode
  ? Foam::min
    (
      scalar(1.0),
      scalar(proModeAdjointRollbackFreeStreak)
     /Foam::max(scalar(proModeAdjointRampRollbackFreeWindow), scalar(1.0))
    )
  : scalar(1.0);
const scalar proModeUbRamp =
    proModeHeatAdjointMode
  ? Foam::min
    (
      scalar(1.0),
      scalar(0.5)*(proModeWarmupProgressUb + proModeRollbackFreeProgressUb)
    )
  : scalar(1.0);
const scalar ubEqnRelaxFactor =
    proModeHeatAdjointMode
  ? proModeUbEqnRelaxMin + (proModeUbEqnRelaxMax - proModeUbEqnRelaxMin)*proModeUbRamp
  : scalar(1.0);
const scalar ubPseudoDecay =
    proModeHeatAdjointMode
  ? Foam::max
    (
      proModePseudoInvTauFloorRatio,
      scalar(1.0)
    - (scalar(1.0) - proModePseudoInvTauFloorRatio)*proModeWarmupProgressUb
    )
  : scalar(1.0);
const scalar ubPseudoInvTauValue =
    proModeHeatAdjointMode
  ? proModeUbPseudoInvTauBase*ubPseudoDecay
  : (optimizationCycleObjectiveAdjointMode ? 5.0 : 0.0);
const dimensionedScalar ubPseudoInvTau
(
  "ubPseudoInvTau",
  dimless/dimTime,
  ubPseudoInvTauValue
);
const scalar ubSourceSplitStrength =
    optimizationLikeHeatAdjointMode
  ? (
        proModeHeatAdjointMode
      ? (
            proModeSAHeatAdjointMode
          ? (0.95 - 0.40*proModeUbRamp)
          : (0.80 - 0.30*proModeUbRamp)
        )
      : 0.50
    )
  : 0.0;
const scalar ubThermalSourceRamp =
    proModeHeatAdjointMode
  ? (
        proModeSAHeatAdjointMode
      ? (0.10 + 0.90*proModeWarmupProgressUb)
      : (0.20 + 0.80*proModeWarmupProgressUb)
    )
  : scalar(1.0);
const bool useProModeRobustHeatLinearSolvers =
    proModeHeatAdjointMode
 && proModeRobustLinearSolverActive;

vectorField optimizationCycleBestUbInternal;
scalarField optimizationCycleBestPbInternal;
scalarField optimizationCycleBestPhibInternal;
scalar optimizationCycleBestCombinedMetric = GREAT;
scalar optimizationCycleBestUbMax = GREAT;
scalar optimizationCycleBestPbRange = GREAT;
label optimizationCycleNoImproveStreak = 0;
label optimizationCycleBestCorrector = 0;
label optimizationCycleAdjointCorrectorsUsed = 0;
word optimizationCycleHeatStopReason(optimizationLikeHeatAdjointMode ? "max-correctors" : "legacy");
optimizationCycleHeatAdjointRolledBack = false;
optimizationCycleHeatAdjointControlledStop = false;
optimizationCycleHeatAdjointRollbackCorrector = 0;
optimizationCycleHeatAdjointStopReason = optimizationCycleHeatStopReason;

if (optimizationLikeHeatAdjointMode)
{
  if (resetHeatAdjointState)
  {
    Ub = dimensionedVector("zeroOptimizationCycleUb", Ub.dimensions(), vector::zero);
    pb = dimensionedScalar("zeroOptimizationCyclePb", pb.dimensions(), 0.0);
    phib = dimensionedScalar("zeroOptimizationCyclePhib", phib.dimensions(), 0.0);
    Ub.correctBoundaryConditions();
    pb.correctBoundaryConditions();
  }

  optimizationCycleBestUbInternal = Ub.internalField();
  optimizationCycleBestPbInternal = pb.internalField();
  optimizationCycleBestPhibInternal = phib.internalField();

  if (Pstream::master())
  {
    std::ofstream debugAdjointCycleLogFileStream(debugAdjointCycleLogFile, std::ios::app);
    debugAdjointCycleLogFileStream
      << "\n"
      << (optimizationCycleObjectiveAdjointMode
          ? "[Optimization-Cycle Heat Adjoint Corrector Trace]\n"
          : "[Pro-Mode Heat Adjoint Corrector Trace]\n")
      << "  Iteration           : " << opt << "\n"
      << "  Config              : reset="
      << (resetHeatAdjointState ? "true" : "false")
      << " viscosity=" << (useNuEffInUbEquation ? "nuEff" : "nu")
      << " legacyPbRelax="
      << (useLegacyPbRelax ? "true" : "false")
      << " noImproveGuard="
      << (optimizationCycleHeatNoImproveCheckEnabled ? "true" : "false")
      << " minPass=" << optimizationCycleHeatNoImproveMinCorrectors
      << " patience=" << optimizationCycleHeatNoImprovePatience
      << " pseudoInvTau=" << ubPseudoInvTauValue
      << " UbRelax=" << ubEqnRelaxFactor
      << " sourceSplit=" << ubSourceSplitStrength
      << " sourceRamp=" << ubThermalSourceRamp
      << " |Ub|maxLimit=" << optimizationCycleHeatAdjointMagnitudeLimit
      << "\n";
  }
}

for (i = 0; i < adjointHeatCorrectors; i++)
{    
vectorField optimizationCycleUbPrevious;
scalarField optimizationCyclePbPrevious;
scalar optimizationCyclePbPassInitialResidual = 0.0;
scalar optimizationCyclePbPassFinalResidual = 0.0;
label optimizationCyclePbPassIterations = 0;
if (optimizationLikeHeatAdjointMode)
{
  optimizationCycleUbPrevious = Ub.internalField();
  optimizationCyclePbPrevious = pb.internalField();
}
const volVectorField UbOldState(Ub);

volVectorField adjointTransposeConvectionb((fvc::grad(U) & Ub));
const volVectorField ubThermalSource
(
  ubThermalSourceRamp*b1*thermalAdjointSourceField*fvc::grad(T)
);
volScalarField ubExplicitSourceMagnitude
(
  mag(adjointTransposeConvectionb) + mag(ubThermalSource)
);

tmp<fvVectorMatrix> tUbEqn
(
  useNuEffInUbEquation
? tmp<fvVectorMatrix>
  (
      new fvVectorMatrix
      (
          fvm::div(-phi, Ub)
        + adjointTransposeConvectionb
        - fvm::laplacian(nuEff, Ub)
        + fvm::Sp(alpha, Ub)
       ==
          fvOptions(Ub)
      )
  )
: tmp<fvVectorMatrix>
  (
      new fvVectorMatrix
      (
          fvm::div(-phi, Ub)
        + adjointTransposeConvectionb
        - fvm::laplacian(nu, Ub)
        + fvm::Sp(alpha, Ub)
       ==
          fvOptions(Ub)
      )
  )
);
fvVectorMatrix& UbEqn = tUbEqn.ref();

UbEqn += ubThermalSource;

if (optimizationLikeHeatAdjointMode)
{
  const volScalarField ubPatankarRef
  (
    Foam::max
    (
      mag(UbOldState),
      dimensionedScalar("ubPatankarRefFloor", Ub.dimensions(), 1.0)
    )
  );
  const volScalarField ubSourceDampingCoeff
  (
    ubSourceSplitStrength*ubExplicitSourceMagnitude/ubPatankarRef
  );
  UbEqn += fvm::Sp(ubPseudoInvTau + ubSourceDampingCoeff, Ub);
  UbEqn -=
      ubPseudoInvTau*UbOldState
    + ubSourceDampingCoeff*UbOldState;
}

if (proModeHeatAdjointMode)
{
  UbEqn.relax(ubEqnRelaxFactor);
}
else
{
  UbEqn.relax();
}

fvOptions.constrain(UbEqn);

tmp<fvVectorMatrix> tUbSolveEqn(UbEqn == -fvc::grad(pb));
fvVectorMatrix& UbSolveEqn = tUbSolveEqn.ref();
const auto UbLinearSolvePerformance =
  useProModeRobustHeatLinearSolvers
? UbSolveEqn.solve(proModeRobustNonsymmetricSolverControls)
: UbSolveEqn.solve();
recordLinearSolve
(
    UbLinearSolveStats,
    cmptMax(UbLinearSolvePerformance.initialResidual()),
    cmptMax(UbLinearSolvePerformance.finalResidual()),
    cmptMax(UbLinearSolvePerformance.nIterations())
);

fvOptions.correct(Ub);
//****************************************
volScalarField rAUb(1.0/UbEqn.A());
volVectorField HbyAb("HbyAb", Ub);
HbyAb = rAUb*UbEqn.H();
surfaceScalarField phiHbyAb("phiHbyAb", fvc::flux(HbyAb));
adjustPhi(phiHbyAb, Ub, pb);

tmp<volScalarField> rAtUb(rAUb);

 if (simple.consistent())
{
    rAtUb = 1.0/(1.0/rAUb - UbEqn.H1());
    phiHbyAb +=
        fvc::interpolate(rAtUb() - rAUb)*fvc::snGrad(pb)*mesh.magSf();
    HbyAb -= (rAUb - rAtUb())*fvc::grad(pb);
}
tUbEqn.clear();
// Update the pressure BCs to ensure flux consistency
constrainPressure(pb, Ub, phiHbyAb, rAtUb(), MRF);  
// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pbEqn
    (
        turbulenceProModeActive
      ? fvm::laplacian(rAUb, pb) == fvc::div(phiHbyAb)
      : fvm::laplacian(rAtUb(), pb) == fvc::div(phiHbyAb)
    );

    pbEqn.setReference(pbRefCell, pbRefValue);
    const auto pbLinearSolvePerformance =
      useProModeRobustHeatLinearSolvers
    ? pbEqn.solve(proModeRobustSymmetricSolverControls)
    : pbEqn.solve();
    recordLinearSolve
    (
        pbLinearSolveStats,
        pbLinearSolvePerformance.initialResidual(),
        pbLinearSolvePerformance.finalResidual(),
        pbLinearSolvePerformance.nIterations()
    );
    optimizationCyclePbPassInitialResidual =
      Foam::max(optimizationCyclePbPassInitialResidual, pbLinearSolvePerformance.initialResidual());
    optimizationCyclePbPassFinalResidual =
      Foam::max(optimizationCyclePbPassFinalResidual, pbLinearSolvePerformance.finalResidual());
    optimizationCyclePbPassIterations =
      Foam::max(optimizationCyclePbPassIterations, pbLinearSolvePerformance.nIterations());

    if (simple.finalNonOrthogonalIter())
    {
        phib = phiHbyAb - pbEqn.flux();
    }
}

//#include "adjointContinuityErrs.H"

if (!optimizationLikeHeatAdjointMode || useLegacyPbRelax)
{
  pb.relax();
}

// Adjoint momentum corrector
Ub =
  turbulenceProModeActive
? HbyAb - rAUb*fvc::grad(pb)
: HbyAb - rAtUb()*fvc::grad(pb);
Ub.correctBoundaryConditions();
fvOptions.correct(Ub);

if (optimizationLikeHeatAdjointMode)
{
  const vectorField& optimizationCycleUbCurrent = Ub.internalField();
  const scalarField& optimizationCyclePbCurrent = pb.internalField();
  scalar localUbDeltaSqSum = 0.0;
  scalar localUbScaleSqSum = 0.0;
  scalar localPbDeltaSqSum = 0.0;
  scalar localPbScaleSqSum = 0.0;
  scalar localPbSum = 0.0;
  scalar localPhibAbsSum = 0.0;
  scalar localPhibNet = 0.0;

  forAll(optimizationCyclePbCurrent, celli)
  {
    const scalar ubReferenceScale =
      Foam::max
      (
        Foam::max(Foam::mag(optimizationCycleUbPrevious[celli]), Foam::mag(optimizationCycleUbCurrent[celli])),
        scalar(1.0)
      );
    const scalar pbReferenceScale =
      Foam::max
      (
        Foam::max(Foam::mag(optimizationCyclePbPrevious[celli]), Foam::mag(optimizationCyclePbCurrent[celli])),
        scalar(1.0)
      );
    const vector ubDelta = optimizationCycleUbCurrent[celli] - optimizationCycleUbPrevious[celli];
    const scalar pbDelta = optimizationCyclePbCurrent[celli] - optimizationCyclePbPrevious[celli];

    localUbDeltaSqSum += magSqr(ubDelta);
    localUbScaleSqSum += sqr(ubReferenceScale);
    localPbDeltaSqSum += sqr(pbDelta);
    localPbScaleSqSum += sqr(pbReferenceScale);
    localPbSum += optimizationCyclePbCurrent[celli]*mesh.V()[celli];
  }

  const scalarField& optimizationCyclePhibInternal = phib.internalField();
  forAll(optimizationCyclePhibInternal, facei)
  {
    localPhibAbsSum += Foam::mag(optimizationCyclePhibInternal[facei]);
    localPhibNet += optimizationCyclePhibInternal[facei];
  }

  reduce(localUbDeltaSqSum, sumOp<scalar>());
  reduce(localUbScaleSqSum, sumOp<scalar>());
  reduce(localPbDeltaSqSum, sumOp<scalar>());
  reduce(localPbScaleSqSum, sumOp<scalar>());
  reduce(localPbSum, sumOp<scalar>());
  reduce(localPhibAbsSum, sumOp<scalar>());
  reduce(localPhibNet, sumOp<scalar>());

  optimizationCycleAdjointCorrectorsUsed = i + 1;
  const scalar optimizationCycleUNonlinearChange =
    Foam::sqrt(localUbDeltaSqSum/Foam::max(localUbScaleSqSum, SMALL));
  const scalar optimizationCyclePNonlinearChange =
    Foam::sqrt(localPbDeltaSqSum/Foam::max(localPbScaleSqSum, SMALL));
  const scalar optimizationCycleUTargetRatio =
    optimizationCycleUNonlinearChange/Foam::max(adjointCycleAdjointUNonlinearTol, SMALL);
  const scalar optimizationCyclePTargetRatio =
    optimizationCyclePNonlinearChange/Foam::max(adjointCycleAdjointPNonlinearTol, SMALL);
  const scalar optimizationCycleCombinedMetric = Foam::max(optimizationCycleUTargetRatio, optimizationCyclePTargetRatio);
  const scalar optimizationCycleBestMetricSafe =
    Foam::max(optimizationCycleBestCombinedMetric, SMALL);
  const volScalarField optimizationCycleUbMagnitude(mag(Ub));
  const scalar optimizationCycleUbMax = gMax(optimizationCycleUbMagnitude);
  const scalar optimizationCyclePbRange = gMax(pb) - gMin(pb);
  const scalar optimizationCyclePbMean = localPbSum/area;
  const scalar optimizationCyclePhibAbsSum = localPhibAbsSum;
  const scalar optimizationCyclePhibNet = localPhibNet;

  bool optimizationCycleIsBestPass(false);
  bool optimizationCyclePendingRollback(false);
  bool optimizationCycleMagnitudeGuardTriggered(false);
  bool optimizationCycleGrowthDiverged(false);
  bool optimizationCycleNoImproveDiverged(false);
  bool optimizationCycleNoImproveControlledStop(false);
  if
  (
      optimizationCycleHeatRollbackGuardsEnabled
   && (optimizationCycleUbMax != optimizationCycleUbMax || optimizationCycleUbMax > optimizationCycleHeatAdjointMagnitudeLimit)
  )
  {
    optimizationCycleMagnitudeGuardTriggered = true;
    optimizationCyclePendingRollback = true;
    optimizationCycleHeatStopReason = "magnitude-guard";
  }

  if
  (
      !optimizationCycleMagnitudeGuardTriggered
   && (
        optimizationCycleBestCorrector == 0
     || optimizationCycleCombinedMetric
        < (scalar(1.0) - optimizationCycleHeatNoImproveTolerance)*optimizationCycleBestMetricSafe
      )
  )
  {
    optimizationCycleBestCombinedMetric = optimizationCycleCombinedMetric;
    optimizationCycleBestUbMax = optimizationCycleUbMax;
    optimizationCycleBestPbRange = optimizationCyclePbRange;
    optimizationCycleBestUbInternal = optimizationCycleUbCurrent;
    optimizationCycleBestPbInternal = optimizationCyclePbCurrent;
    optimizationCycleBestPhibInternal = phib.internalField();
    optimizationCycleBestCorrector = optimizationCycleAdjointCorrectorsUsed;
    optimizationCycleNoImproveStreak = 0;
    optimizationCycleIsBestPass = true;
  }
  else if (!optimizationCycleMagnitudeGuardTriggered)
  {
    optimizationCycleNoImproveStreak++;
  }

  const bool optimizationCycleConverged =
  (
      (i + 1) >= 2
   && optimizationCycleUNonlinearChange <= adjointCycleAdjointUNonlinearTol
   && optimizationCyclePNonlinearChange <= adjointCycleAdjointPNonlinearTol
  );

  if
  (
      !optimizationCycleMagnitudeGuardTriggered
   && !optimizationCycleConverged
   && optimizationCycleHeatRollbackGuardsEnabled
   && optimizationCycleBestCorrector > 0
   && optimizationCycleAdjointCorrectorsUsed > optimizationCycleBestCorrector
   && optimizationCycleAdjointCorrectorsUsed >= optimizationCycleHeatNoImproveMinCorrectors
   && (
        optimizationCycleUbMax
        > optimizationCycleHeatGrowthFactor*Foam::max(optimizationCycleBestUbMax, scalar(1.0))
     || optimizationCyclePbRange
        > optimizationCycleHeatGrowthFactor*Foam::max(optimizationCycleBestPbRange, scalar(1.0))
      )
  )
  {
    optimizationCycleGrowthDiverged = true;
    optimizationCyclePendingRollback = true;
    optimizationCycleHeatStopReason = "growth-diverged";
  }

  if
  (
      !optimizationCycleMagnitudeGuardTriggered
   && !optimizationCycleGrowthDiverged
   && !optimizationCycleConverged
   && optimizationCycleHeatNoImproveCheckEnabled
   && optimizationCycleBestCorrector > 0
   && optimizationCycleAdjointCorrectorsUsed > optimizationCycleBestCorrector
   && optimizationCycleAdjointCorrectorsUsed >= optimizationCycleHeatNoImproveMinCorrectors
   && optimizationCycleNoImproveStreak >= optimizationCycleHeatNoImprovePatience
   && optimizationCycleCombinedMetric
      > (scalar(1.0) + optimizationCycleHeatNoImproveTolerance)*optimizationCycleBestMetricSafe
  )
  {
    optimizationCycleNoImproveDiverged = true;
    optimizationCycleNoImproveControlledStop = true;
    optimizationCyclePendingRollback = true;
    optimizationCycleHeatStopReason = "no-improve";
  }

  if (optimizationCycleConverged)
  {
    optimizationCycleHeatStopReason = "converged";
  }
  else if (!optimizationCyclePendingRollback && optimizationCycleAdjointCorrectorsUsed >= adjointHeatCorrectors)
  {
    optimizationCycleHeatStopReason = "max-correctors";
  }
  optimizationCycleHeatAdjointStopReason = optimizationCycleHeatStopReason;

  if (Pstream::master())
  {
    std::ofstream debugAdjointCycleLogFileStream(debugAdjointCycleLogFile, std::ios::app);
    auto fmtTraceScalar =
    [&](const scalar value)
    {
      std::ostringstream os;
      os << std::scientific << std::setprecision(3) << value;
      return os.str();
    };

    std::string optimizationCyclePassState("run");
    if (optimizationCycleConverged)
    {
      optimizationCyclePassState = "converged";
    }
    else if (optimizationCycleMagnitudeGuardTriggered || optimizationCycleGrowthDiverged)
    {
      optimizationCyclePassState = "rollback";
    }
    else if (optimizationCycleNoImproveControlledStop)
    {
      optimizationCyclePassState = "controlled-stop";
    }
    else if (optimizationCycleIsBestPass)
    {
      optimizationCyclePassState = "best";
    }

    debugAdjointCycleLogFileStream
      << "  pass=" << optimizationCycleAdjointCorrectorsUsed
      << " state=" << optimizationCyclePassState
      << " stop=" << optimizationCycleHeatStopReason
      << " UbInit=" << fmtTraceScalar(cmptMax(UbLinearSolvePerformance.initialResidual()))
      << " UbFinal=" << fmtTraceScalar(cmptMax(UbLinearSolvePerformance.finalResidual()))
      << " UbItr=" << cmptMax(UbLinearSolvePerformance.nIterations())
      << " pbInit=" << fmtTraceScalar(optimizationCyclePbPassInitialResidual)
      << " pbFinal=" << fmtTraceScalar(optimizationCyclePbPassFinalResidual)
      << " pbItr=" << optimizationCyclePbPassIterations
      << " dU=" << fmtTraceScalar(optimizationCycleUNonlinearChange)
      << " dP=" << fmtTraceScalar(optimizationCyclePNonlinearChange)
      << " |Ub|max=" << fmtTraceScalar(optimizationCycleUbMax)
      << " pbRange=" << fmtTraceScalar(optimizationCyclePbRange)
      << " pbMean=" << fmtTraceScalar(optimizationCyclePbMean)
      << " |phib|sum=" << fmtTraceScalar(optimizationCyclePhibAbsSum)
      << " phibNet=" << fmtTraceScalar(optimizationCyclePhibNet)
      << "\n";
  }

  if (optimizationCycleConverged)
  {
    break;
  }
  if (optimizationCyclePendingRollback)
  {
    Ub.primitiveFieldRef() = optimizationCycleBestUbInternal;
    pb.primitiveFieldRef() = optimizationCycleBestPbInternal;
    phib.primitiveFieldRef() = optimizationCycleBestPhibInternal;
    Ub.correctBoundaryConditions();
    pb.correctBoundaryConditions();
    fvOptions.correct(Ub);
    const bool optimizationCycleHardRollback =
        optimizationCycleMagnitudeGuardTriggered
     || optimizationCycleGrowthDiverged;
    if (optimizationCycleMagnitudeGuardTriggered)
    {
      optimizationCycleHeatStopReason = "magnitude-guard-rollback";
    }
    else if (optimizationCycleGrowthDiverged)
    {
      optimizationCycleHeatStopReason = "growth-diverged-rollback";
    }
    else if (optimizationCycleNoImproveDiverged)
    {
      optimizationCycleHeatStopReason =
        proModeHeatAdjointMode && proModeWarmupProgressUb < scalar(1.0)
      ? "no-improve-warmup-stop"
      : "no-improve-stop";
    }
    if (optimizationCycleHardRollback)
    {
      optimizationCycleHeatAdjointRolledBack = true;
      optimizationCycleHeatAdjointControlledStop = false;
      optimizationCycleHeatAdjointRollbackCorrector = optimizationCycleBestCorrector;
    }
    else if (optimizationCycleNoImproveControlledStop)
    {
      optimizationCycleHeatAdjointRolledBack = false;
      optimizationCycleHeatAdjointControlledStop = true;
      optimizationCycleHeatAdjointRollbackCorrector = optimizationCycleBestCorrector;
    }
    else
    {
      optimizationCycleHeatAdjointRolledBack = false;
      optimizationCycleHeatAdjointControlledStop = false;
      optimizationCycleHeatAdjointRollbackCorrector = 0;
    }
    optimizationCycleHeatAdjointStopReason = optimizationCycleHeatStopReason;
    break;
  }
}

Ub.storePrevIter();
pb.storePrevIter();
phib.storePrevIter();
if
(
    !optimizationCycleHeatAdjointRolledBack
 && !optimizationCycleHeatAdjointControlledStop
)
{
  optimizationCycleHeatAdjointRollbackCorrector = 0;
}
optimizationCycleHeatAdjointStopReason = optimizationCycleHeatStopReason;
}
