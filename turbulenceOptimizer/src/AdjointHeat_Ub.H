const bool optimizationCycleObjectiveAdjointMode =
    customTurbulenceActive && optimizationCycleTurbulentOptimizationMode;
const bool turbulenceProHeatAdjointMode =
    customTurbulenceActive && turbulenceProModeActive;
const bool optimizationCycleHeatRollbackGuardsEnabled =
    static_cast<bool>(optAdjointRollbackGuardsEnable);
const volScalarField& thermalAdjointSourceField =
    optimizationCycleObjectiveAdjointMode ? optimizationCycleThermalAdjointField : Tb;
const label adjointHeatCorrectors =
    (optimizationCycleObjectiveAdjointMode || turbulenceProHeatAdjointMode)
  ? optimizationCycleHeatAdjointCorrectors
  : 40;
const scalar optimizationCycleHeatAdjointMagnitudeLimit =
    (optimizationCycleObjectiveAdjointMode || turbulenceProHeatAdjointMode)
  ? optimizationCycleAdjointMaxMagnitude
  : adjointCycleAdjointMaxMagnitude;
const bool optimizationCycleHeatNoImproveCheckEnabled =
    optimizationCycleHeatRollbackGuardsEnabled
 && (
      (optimizationCycleObjectiveAdjointMode || turbulenceProHeatAdjointMode)
      ? static_cast<bool>(optimizationCycleNoImproveCheck)
      : static_cast<bool>(adjointCycleStagnationCheck)
    );
const label optimizationCycleHeatNoImprovePatience =
    (optimizationCycleObjectiveAdjointMode || turbulenceProHeatAdjointMode)
  ? optimizationCycleNoImprovePatience
  : adjointCycleStagnationPatience;
const label optimizationCycleHeatNoImproveMinCorrectors =
    (optimizationCycleObjectiveAdjointMode || turbulenceProHeatAdjointMode)
  ? optimizationCycleNoImproveMinCorrectors
  : label(1);
const bool resetHeatAdjointState =
    optimizationCycleObjectiveAdjointMode && static_cast<bool>(optimizationCycleResetAdjointState);

vectorField optimizationCycleBestUbInternal;
scalarField optimizationCycleBestPbInternal;
scalarField optimizationCycleBestPhibInternal;
scalar optimizationCycleBestCombinedMetric = GREAT;
label optimizationCycleNoImproveStreak = 0;
label optimizationCycleBestCorrector = 0;
label optimizationCycleAdjointCorrectorsUsed = 0;
word optimizationCycleHeatStopReason(optimizationCycleObjectiveAdjointMode ? "max-correctors" : "legacy");

if (optimizationCycleObjectiveAdjointMode)
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
      << "\n[Optimization-Cycle Heat Adjoint Corrector Trace]\n"
      << "  Iteration           : " << opt << "\n"
      << "  Config              : reset="
      << (resetHeatAdjointState ? "true" : "false")
      << " viscosity=" << (optimizationCycleObjectiveAdjointMode ? "nuEff" : "nu")
      << " legacyPbRelax="
      << (adjointCycleUseLegacyPressureRelax ? "true" : "false")
      << " noImproveGuard="
      << (optimizationCycleHeatNoImproveCheckEnabled ? "true" : "false")
      << " minPass=" << optimizationCycleHeatNoImproveMinCorrectors
      << " patience=" << optimizationCycleHeatNoImprovePatience
      << " |Ub|maxLimit=" << optimizationCycleHeatAdjointMagnitudeLimit
      << "\n";
  }
}

if (turbulenceProModeActive)
{
  // Prevent Pro-mode warm-start amplification across optimization iterations.
  Ub = dimensionedVector("zeroProUb", Ub.dimensions(), vector::zero);
  pb = dimensionedScalar("zeroProPb", pb.dimensions(), 0.0);
  phib = dimensionedScalar("zeroProPhib", phib.dimensions(), 0.0);
  Ub.correctBoundaryConditions();
  pb.correctBoundaryConditions();
}

for (i = 0; i < adjointHeatCorrectors; i++)
{    
vectorField optimizationCycleUbPrevious;
scalarField optimizationCyclePbPrevious;
scalar optimizationCyclePbPassInitialResidual = 0.0;
scalar optimizationCyclePbPassFinalResidual = 0.0;
label optimizationCyclePbPassIterations = 0;
if (optimizationCycleObjectiveAdjointMode)
{
  optimizationCycleUbPrevious = Ub.internalField();
  optimizationCyclePbPrevious = pb.internalField();
}

volVectorField adjointTransposeConvectionb((fvc::grad(U) & Ub));

tmp<fvVectorMatrix> tUbEqn
(
  optimizationCycleObjectiveAdjointMode
? tmp<fvVectorMatrix>
  (
      new fvVectorMatrix
      (
          fvm::div(-phi, Ub)
        + adjointTransposeConvectionb
        - fvm::laplacian(nuEff, Ub)
        + fvm::Sp(alpha, Ub)
        + b1*thermalAdjointSourceField*fvc::grad(T)
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
        + b1*thermalAdjointSourceField*fvc::grad(T)
       ==
          fvOptions(Ub)
      )
  )
);
fvVectorMatrix& UbEqn = tUbEqn.ref();

UbEqn.relax();

fvOptions.constrain(UbEqn);

const auto UbLinearSolvePerformance = solve(UbEqn == -fvc::grad(pb));
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
    const auto pbLinearSolvePerformance = pbEqn.solve();
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

if (!optimizationCycleObjectiveAdjointMode || adjointCycleUseLegacyPressureRelax)
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

if (turbulenceProModeActive)
{
  const volScalarField proUbMagnitude(mag(Ub));
  const scalar ubAbsMax = gMax(proUbMagnitude);
  if (!(ubAbsMax < optimizationCycleAdjointMaxMagnitude))
  {
    const scalar proAdjointScale =
      optimizationCycleAdjointMaxMagnitude/Foam::max(ubAbsMax, SMALL);
    Ub.primitiveFieldRef() *= proAdjointScale;
    pb.primitiveFieldRef() *= proAdjointScale;
    phib.primitiveFieldRef() *= proAdjointScale;
    Ub.correctBoundaryConditions();
    pb.correctBoundaryConditions();
    fvOptions.correct(Ub);
  }
}

if (optimizationCycleObjectiveAdjointMode)
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
  const volScalarField optimizationCycleUbMagnitude(mag(Ub));
  const scalar optimizationCycleUbMax = gMax(optimizationCycleUbMagnitude);
  const scalar optimizationCyclePbRange = gMax(pb) - gMin(pb);
  const scalar optimizationCyclePbMean = localPbSum/area;
  const scalar optimizationCyclePhibAbsSum = localPhibAbsSum;
  const scalar optimizationCyclePhibNet = localPhibNet;

  bool optimizationCycleIsBestPass(false);
  bool optimizationCyclePendingRollback(false);
  bool optimizationCycleMagnitudeGuardTriggered(false);
  bool optimizationCycleNoImproveDiverged(false);
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
   && (optimizationCycleBestCorrector == 0 || optimizationCycleCombinedMetric < optimizationCycleBestCombinedMetric)
  )
  {
    optimizationCycleBestCombinedMetric = optimizationCycleCombinedMetric;
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
   && optimizationCycleHeatNoImproveCheckEnabled
   && optimizationCycleBestCorrector > 0
   && optimizationCycleAdjointCorrectorsUsed > optimizationCycleBestCorrector
   && optimizationCycleAdjointCorrectorsUsed >= optimizationCycleHeatNoImproveMinCorrectors
   && optimizationCycleNoImproveStreak >= optimizationCycleHeatNoImprovePatience
  )
  {
    optimizationCycleNoImproveDiverged = true;
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
    else if (optimizationCycleMagnitudeGuardTriggered || optimizationCycleNoImproveDiverged)
    {
      optimizationCyclePassState = "rollback";
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
    if (optimizationCycleMagnitudeGuardTriggered)
    {
      optimizationCycleHeatStopReason = "magnitude-guard-rollback";
    }
    else if (optimizationCycleNoImproveDiverged)
    {
      optimizationCycleHeatStopReason = "no-improve-rollback";
    }
    break;
  }
}

Ub.storePrevIter();
pb.storePrevIter();
phib.storePrevIter();
}
