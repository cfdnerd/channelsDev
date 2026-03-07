const bool stageDObjectiveAdjointMode =
    customTurbulenceActive && stageDTurbulentOptimizationMode;
const volScalarField& thermalAdjointSourceField = stageDObjectiveAdjointMode ? Ta : Tb;
const label adjointHeatCorrectors = stageDObjectiveAdjointMode ? stageDHeatAdjointCorrectors : 40;
const scalar stageDHeatAdjointMagnitudeLimit =
    stageDObjectiveAdjointMode ? stageDAdjointMaxMagnitude : stageCAdjointMaxMagnitude;
const bool stageDHeatNoImproveCheckEnabled =
    stageDObjectiveAdjointMode
  ? static_cast<bool>(stageDNoImproveCheck)
  : static_cast<bool>(stageCStagnationCheck);
const label stageDHeatNoImprovePatience =
    stageDObjectiveAdjointMode ? stageDNoImprovePatience : stageCStagnationPatience;
const label stageDHeatNoImproveMinCorrectors =
    stageDObjectiveAdjointMode ? stageDNoImproveMinCorrectors : label(1);
const bool resetHeatAdjointState =
    stageDObjectiveAdjointMode && static_cast<bool>(stageDResetAdjointState);

vectorField stageDBestUbInternal;
scalarField stageDBestPbInternal;
scalarField stageDBestPhibInternal;
scalar stageDBestCombinedMetric = GREAT;
label stageDNoImproveStreak = 0;
label stageDBestCorrector = 0;
label stageDAdjointCorrectorsUsed = 0;
word stageDHeatStopReason(stageDObjectiveAdjointMode ? "max-correctors" : "legacy");

if (stageDObjectiveAdjointMode)
{
  if (resetHeatAdjointState)
  {
    Ub = dimensionedVector("zeroStageDUb", Ub.dimensions(), vector::zero);
    pb = dimensionedScalar("zeroStageDPb", pb.dimensions(), 0.0);
    phib = dimensionedScalar("zeroStageDPhib", phib.dimensions(), 0.0);
    Ub.correctBoundaryConditions();
    pb.correctBoundaryConditions();
  }

  stageDBestUbInternal = Ub.internalField();
  stageDBestPbInternal = pb.internalField();
  stageDBestPhibInternal = phib.internalField();

  if (Pstream::master())
  {
    std::ofstream debugStageCLogFileStream(debugStageCLogFile, std::ios::app);
    debugStageCLogFileStream
      << "\n[Stage-D Heat Adjoint Corrector Trace]\n"
      << "  Iteration           : " << opt << "\n"
      << "  Config              : reset="
      << (resetHeatAdjointState ? "true" : "false")
      << " viscosity=" << (stageDObjectiveAdjointMode ? "nuEff" : "nu")
      << " legacyPbRelax="
      << (stageCUseLegacyPressureRelax ? "true" : "false")
      << " noImproveGuard="
      << (stageDHeatNoImproveCheckEnabled ? "true" : "false")
      << " minPass=" << stageDHeatNoImproveMinCorrectors
      << " patience=" << stageDHeatNoImprovePatience
      << " |Ub|maxLimit=" << stageDHeatAdjointMagnitudeLimit
      << "\n";
  }
}

for (i = 0; i < adjointHeatCorrectors; i++)
{    
vectorField stageDUbPrevious;
scalarField stageDPbPrevious;
scalar stageDPbPassInitialResidual = 0.0;
scalar stageDPbPassFinalResidual = 0.0;
label stageDPbPassIterations = 0;
if (stageDObjectiveAdjointMode)
{
  stageDUbPrevious = Ub.internalField();
  stageDPbPrevious = pb.internalField();
}

volVectorField adjointTransposeConvectionb((fvc::grad(U) & Ub));

tmp<fvVectorMatrix> tUbEqn
(
  stageDObjectiveAdjointMode
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
        fvm::laplacian(rAtUb(), pb) == fvc::div(phiHbyAb)
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
    stageDPbPassInitialResidual =
      Foam::max(stageDPbPassInitialResidual, pbLinearSolvePerformance.initialResidual());
    stageDPbPassFinalResidual =
      Foam::max(stageDPbPassFinalResidual, pbLinearSolvePerformance.finalResidual());
    stageDPbPassIterations =
      Foam::max(stageDPbPassIterations, pbLinearSolvePerformance.nIterations());

    if (simple.finalNonOrthogonalIter())
    {
        phib = phiHbyAb - pbEqn.flux();
    }
}

//#include "adjointContinuityErrs.H"

if (!stageDObjectiveAdjointMode || stageCUseLegacyPressureRelax)
{
  pb.relax();
}

// Adjoint momentum corrector
Ub = HbyAb - rAtUb()*fvc::grad(pb);
Ub.correctBoundaryConditions();
fvOptions.correct(Ub);

if (stageDObjectiveAdjointMode)
{
  const vectorField& stageDUbCurrent = Ub.internalField();
  const scalarField& stageDPbCurrent = pb.internalField();
  scalar localUbDeltaSqSum = 0.0;
  scalar localUbScaleSqSum = 0.0;
  scalar localPbDeltaSqSum = 0.0;
  scalar localPbScaleSqSum = 0.0;
  scalar localPbSum = 0.0;
  scalar localPhibAbsSum = 0.0;
  scalar localPhibNet = 0.0;

  forAll(stageDPbCurrent, celli)
  {
    const scalar ubReferenceScale =
      Foam::max
      (
        Foam::max(Foam::mag(stageDUbPrevious[celli]), Foam::mag(stageDUbCurrent[celli])),
        scalar(1.0)
      );
    const scalar pbReferenceScale =
      Foam::max
      (
        Foam::max(Foam::mag(stageDPbPrevious[celli]), Foam::mag(stageDPbCurrent[celli])),
        scalar(1.0)
      );
    const vector ubDelta = stageDUbCurrent[celli] - stageDUbPrevious[celli];
    const scalar pbDelta = stageDPbCurrent[celli] - stageDPbPrevious[celli];

    localUbDeltaSqSum += magSqr(ubDelta);
    localUbScaleSqSum += sqr(ubReferenceScale);
    localPbDeltaSqSum += sqr(pbDelta);
    localPbScaleSqSum += sqr(pbReferenceScale);
    localPbSum += stageDPbCurrent[celli]*mesh.V()[celli];
  }

  const scalarField& stageDPhibInternal = phib.internalField();
  forAll(stageDPhibInternal, facei)
  {
    localPhibAbsSum += Foam::mag(stageDPhibInternal[facei]);
    localPhibNet += stageDPhibInternal[facei];
  }

  reduce(localUbDeltaSqSum, sumOp<scalar>());
  reduce(localUbScaleSqSum, sumOp<scalar>());
  reduce(localPbDeltaSqSum, sumOp<scalar>());
  reduce(localPbScaleSqSum, sumOp<scalar>());
  reduce(localPbSum, sumOp<scalar>());
  reduce(localPhibAbsSum, sumOp<scalar>());
  reduce(localPhibNet, sumOp<scalar>());

  stageDAdjointCorrectorsUsed = i + 1;
  const scalar stageDUNonlinearChange =
    Foam::sqrt(localUbDeltaSqSum/Foam::max(localUbScaleSqSum, SMALL));
  const scalar stageDPNonlinearChange =
    Foam::sqrt(localPbDeltaSqSum/Foam::max(localPbScaleSqSum, SMALL));
  const scalar stageDUTargetRatio =
    stageDUNonlinearChange/Foam::max(stageCAdjointUNonlinearTol, SMALL);
  const scalar stageDPTargetRatio =
    stageDPNonlinearChange/Foam::max(stageCAdjointPNonlinearTol, SMALL);
  const scalar stageDCombinedMetric = Foam::max(stageDUTargetRatio, stageDPTargetRatio);
  const volScalarField stageDUbMagnitude(mag(Ub));
  const scalar stageDUbMax = gMax(stageDUbMagnitude);
  const scalar stageDPbRange = gMax(pb) - gMin(pb);
  const scalar stageDPbMean = localPbSum/area;
  const scalar stageDPhibAbsSum = localPhibAbsSum;
  const scalar stageDPhibNet = localPhibNet;

  bool stageDIsBestPass(false);
  bool stageDPendingRollback(false);
  bool stageDMagnitudeGuardTriggered(false);
  bool stageDNoImproveDiverged(false);
  if (stageDUbMax != stageDUbMax || stageDUbMax > stageDHeatAdjointMagnitudeLimit)
  {
    stageDMagnitudeGuardTriggered = true;
    stageDPendingRollback = true;
    stageDHeatStopReason = "magnitude-guard";
  }

  if
  (
      !stageDMagnitudeGuardTriggered
   && (stageDBestCorrector == 0 || stageDCombinedMetric < stageDBestCombinedMetric)
  )
  {
    stageDBestCombinedMetric = stageDCombinedMetric;
    stageDBestUbInternal = stageDUbCurrent;
    stageDBestPbInternal = stageDPbCurrent;
    stageDBestPhibInternal = phib.internalField();
    stageDBestCorrector = stageDAdjointCorrectorsUsed;
    stageDNoImproveStreak = 0;
    stageDIsBestPass = true;
  }
  else if (!stageDMagnitudeGuardTriggered)
  {
    stageDNoImproveStreak++;
  }

  const bool stageDConverged =
  (
      (i + 1) >= 2
   && stageDUNonlinearChange <= stageCAdjointUNonlinearTol
   && stageDPNonlinearChange <= stageCAdjointPNonlinearTol
  );

  if
  (
      !stageDMagnitudeGuardTriggered
   && !stageDConverged
   && stageDHeatNoImproveCheckEnabled
   && stageDBestCorrector > 0
   && stageDAdjointCorrectorsUsed > stageDBestCorrector
   && stageDAdjointCorrectorsUsed >= stageDHeatNoImproveMinCorrectors
   && stageDNoImproveStreak >= stageDHeatNoImprovePatience
  )
  {
    stageDNoImproveDiverged = true;
    stageDPendingRollback = true;
    stageDHeatStopReason = "no-improve";
  }

  if (stageDConverged)
  {
    stageDHeatStopReason = "converged";
  }
  else if (!stageDPendingRollback && stageDAdjointCorrectorsUsed >= adjointHeatCorrectors)
  {
    stageDHeatStopReason = "max-correctors";
  }

  if (Pstream::master())
  {
    std::ofstream debugStageCLogFileStream(debugStageCLogFile, std::ios::app);
    auto fmtTraceScalar =
    [&](const scalar value)
    {
      std::ostringstream os;
      os << std::scientific << std::setprecision(3) << value;
      return os.str();
    };

    std::string stageDPassState("run");
    if (stageDConverged)
    {
      stageDPassState = "converged";
    }
    else if (stageDMagnitudeGuardTriggered || stageDNoImproveDiverged)
    {
      stageDPassState = "rollback";
    }
    else if (stageDIsBestPass)
    {
      stageDPassState = "best";
    }

    debugStageCLogFileStream
      << "  pass=" << stageDAdjointCorrectorsUsed
      << " state=" << stageDPassState
      << " stop=" << stageDHeatStopReason
      << " UbInit=" << fmtTraceScalar(cmptMax(UbLinearSolvePerformance.initialResidual()))
      << " UbFinal=" << fmtTraceScalar(cmptMax(UbLinearSolvePerformance.finalResidual()))
      << " UbItr=" << cmptMax(UbLinearSolvePerformance.nIterations())
      << " pbInit=" << fmtTraceScalar(stageDPbPassInitialResidual)
      << " pbFinal=" << fmtTraceScalar(stageDPbPassFinalResidual)
      << " pbItr=" << stageDPbPassIterations
      << " dU=" << fmtTraceScalar(stageDUNonlinearChange)
      << " dP=" << fmtTraceScalar(stageDPNonlinearChange)
      << " |Ub|max=" << fmtTraceScalar(stageDUbMax)
      << " pbRange=" << fmtTraceScalar(stageDPbRange)
      << " pbMean=" << fmtTraceScalar(stageDPbMean)
      << " |phib|sum=" << fmtTraceScalar(stageDPhibAbsSum)
      << " phibNet=" << fmtTraceScalar(stageDPhibNet)
      << "\n";
  }

  if (stageDConverged)
  {
    break;
  }
  if (stageDPendingRollback)
  {
    Ub.primitiveFieldRef() = stageDBestUbInternal;
    pb.primitiveFieldRef() = stageDBestPbInternal;
    phib.primitiveFieldRef() = stageDBestPhibInternal;
    Ub.correctBoundaryConditions();
    pb.correctBoundaryConditions();
    fvOptions.correct(Ub);
    if (stageDMagnitudeGuardTriggered)
    {
      stageDHeatStopReason = "magnitude-guard-rollback";
    }
    else if (stageDNoImproveDiverged)
    {
      stageDHeatStopReason = "no-improve-rollback";
    }
    break;
  }
}

Ub.storePrevIter();
pb.storePrevIter();
phib.storePrevIter();
}
