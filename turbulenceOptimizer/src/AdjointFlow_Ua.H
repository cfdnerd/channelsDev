const bool adjointCycleFlowAdjointMode = customTurbulenceActive && adjointCycleFullAdjointVerification;
const bool optimizationCycleFlowAdjointMode =
    customTurbulenceActive && optimizationCycleTurbulentOptimizationMode;
const bool proModeFlowAdjointMode =
    customTurbulenceActive
 && turbulenceProModeActive
 && proModeAdjointStabilityPackActive;
const bool optimizationLikeFlowAdjointMode =
    optimizationCycleFlowAdjointMode || proModeFlowAdjointMode;
const bool boundedFlowAdjointMode = adjointCycleFlowAdjointMode || optimizationLikeFlowAdjointMode;
const bool adjointRollbackGuardsEnabled =
    adjointRollbackGuardsEffective;
const bool resetAdjointState =
    (adjointCycleFlowAdjointMode && static_cast<bool>(adjointCycleResetAdjointState))
 || (optimizationCycleFlowAdjointMode && static_cast<bool>(optimizationCycleResetAdjointState))
 || (proModeFlowAdjointMode && (static_cast<bool>(proModeResetAdjointState) || proModeForceAdjointResetNextIter));
const bool stagnationCheckEnabled =
    adjointRollbackGuardsEnabled
 && (
      (adjointCycleFlowAdjointMode && static_cast<bool>(adjointCycleStagnationCheck))
   || (optimizationCycleFlowAdjointMode && static_cast<bool>(optimizationCycleNoImproveCheck))
   || (proModeFlowAdjointMode && static_cast<bool>(proModeAdjointNoImproveCheck))
    );
const label noImprovePatience =
    optimizationCycleFlowAdjointMode
  ? optimizationCycleNoImprovePatience
  : (proModeFlowAdjointMode ? proModeAdjointNoImprovePatience : adjointCycleStagnationPatience);
const label noImproveMinCorrectors =
    optimizationCycleFlowAdjointMode
  ? optimizationCycleNoImproveMinCorrectors
  : (proModeFlowAdjointMode ? proModeAdjointNoImproveMinCorrectors : label(1));
const bool divergenceCheckEnabled =
    adjointRollbackGuardsEnabled
 && adjointCycleFlowAdjointMode
 && adjointCycleDivergenceCheck;
const scalar flowAdjointMagnitudeLimit =
    optimizationCycleFlowAdjointMode
  ? optimizationCycleAdjointMaxMagnitude
  : (proModeFlowAdjointMode ? proModeAdjointMaxMagnitude : adjointCycleAdjointMaxMagnitude);
const bool useNuEffInUaEquation =
    optimizationCycleFlowAdjointMode
 || (proModeFlowAdjointMode && static_cast<bool>(proModeAdjointUseNuEff))
 || (adjointCycleFlowAdjointMode && adjointCycleUseNuEff);
const bool includeTaSourceInUaEquation =
    optimizationCycleFlowAdjointMode || (adjointCycleFlowAdjointMode && adjointCycleIncludeThermalMomentumSource);
const bool useLegacyPressureRelax =
    optimizationCycleFlowAdjointMode
  ? static_cast<bool>(adjointCycleUseLegacyPressureRelax)
  : (proModeFlowAdjointMode ? static_cast<bool>(proModeAdjointUseLegacyPressureRelax) : static_cast<bool>(adjointCycleUseLegacyPressureRelax));
const label adjointFlowCorrectors =
    optimizationCycleFlowAdjointMode
  ? optimizationCycleFlowAdjointCorrectors
  : (proModeFlowAdjointMode
    ? proModeFlowAdjointCorrectors
    : (adjointCycleFlowAdjointMode ? adjointCycleAdjointCorrectors : 40));
const scalar proModeWarmupProgressUa =
    proModeFlowAdjointMode
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
const scalar proModeRollbackFreeProgressUa =
    proModeFlowAdjointMode
  ? Foam::min
    (
      scalar(1.0),
      scalar(proModeAdjointRollbackFreeStreak)
     /Foam::max(scalar(proModeAdjointRampRollbackFreeWindow), scalar(1.0))
    )
  : scalar(1.0);
const scalar proModeUaRamp =
    proModeFlowAdjointMode
  ? Foam::min
    (
      scalar(1.0),
      scalar(0.5)*(proModeWarmupProgressUa + proModeRollbackFreeProgressUa)
    )
  : scalar(1.0);
const scalar uaEqnRelaxFactor =
    proModeFlowAdjointMode
  ? proModeUaEqnRelaxMin + (proModeUaEqnRelaxMax - proModeUaEqnRelaxMin)*proModeUaRamp
  : scalar(1.0);
const scalar uaPseudoDecay =
    proModeFlowAdjointMode
  ? Foam::max
    (
      Foam::max(proModePseudoInvTauFloorRatio, scalar(0.40)),
      scalar(1.0)
    - (scalar(1.0) - proModePseudoInvTauFloorRatio)*proModeWarmupProgressUa
    )
  : scalar(1.0);
const scalar uaPseudoInvTauValue =
    proModeFlowAdjointMode
  ? proModeUaPseudoInvTauBase*uaPseudoDecay
  : (optimizationCycleFlowAdjointMode ? 4.0 : (adjointCycleFlowAdjointMode ? 2.0 : 0.0));
const dimensionedScalar uaPseudoInvTau
(
  "uaPseudoInvTau",
  dimless/dimTime,
  uaPseudoInvTauValue
);
const scalar uaSourceSplitStrength =
    boundedFlowAdjointMode
  ? (
      proModeFlowAdjointMode
      ? (0.55 - 0.20*proModeUaRamp)
      : 0.30
    )
  : 0.0;
const bool useProModeRobustFlowLinearSolvers =
    proModeFlowAdjointMode
 && proModeRobustLinearSolverActive;
adjointCycleAdjointUNonlinearChange = 0.0;
adjointCycleAdjointPNonlinearChange = 0.0;
adjointCycleAdjointFlowConverged = !boundedFlowAdjointMode;
adjointCycleAdjointCorrectorsUsed = 0;
adjointCycleAdjointCombinedMetric = 0.0;
adjointCycleAdjointBestCombinedMetric = GREAT;
adjointCycleAdjointBestUNonlinearChange = GREAT;
adjointCycleAdjointBestPNonlinearChange = GREAT;
adjointCycleAdjointBestUaMax = GREAT;
adjointCycleAdjointBestPaRange = GREAT;
adjointCycleAdjointBestPaMean = 0.0;
adjointCycleAdjointBestPhiaAbsSum = 0.0;
adjointCycleAdjointBestPhiaNet = 0.0;
adjointCycleAdjointBestCorrector = 0;
adjointCycleAdjointUaMax = 0.0;
adjointCycleAdjointPaRange = 0.0;
adjointCycleAdjointPaMean = 0.0;
adjointCycleAdjointPhiaAbsSum = 0.0;
adjointCycleAdjointPhiaNet = 0.0;
adjointCycleAdjointFlowDiverged = false;
adjointCycleAdjointDivergedCorrector = 0;
adjointCycleAdjointRolledBack = false;
adjointCycleAdjointFlowControlledStop = false;
adjointCycleAdjointRollbackCorrector = 0;
adjointCycleAdjointStopReason = boundedFlowAdjointMode ? "max-correctors" : "not-bounded";

vectorField adjointCycleBestUaInternal;
scalarField adjointCycleBestPaInternal;
scalarField adjointCycleBestPhiaInternal;

if (boundedFlowAdjointMode)
{
  if (resetAdjointState)
  {
    // Default Adjoint-Cycle behavior: zero the adjoint state so each verification
    // sweep is isolated. Warm-starting remains available as an A/B toggle.
    Ua = dimensionedVector("zeroAdjointCycleUa", Ua.dimensions(), vector::zero);
    pa = dimensionedScalar("zeroAdjointCyclePa", pa.dimensions(), 0.0);
    phia = dimensionedScalar("zeroAdjointCyclePhia", phia.dimensions(), 0.0);
    Ua.correctBoundaryConditions();
    pa.correctBoundaryConditions();
    if (proModeFlowAdjointMode && proModeForceAdjointResetNextIter)
    {
      proModeForceAdjointResetNextIter = false;
    }
  }

  if (Pstream::master())
  {
      std::ofstream debugAdjointCycleLogFileStream(debugAdjointCycleLogFile, std::ios::app);
      const char* adjointCycleTraceTag =
          optimizationCycleFlowAdjointMode
        ? "[Optimization-Cycle Flow Adjoint Corrector Trace]"
        : (proModeFlowAdjointMode
          ? "[Pro-Mode Flow Adjoint Corrector Trace]"
          : "[Adjoint-Cycle Flow Adjoint Corrector Trace]");
      debugAdjointCycleLogFileStream
        << "\n" << adjointCycleTraceTag << "\n"
        << "  Iteration           : " << opt << "\n"
        << "  Config              : reset="
        << (resetAdjointState ? "true" : "false")
        << " viscosity=" << (useNuEffInUaEquation ? "nuEff" : "nu")
        << " TaSource="
        << (includeTaSourceInUaEquation ? "true" : "false")
        << " legacyPaRelax="
        << (useLegacyPressureRelax ? "true" : "false")
        << " divergenceGuard="
        << (divergenceCheckEnabled ? "true" : "false")
        << " stagnationGuard="
        << (stagnationCheckEnabled ? "true" : "false")
        << " |Ua|maxLimit=" << flowAdjointMagnitudeLimit
        << " pseudoInvTau=" << uaPseudoInvTauValue
        << " UaRelax=" << uaEqnRelaxFactor
        << " sourceSplit=" << uaSourceSplitStrength
        << " growth=" << adjointCycleDivergenceGrowthFactor
        << " minPass="
        << (optimizationLikeFlowAdjointMode ? noImproveMinCorrectors : adjointCycleDivergenceMinCorrectors)
        << " patience=" << noImprovePatience
        << "\n"
        << "  "
        << std::left
        << std::setw(6) << "Pass"
        << std::setw(12) << "State"
        << std::setw(14) << "UaInit"
        << std::setw(14) << "UaFinal"
        << std::setw(8) << "UaItr"
        << std::setw(14) << "paInit"
        << std::setw(14) << "paFinal"
        << std::setw(8) << "paItr"
        << std::setw(14) << "dU"
        << std::setw(14) << "dP"
        << std::setw(14) << "|Ua|max"
        << std::setw(14) << "paRange"
        << std::setw(14) << "paMean"
        << std::setw(14) << "|phia|sum"
        << std::setw(14) << "phiaNet"
        << "\n"
        << "  "
        << std::string(6 + 12 + 14 + 14 + 8 + 14 + 14 + 8 + 14 + 14 + 14 + 14 + 14 + 14 + 14, '-')
        << "\n";
  }

  adjointCycleBestUaInternal = Ua.internalField();
  adjointCycleBestPaInternal = pa.internalField();
  adjointCycleBestPhiaInternal = phia.internalField();
}

label adjointCycleNoImproveStreak = 0;
vectorField adjointCycleUaDeltaPrevious;
bool adjointCycleUaHasDeltaPrevious = false;
scalar adjointCycleUaAitkenOmega = 1.0;

for (i = 0; i < adjointFlowCorrectors; i++)
{
  vectorField UaPrevious;
  scalarField paPrevious;
  scalar adjointCyclePaPassInitialResidual = 0.0;
  scalar adjointCyclePaPassFinalResidual = 0.0;
  label adjointCyclePaPassIterations = 0;
  if (boundedFlowAdjointMode)
  {
      UaPrevious = Ua.internalField();
      paPrevious = pa.internalField();
  }
  const volVectorField UaOldState(Ua);

  volVectorField adjointTransposeConvectiona((fvc::grad(U) & Ua));
  volScalarField uaExplicitSourceMagnitude(mag(adjointTransposeConvectiona));

  tmp<fvVectorMatrix> tUaEqn
  (
    useNuEffInUaEquation
  ? tmp<fvVectorMatrix>
    (
        new fvVectorMatrix
        (
            fvm::div(-phi, Ua)
          + adjointTransposeConvectiona
          - fvm::laplacian(nuEff, Ua)
          + fvm::Sp(alpha, Ua)
         ==
            fvOptions(Ua)
        )
    )
  : tmp<fvVectorMatrix>
    (
        new fvVectorMatrix
        (
            fvm::div(-phi, Ua)
          + adjointTransposeConvectiona
          - fvm::laplacian(nu, Ua)
          + fvm::Sp(alpha, Ua)
         ==
            fvOptions(Ua)
        )
    )
  );
  fvVectorMatrix& UaEqn = tUaEqn.ref();

  if (includeTaSourceInUaEquation)
  {
      const volVectorField uaThermalSource
      (
        b1*optimizationCycleThermalAdjointField*fvc::grad(T)
      );
      UaEqn += uaThermalSource;
      uaExplicitSourceMagnitude += mag(uaThermalSource);
  }

  if (boundedFlowAdjointMode)
  {
      const volScalarField uaPatankarRef
      (
        Foam::max
        (
          mag(UaOldState),
          dimensionedScalar("uaPatankarRefFloor", Ua.dimensions(), 1.0)
        )
      );
      const volScalarField uaSourceDampingCoeff
      (
        uaSourceSplitStrength*uaExplicitSourceMagnitude/uaPatankarRef
      );
      UaEqn += fvm::Sp(uaPseudoInvTau + uaSourceDampingCoeff, Ua);
      UaEqn -=
          uaPseudoInvTau*UaOldState
        + uaSourceDampingCoeff*UaOldState;
  }

  if (proModeFlowAdjointMode)
  {
      UaEqn.relax(uaEqnRelaxFactor);
  }
  else
  {
      UaEqn.relax();
  }

  fvOptions.constrain(UaEqn);

  tmp<fvVectorMatrix> tUaSolveEqn(UaEqn == -fvc::grad(pa));
  fvVectorMatrix& UaSolveEqn = tUaSolveEqn.ref();
  const auto UaLinearSolvePerformance =
      useProModeRobustFlowLinearSolvers
    ? UaSolveEqn.solve(proModeRobustNonsymmetricSolverControls)
    : UaSolveEqn.solve();
  recordLinearSolve
  (
      UaLinearSolveStats,
      cmptMax(UaLinearSolvePerformance.initialResidual()),
      cmptMax(UaLinearSolvePerformance.finalResidual()),
      cmptMax(UaLinearSolvePerformance.nIterations())
  );

  fvOptions.correct(Ua);

  volScalarField rAUa(1.0/UaEqn.A());
  volVectorField HbyAa("HbyAa", Ua);
  HbyAa = rAUa*UaEqn.H();
  surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
  adjustPhi(phiHbyAa, Ua, pa);

  tmp<volScalarField> rAtUa(rAUa);

  if (simple.consistent())
  {
      rAtUa = 1.0/(1.0/rAUa - UaEqn.H1());
      phiHbyAa +=
          fvc::interpolate(rAtUa() - rAUa)*fvc::snGrad(pa)*mesh.magSf();
      HbyAa -= (rAUa - rAtUa())*fvc::grad(pa);
  }

  tUaEqn.clear();

  constrainPressure(pa, Ua, phiHbyAa, rAtUa(), MRF);
  while (simple.correctNonOrthogonal())
  {
      fvScalarMatrix paEqn
      (
          turbulenceProModeActive
        ? fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
        : fvm::laplacian(rAtUa(), pa) == fvc::div(phiHbyAa)
      );

      paEqn.setReference(paRefCell, paRefValue);
      const auto paLinearSolvePerformance =
          useProModeRobustFlowLinearSolvers
        ? paEqn.solve(proModeRobustSymmetricSolverControls)
        : paEqn.solve();
      recordLinearSolve
      (
          paLinearSolveStats,
          paLinearSolvePerformance.initialResidual(),
          paLinearSolvePerformance.finalResidual(),
          paLinearSolvePerformance.nIterations()
      );
      adjointCyclePaPassInitialResidual =
        Foam::max(adjointCyclePaPassInitialResidual, paLinearSolvePerformance.initialResidual());
      adjointCyclePaPassFinalResidual =
        Foam::max(adjointCyclePaPassFinalResidual, paLinearSolvePerformance.finalResidual());
      adjointCyclePaPassIterations =
        Foam::max(adjointCyclePaPassIterations, paLinearSolvePerformance.nIterations());

      if (simple.finalNonOrthogonalIter())
      {
          phia = phiHbyAa - paEqn.flux();
      }
  }

  if (!boundedFlowAdjointMode || useLegacyPressureRelax)
  {
      pa.relax();
  }

  Ua =
    turbulenceProModeActive
  ? HbyAa - rAUa*fvc::grad(pa)
  : HbyAa - rAtUa()*fvc::grad(pa);
Ua.correctBoundaryConditions();
fvOptions.correct(Ua);
if
(
    boundedFlowAdjointMode
 && static_cast<bool>(optAdjointAitkenEnable)
 && (i + 1) >= optAdjointAitkenStartCorrector
)
{
    const vectorField adjointCycleUaCurrentRaw(Ua.internalField());
    vectorField adjointCycleUaDelta(adjointCycleUaCurrentRaw - UaPrevious);
    if (adjointCycleUaHasDeltaPrevious)
    {
        scalar adjointCycleAitkenNumerator = 0.0;
        scalar adjointCycleAitkenDenominator = 0.0;
        forAll(adjointCycleUaDelta, celli)
        {
            const vector adjointCycleUaDeltaDiff =
              adjointCycleUaDelta[celli] - adjointCycleUaDeltaPrevious[celli];
            adjointCycleAitkenNumerator +=
              (adjointCycleUaDeltaPrevious[celli] & adjointCycleUaDeltaDiff);
            adjointCycleAitkenDenominator +=
              (adjointCycleUaDeltaDiff & adjointCycleUaDeltaDiff);
        }
        reduce(adjointCycleAitkenNumerator, sumOp<scalar>());
        reduce(adjointCycleAitkenDenominator, sumOp<scalar>());
        if (adjointCycleAitkenDenominator > optAdjointAitkenEpsilon)
        {
            scalar adjointCycleAitkenOmegaNew =
              -adjointCycleUaAitkenOmega
              *adjointCycleAitkenNumerator
              /Foam::max(adjointCycleAitkenDenominator, optAdjointAitkenEpsilon);
            adjointCycleAitkenOmegaNew =
              Foam::max
              (
                Foam::min(adjointCycleAitkenOmegaNew, optAdjointAitkenOmegaMax),
                optAdjointAitkenOmegaMin
              );
            vectorField& adjointCycleUaInternal = Ua.primitiveFieldRef();
            forAll(adjointCycleUaInternal, celli)
            {
                adjointCycleUaInternal[celli] =
                  UaPrevious[celli] + adjointCycleAitkenOmegaNew*adjointCycleUaDelta[celli];
            }
            Ua.correctBoundaryConditions();
            fvOptions.correct(Ua);
            adjointCycleUaAitkenOmega = adjointCycleAitkenOmegaNew;
            adjointCycleUaDelta = Ua.internalField() - UaPrevious;
        }
    }
    adjointCycleUaDeltaPrevious = adjointCycleUaDelta;
    adjointCycleUaHasDeltaPrevious = true;
}
else if (boundedFlowAdjointMode)
{
    adjointCycleUaHasDeltaPrevious = false;
    adjointCycleUaAitkenOmega = 1.0;
}

if (boundedFlowAdjointMode)
{
      const vectorField& UaCurrent = Ua.internalField();
      const scalarField& paCurrent = pa.internalField();
      scalar localUaDeltaSqSum = 0.0;
      scalar localUaScaleSqSum = 0.0;
      scalar localPaDeltaSqSum = 0.0;
      scalar localPaScaleSqSum = 0.0;
      scalar localPaSum = 0.0;
      scalar localPhiaAbsSum = 0.0;
      scalar localPhiaNet = 0.0;

      forAll(paCurrent, celli)
      {
          const scalar uaReferenceScale =
            Foam::max
            (
              Foam::max(Foam::mag(UaPrevious[celli]), Foam::mag(UaCurrent[celli])),
              scalar(1.0)
            );
          const scalar paReferenceScale =
            Foam::max
            (
              Foam::max(Foam::mag(paPrevious[celli]), Foam::mag(paCurrent[celli])),
              scalar(1.0)
            );
          const vector uaDelta = UaCurrent[celli] - UaPrevious[celli];
          const scalar paDelta = paCurrent[celli] - paPrevious[celli];

          localUaDeltaSqSum += magSqr(uaDelta);
          localUaScaleSqSum += sqr(uaReferenceScale);
          localPaDeltaSqSum += sqr(paDelta);
          localPaScaleSqSum += sqr(paReferenceScale);
          localPaSum += paCurrent[celli]*mesh.V()[celli];
      }

      const scalarField& phiaInternal = phia.internalField();
      forAll(phiaInternal, facei)
      {
          localPhiaAbsSum += Foam::mag(phiaInternal[facei]);
          localPhiaNet += phiaInternal[facei];
      }

      reduce(localUaDeltaSqSum, sumOp<scalar>());
      reduce(localUaScaleSqSum, sumOp<scalar>());
      reduce(localPaDeltaSqSum, sumOp<scalar>());
      reduce(localPaScaleSqSum, sumOp<scalar>());
      reduce(localPaSum, sumOp<scalar>());
      reduce(localPhiaAbsSum, sumOp<scalar>());
      reduce(localPhiaNet, sumOp<scalar>());

      adjointCycleAdjointUNonlinearChange =
        Foam::sqrt(localUaDeltaSqSum/Foam::max(localUaScaleSqSum, SMALL));
      adjointCycleAdjointPNonlinearChange =
        Foam::sqrt(localPaDeltaSqSum/Foam::max(localPaScaleSqSum, SMALL));
      adjointCycleAdjointCorrectorsUsed = i + 1;
      const volScalarField adjointCycleUaMag(mag(Ua));
      const scalar adjointCycleUaMax = gMax(adjointCycleUaMag);
      const scalar adjointCyclePaRange = gMax(pa) - gMin(pa);
      adjointCycleAdjointUaMax = adjointCycleUaMax;
      adjointCycleAdjointPaRange = adjointCyclePaRange;
      adjointCycleAdjointPaMean = localPaSum/area;
      adjointCycleAdjointPhiaAbsSum = localPhiaAbsSum;
      adjointCycleAdjointPhiaNet = localPhiaNet;
      const scalar adjointCycleUaTargetRatio =
        adjointCycleAdjointUNonlinearChange/Foam::max(adjointCycleAdjointUNonlinearTol, SMALL);
      const scalar adjointCyclePaTargetRatio =
        adjointCycleAdjointPNonlinearChange/Foam::max(adjointCycleAdjointPNonlinearTol, SMALL);
      const char* adjointCycleCurrentBlocker =
          (adjointCycleUaTargetRatio > 1.1*adjointCyclePaTargetRatio)
        ? "Ua"
        : (adjointCyclePaTargetRatio > 1.1*adjointCycleUaTargetRatio)
          ? "pa"
          : "Ua+pa";
      adjointCycleAdjointCombinedMetric =
        Foam::max
        (
            adjointCycleUaTargetRatio,
            adjointCyclePaTargetRatio
        );

      bool adjointCycleIsBestPass = false;
      bool adjointCyclePendingRollback(false);

      if (adjointCycleUaMax != adjointCycleUaMax)
      {
          FatalErrorInFunction
              << "Adjoint-Cycle Ua adjoint magnitude became invalid (|Ua|max="
              << adjointCycleUaMax
              << "). This indicates verification-loop divergence."
              << exit(FatalError);
      }

      bool adjointCycleMagnitudeGuardTriggered(false);
      if (adjointRollbackGuardsEnabled && adjointCycleUaMax > flowAdjointMagnitudeLimit)
      {
          adjointCycleMagnitudeGuardTriggered = true;
          adjointCycleAdjointFlowDiverged = true;
          adjointCycleAdjointDivergedCorrector = adjointCycleAdjointCorrectorsUsed;
          adjointCycleAdjointStopReason = "magnitude-guard";
          adjointCyclePendingRollback = true;
      }

      if
      (
          !adjointCycleMagnitudeGuardTriggered
       && (
              adjointCycleAdjointBestCorrector == 0
           || adjointCycleAdjointCombinedMetric < adjointCycleAdjointBestCombinedMetric
          )
      )
      {
          adjointCycleAdjointBestCombinedMetric = adjointCycleAdjointCombinedMetric;
          adjointCycleAdjointBestUNonlinearChange = adjointCycleAdjointUNonlinearChange;
          adjointCycleAdjointBestPNonlinearChange = adjointCycleAdjointPNonlinearChange;
          adjointCycleAdjointBestUaMax = adjointCycleUaMax;
          adjointCycleAdjointBestPaRange = adjointCyclePaRange;
          adjointCycleAdjointBestPaMean = adjointCycleAdjointPaMean;
          adjointCycleAdjointBestPhiaAbsSum = adjointCycleAdjointPhiaAbsSum;
          adjointCycleAdjointBestPhiaNet = adjointCycleAdjointPhiaNet;
          adjointCycleAdjointBestCorrector = adjointCycleAdjointCorrectorsUsed;
          adjointCycleBestUaInternal = UaCurrent;
          adjointCycleBestPaInternal = paCurrent;
          adjointCycleBestPhiaInternal = phia.internalField();
          adjointCycleIsBestPass = true;
      }
      if (adjointCycleIsBestPass)
      {
          adjointCycleNoImproveStreak = 0;
      }
      else if (!adjointCycleMagnitudeGuardTriggered)
      {
          adjointCycleNoImproveStreak++;
      }

      adjointCycleAdjointFlowConverged =
      (
          adjointCycleAdjointCorrectorsUsed >= 2
       && adjointCycleAdjointUNonlinearChange <= adjointCycleAdjointUNonlinearTol
       && adjointCycleAdjointPNonlinearChange <= adjointCycleAdjointPNonlinearTol
      );

      bool adjointCycleGrowthDiverged(false);
      bool adjointCycleNoImproveDiverged(false);
      bool adjointCycleNoImproveControlledStop(false);
      if
      (
          !adjointCycleMagnitudeGuardTriggered
       && !adjointCycleAdjointFlowConverged
       && divergenceCheckEnabled
       && adjointCycleAdjointCorrectorsUsed >= adjointCycleDivergenceMinCorrectors
       && adjointCycleAdjointBestCorrector > 0
       && adjointCycleAdjointCorrectorsUsed > adjointCycleAdjointBestCorrector
       && adjointCycleAdjointUNonlinearChange
          > adjointCycleDivergenceGrowthFactor
           *Foam::max(adjointCycleAdjointBestUNonlinearChange, adjointCycleAdjointUNonlinearTol)
       && adjointCycleAdjointPNonlinearChange
          > adjointCycleDivergenceGrowthFactor
           *Foam::max(adjointCycleAdjointBestPNonlinearChange, adjointCycleAdjointPNonlinearTol)
       && adjointCycleUaMax
          > adjointCycleDivergenceGrowthFactor
           *Foam::max(adjointCycleAdjointBestUaMax, scalar(1.0))
       && adjointCyclePaRange
          > adjointCycleDivergenceGrowthFactor
           *Foam::max(adjointCycleAdjointBestPaRange, scalar(1.0))
      )
      {
          adjointCycleGrowthDiverged = true;
          adjointCycleAdjointFlowDiverged = true;
          adjointCycleAdjointDivergedCorrector = adjointCycleAdjointCorrectorsUsed;
          adjointCycleAdjointStopReason = "growth-diverged";
          if (adjointCycleAdjointBestCorrector > 0)
          {
              adjointCyclePendingRollback = true;
          }
      }
      if
      (
          !adjointCycleMagnitudeGuardTriggered
       && !adjointCycleGrowthDiverged
       && !adjointCycleAdjointFlowConverged
       && stagnationCheckEnabled
       && adjointCycleAdjointBestCorrector > 0
       && adjointCycleAdjointCorrectorsUsed > adjointCycleAdjointBestCorrector
       && adjointCycleAdjointCorrectorsUsed >= noImproveMinCorrectors
       && adjointCycleNoImproveStreak >= noImprovePatience
      )
      {
          adjointCycleNoImproveDiverged = true;
          adjointCycleNoImproveControlledStop = true;
          adjointCycleAdjointStopReason =
            std::string("no-improve-") + adjointCycleCurrentBlocker;
          adjointCyclePendingRollback = true;
      }

      if (adjointCycleAdjointFlowConverged)
      {
          adjointCycleAdjointStopReason = "converged";
      }
      else if
      (
          !adjointCycleAdjointFlowDiverged
       && adjointCycleAdjointCorrectorsUsed >= adjointFlowCorrectors
      )
      {
          adjointCycleAdjointStopReason = "max-correctors";
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

          std::string adjointCyclePassState("run");
          if (adjointCycleAdjointFlowConverged)
          {
              adjointCyclePassState = "converged";
          }
          else if (adjointCycleNoImproveControlledStop)
          {
              adjointCyclePassState = "controlled-stop";
          }
          else if (adjointCycleNoImproveDiverged)
          {
              adjointCyclePassState = "stalled";
          }
          else if (adjointCycleMagnitudeGuardTriggered || adjointCycleGrowthDiverged || adjointCycleNoImproveDiverged)
          {
              adjointCyclePassState = "diverged";
          }
          else if (adjointCycleIsBestPass)
          {
              adjointCyclePassState = "best";
          }

          debugAdjointCycleLogFileStream
            << "  "
            << std::left
            << std::setw(6) << adjointCycleAdjointCorrectorsUsed
            << std::setw(12) << adjointCyclePassState
            << std::setw(14) << fmtTraceScalar(cmptMax(UaLinearSolvePerformance.initialResidual()))
            << std::setw(14) << fmtTraceScalar(cmptMax(UaLinearSolvePerformance.finalResidual()))
            << std::setw(8) << cmptMax(UaLinearSolvePerformance.nIterations())
            << std::setw(14) << fmtTraceScalar(adjointCyclePaPassInitialResidual)
            << std::setw(14) << fmtTraceScalar(adjointCyclePaPassFinalResidual)
            << std::setw(8) << adjointCyclePaPassIterations
            << std::setw(14) << fmtTraceScalar(adjointCycleAdjointUNonlinearChange)
            << std::setw(14) << fmtTraceScalar(adjointCycleAdjointPNonlinearChange)
            << std::setw(14) << fmtTraceScalar(adjointCycleUaMax)
            << std::setw(14) << fmtTraceScalar(adjointCyclePaRange)
            << std::setw(14) << fmtTraceScalar(adjointCycleAdjointPaMean)
            << std::setw(14) << fmtTraceScalar(adjointCycleAdjointPhiaAbsSum)
            << std::setw(14) << fmtTraceScalar(adjointCycleAdjointPhiaNet)
            << "\n";
      }

      if (adjointCycleAdjointFlowConverged)
      {
          break;
      }
      if (adjointCyclePendingRollback)
      {
          Ua.primitiveFieldRef() = adjointCycleBestUaInternal;
          pa.primitiveFieldRef() = adjointCycleBestPaInternal;
          phia.primitiveFieldRef() = adjointCycleBestPhiaInternal;
          Ua.correctBoundaryConditions();
          pa.correctBoundaryConditions();
          fvOptions.correct(Ua);
          const bool adjointCycleHardRollback =
              adjointCycleMagnitudeGuardTriggered
           || adjointCycleGrowthDiverged;
          if (adjointCycleHardRollback)
          {
              adjointCycleAdjointRolledBack = true;
              adjointCycleAdjointFlowControlledStop = false;
              adjointCycleAdjointRollbackCorrector = adjointCycleAdjointBestCorrector;
          }
          else if (adjointCycleNoImproveControlledStop)
          {
              adjointCycleAdjointRolledBack = false;
              adjointCycleAdjointFlowControlledStop = true;
              adjointCycleAdjointRollbackCorrector = adjointCycleAdjointBestCorrector;
          }
          else
          {
              adjointCycleAdjointRolledBack = false;
              adjointCycleAdjointFlowControlledStop = false;
              adjointCycleAdjointRollbackCorrector = 0;
          }
          if (adjointCycleAdjointBestCorrector > 0)
          {
              adjointCycleAdjointUNonlinearChange = adjointCycleAdjointBestUNonlinearChange;
              adjointCycleAdjointPNonlinearChange = adjointCycleAdjointBestPNonlinearChange;
              adjointCycleAdjointCombinedMetric = adjointCycleAdjointBestCombinedMetric;
              adjointCycleAdjointUaMax = adjointCycleAdjointBestUaMax;
              adjointCycleAdjointPaRange = adjointCycleAdjointBestPaRange;
              adjointCycleAdjointPaMean = adjointCycleAdjointBestPaMean;
              adjointCycleAdjointPhiaAbsSum = adjointCycleAdjointBestPhiaAbsSum;
              adjointCycleAdjointPhiaNet = adjointCycleAdjointBestPhiaNet;
          }
          if (adjointCycleGrowthDiverged)
          {
              adjointCycleAdjointStopReason = "growth-diverged-rollback";
          }
          else if (adjointCycleNoImproveDiverged)
          {
              adjointCycleAdjointStopReason +=
                proModeFlowAdjointMode && proModeWarmupProgressUa < scalar(1.0)
              ? "-warmup-stop"
              : "-stop";
          }
          else if (adjointCycleMagnitudeGuardTriggered)
          {
              adjointCycleAdjointStopReason = "magnitude-guard-rollback";
          }
      }
      if (adjointCycleMagnitudeGuardTriggered || adjointCycleGrowthDiverged || adjointCycleNoImproveDiverged)
      {
          break;
      }
  }

  Ua.storePrevIter();
  pa.storePrevIter();
  phia.storePrevIter();
}
