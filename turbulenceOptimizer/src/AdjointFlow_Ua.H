const bool stageCFlowAdjointMode = customTurbulenceActive && stageCFullAdjointVerification;
const bool stageDFlowAdjointMode =
    customTurbulenceActive && stageDTurbulentOptimizationMode;
const bool boundedFlowAdjointMode = stageCFlowAdjointMode || stageDFlowAdjointMode;
const bool resetAdjointState =
    (stageCFlowAdjointMode && static_cast<bool>(stageCResetAdjointState))
 || (stageDFlowAdjointMode && static_cast<bool>(stageDResetAdjointState));
const bool stagnationCheckEnabled =
    (stageCFlowAdjointMode && static_cast<bool>(stageCStagnationCheck))
 || (stageDFlowAdjointMode && static_cast<bool>(stageDNoImproveCheck));
const label noImprovePatience =
    stageDFlowAdjointMode ? stageDNoImprovePatience : stageCStagnationPatience;
const label noImproveMinCorrectors =
    stageDFlowAdjointMode ? stageDNoImproveMinCorrectors : label(1);
const bool divergenceCheckEnabled = stageCFlowAdjointMode && stageCDivergenceCheck;
const scalar flowAdjointMagnitudeLimit =
    stageDFlowAdjointMode ? stageDAdjointMaxMagnitude : stageCAdjointMaxMagnitude;
const bool useNuEffInUaEquation =
    stageDFlowAdjointMode || (stageCFlowAdjointMode && stageCUseNuEff);
const bool includeTaSourceInUaEquation =
    stageDFlowAdjointMode || (stageCFlowAdjointMode && stageCIncludeThermalMomentumSource);
const label adjointFlowCorrectors =
    stageDFlowAdjointMode
  ? stageDFlowAdjointCorrectors
  : (stageCFlowAdjointMode ? stageCAdjointCorrectors : 40);
stageCAdjointUNonlinearChange = 0.0;
stageCAdjointPNonlinearChange = 0.0;
stageCAdjointFlowConverged = !boundedFlowAdjointMode;
stageCAdjointCorrectorsUsed = 0;
stageCAdjointCombinedMetric = 0.0;
stageCAdjointBestCombinedMetric = GREAT;
stageCAdjointBestUNonlinearChange = GREAT;
stageCAdjointBestPNonlinearChange = GREAT;
stageCAdjointBestUaMax = GREAT;
stageCAdjointBestPaRange = GREAT;
stageCAdjointBestPaMean = 0.0;
stageCAdjointBestPhiaAbsSum = 0.0;
stageCAdjointBestPhiaNet = 0.0;
stageCAdjointBestCorrector = 0;
stageCAdjointUaMax = 0.0;
stageCAdjointPaRange = 0.0;
stageCAdjointPaMean = 0.0;
stageCAdjointPhiaAbsSum = 0.0;
stageCAdjointPhiaNet = 0.0;
stageCAdjointFlowDiverged = false;
stageCAdjointDivergedCorrector = 0;
stageCAdjointRolledBack = false;
stageCAdjointRollbackCorrector = 0;
stageCAdjointStopReason = boundedFlowAdjointMode ? "max-correctors" : "not-bounded";

vectorField stageCBestUaInternal;
scalarField stageCBestPaInternal;
scalarField stageCBestPhiaInternal;

if (boundedFlowAdjointMode)
{
  if (resetAdjointState)
  {
    // Default Stage-C behavior: zero the adjoint state so each verification
    // sweep is isolated. Warm-starting remains available as an A/B toggle.
    Ua = dimensionedVector("zeroStageCUa", Ua.dimensions(), vector::zero);
    pa = dimensionedScalar("zeroStageCPa", pa.dimensions(), 0.0);
    phia = dimensionedScalar("zeroStageCPhia", phia.dimensions(), 0.0);
    Ua.correctBoundaryConditions();
    pa.correctBoundaryConditions();
  }

  if (Pstream::master())
  {
      std::ofstream debugStageCLogFileStream(debugStageCLogFile, std::ios::app);
      const char* stageCTraceTag =
          stageDFlowAdjointMode
        ? "[Stage-D Flow Adjoint Corrector Trace]"
        : "[Stage-C Flow Adjoint Corrector Trace]";
      debugStageCLogFileStream
        << "\n" << stageCTraceTag << "\n"
        << "  Iteration           : " << opt << "\n"
        << "  Config              : reset="
        << (resetAdjointState ? "true" : "false")
        << " viscosity=" << (useNuEffInUaEquation ? "nuEff" : "nu")
        << " TaSource="
        << (includeTaSourceInUaEquation ? "true" : "false")
        << " legacyPaRelax="
        << (stageCUseLegacyPressureRelax ? "true" : "false")
        << " divergenceGuard="
        << (divergenceCheckEnabled ? "true" : "false")
        << " stagnationGuard="
        << (stagnationCheckEnabled ? "true" : "false")
        << " |Ua|maxLimit=" << flowAdjointMagnitudeLimit
        << " growth=" << stageCDivergenceGrowthFactor
        << " minPass="
        << (stageDFlowAdjointMode ? noImproveMinCorrectors : stageCDivergenceMinCorrectors)
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

  stageCBestUaInternal = Ua.internalField();
  stageCBestPaInternal = pa.internalField();
  stageCBestPhiaInternal = phia.internalField();
}

label stageCNoImproveStreak = 0;

for (i = 0; i < adjointFlowCorrectors; i++)
{
  vectorField UaPrevious;
  scalarField paPrevious;
  scalar stageCPaPassInitialResidual = 0.0;
  scalar stageCPaPassFinalResidual = 0.0;
  label stageCPaPassIterations = 0;
  if (boundedFlowAdjointMode)
  {
      UaPrevious = Ua.internalField();
      paPrevious = pa.internalField();
  }

  volVectorField adjointTransposeConvectiona((fvc::grad(U) & Ua));

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
      UaEqn += b1*Ta*fvc::grad(T);
  }

  UaEqn.relax();

  fvOptions.constrain(UaEqn);

  const auto UaLinearSolvePerformance = solve(UaEqn == -fvc::grad(pa));
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
          fvm::laplacian(rAtUa(), pa) == fvc::div(phiHbyAa)
      );

      paEqn.setReference(paRefCell, paRefValue);
      const auto paLinearSolvePerformance = paEqn.solve();
      recordLinearSolve
      (
          paLinearSolveStats,
          paLinearSolvePerformance.initialResidual(),
          paLinearSolvePerformance.finalResidual(),
          paLinearSolvePerformance.nIterations()
      );
      stageCPaPassInitialResidual =
        Foam::max(stageCPaPassInitialResidual, paLinearSolvePerformance.initialResidual());
      stageCPaPassFinalResidual =
        Foam::max(stageCPaPassFinalResidual, paLinearSolvePerformance.finalResidual());
      stageCPaPassIterations =
        Foam::max(stageCPaPassIterations, paLinearSolvePerformance.nIterations());

      if (simple.finalNonOrthogonalIter())
      {
          phia = phiHbyAa - paEqn.flux();
      }
  }

  if (!boundedFlowAdjointMode || stageCUseLegacyPressureRelax)
  {
      pa.relax();
  }

  Ua = HbyAa - rAtUa()*fvc::grad(pa);
  Ua.correctBoundaryConditions();
  fvOptions.correct(Ua);

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

      stageCAdjointUNonlinearChange =
        Foam::sqrt(localUaDeltaSqSum/Foam::max(localUaScaleSqSum, SMALL));
      stageCAdjointPNonlinearChange =
        Foam::sqrt(localPaDeltaSqSum/Foam::max(localPaScaleSqSum, SMALL));
      stageCAdjointCorrectorsUsed = i + 1;
      const volScalarField stageCUaMag(mag(Ua));
      const scalar stageCUaMax = gMax(stageCUaMag);
      const scalar stageCPaRange = gMax(pa) - gMin(pa);
      stageCAdjointUaMax = stageCUaMax;
      stageCAdjointPaRange = stageCPaRange;
      stageCAdjointPaMean = localPaSum/area;
      stageCAdjointPhiaAbsSum = localPhiaAbsSum;
      stageCAdjointPhiaNet = localPhiaNet;
      const scalar stageCUaTargetRatio =
        stageCAdjointUNonlinearChange/Foam::max(stageCAdjointUNonlinearTol, SMALL);
      const scalar stageCPaTargetRatio =
        stageCAdjointPNonlinearChange/Foam::max(stageCAdjointPNonlinearTol, SMALL);
      const char* stageCCurrentBlocker =
          (stageCUaTargetRatio > 1.1*stageCPaTargetRatio)
        ? "Ua"
        : (stageCPaTargetRatio > 1.1*stageCUaTargetRatio)
          ? "pa"
          : "Ua+pa";
      stageCAdjointCombinedMetric =
        Foam::max
        (
            stageCUaTargetRatio,
            stageCPaTargetRatio
        );

      bool stageCIsBestPass = false;
      bool stageCPendingRollback(false);

      if (stageCUaMax != stageCUaMax)
      {
          FatalErrorInFunction
              << "Stage-C Ua adjoint magnitude became invalid (|Ua|max="
              << stageCUaMax
              << "). This indicates verification-loop divergence."
              << exit(FatalError);
      }

      bool stageCMagnitudeGuardTriggered(false);
      if (stageCUaMax > flowAdjointMagnitudeLimit)
      {
          stageCMagnitudeGuardTriggered = true;
          stageCAdjointFlowDiverged = true;
          stageCAdjointDivergedCorrector = stageCAdjointCorrectorsUsed;
          stageCAdjointStopReason = "magnitude-guard";
          stageCPendingRollback = true;
      }

      if
      (
          !stageCMagnitudeGuardTriggered
       && (
              stageCAdjointBestCorrector == 0
           || stageCAdjointCombinedMetric < stageCAdjointBestCombinedMetric
          )
      )
      {
          stageCAdjointBestCombinedMetric = stageCAdjointCombinedMetric;
          stageCAdjointBestUNonlinearChange = stageCAdjointUNonlinearChange;
          stageCAdjointBestPNonlinearChange = stageCAdjointPNonlinearChange;
          stageCAdjointBestUaMax = stageCUaMax;
          stageCAdjointBestPaRange = stageCPaRange;
          stageCAdjointBestPaMean = stageCAdjointPaMean;
          stageCAdjointBestPhiaAbsSum = stageCAdjointPhiaAbsSum;
          stageCAdjointBestPhiaNet = stageCAdjointPhiaNet;
          stageCAdjointBestCorrector = stageCAdjointCorrectorsUsed;
          stageCBestUaInternal = UaCurrent;
          stageCBestPaInternal = paCurrent;
          stageCBestPhiaInternal = phia.internalField();
          stageCIsBestPass = true;
      }
      if (stageCIsBestPass)
      {
          stageCNoImproveStreak = 0;
      }
      else if (!stageCMagnitudeGuardTriggered)
      {
          stageCNoImproveStreak++;
      }

      stageCAdjointFlowConverged =
      (
          stageCAdjointCorrectorsUsed >= 2
       && stageCAdjointUNonlinearChange <= stageCAdjointUNonlinearTol
       && stageCAdjointPNonlinearChange <= stageCAdjointPNonlinearTol
      );

      bool stageCGrowthDiverged(false);
      bool stageCNoImproveDiverged(false);
      if
      (
          !stageCMagnitudeGuardTriggered
       && !stageCAdjointFlowConverged
       && divergenceCheckEnabled
       && stageCAdjointCorrectorsUsed >= stageCDivergenceMinCorrectors
       && stageCAdjointBestCorrector > 0
       && stageCAdjointCorrectorsUsed > stageCAdjointBestCorrector
       && stageCAdjointUNonlinearChange
          > stageCDivergenceGrowthFactor
           *Foam::max(stageCAdjointBestUNonlinearChange, stageCAdjointUNonlinearTol)
       && stageCAdjointPNonlinearChange
          > stageCDivergenceGrowthFactor
           *Foam::max(stageCAdjointBestPNonlinearChange, stageCAdjointPNonlinearTol)
       && stageCUaMax
          > stageCDivergenceGrowthFactor
           *Foam::max(stageCAdjointBestUaMax, scalar(1.0))
       && stageCPaRange
          > stageCDivergenceGrowthFactor
           *Foam::max(stageCAdjointBestPaRange, scalar(1.0))
      )
      {
          stageCGrowthDiverged = true;
          stageCAdjointFlowDiverged = true;
          stageCAdjointDivergedCorrector = stageCAdjointCorrectorsUsed;
          stageCAdjointStopReason = "growth-diverged";
          if (stageCAdjointBestCorrector > 0)
          {
              stageCPendingRollback = true;
          }
      }
      if
      (
          !stageCMagnitudeGuardTriggered
       && !stageCGrowthDiverged
       && !stageCAdjointFlowConverged
       && stagnationCheckEnabled
       && stageCAdjointBestCorrector > 0
       && stageCAdjointCorrectorsUsed > stageCAdjointBestCorrector
       && stageCAdjointCorrectorsUsed >= noImproveMinCorrectors
       && stageCNoImproveStreak >= noImprovePatience
      )
      {
          stageCNoImproveDiverged = true;
          stageCAdjointFlowDiverged = true;
          stageCAdjointDivergedCorrector = stageCAdjointCorrectorsUsed;
          stageCAdjointStopReason =
            std::string("no-improve-") + stageCCurrentBlocker;
          stageCPendingRollback = true;
      }

      if (stageCAdjointFlowConverged)
      {
          stageCAdjointStopReason = "converged";
      }
      else if
      (
          !stageCAdjointFlowDiverged
       && stageCAdjointCorrectorsUsed >= adjointFlowCorrectors
      )
      {
          stageCAdjointStopReason = "max-correctors";
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

          std::string stageCPassState("run");
          if (stageCAdjointFlowConverged)
          {
              stageCPassState = "converged";
          }
          else if (stageCNoImproveDiverged)
          {
              stageCPassState = "stalled";
          }
          else if (stageCMagnitudeGuardTriggered || stageCGrowthDiverged || stageCNoImproveDiverged)
          {
              stageCPassState = "diverged";
          }
          else if (stageCIsBestPass)
          {
              stageCPassState = "best";
          }

          debugStageCLogFileStream
            << "  "
            << std::left
            << std::setw(6) << stageCAdjointCorrectorsUsed
            << std::setw(12) << stageCPassState
            << std::setw(14) << fmtTraceScalar(cmptMax(UaLinearSolvePerformance.initialResidual()))
            << std::setw(14) << fmtTraceScalar(cmptMax(UaLinearSolvePerformance.finalResidual()))
            << std::setw(8) << cmptMax(UaLinearSolvePerformance.nIterations())
            << std::setw(14) << fmtTraceScalar(stageCPaPassInitialResidual)
            << std::setw(14) << fmtTraceScalar(stageCPaPassFinalResidual)
            << std::setw(8) << stageCPaPassIterations
            << std::setw(14) << fmtTraceScalar(stageCAdjointUNonlinearChange)
            << std::setw(14) << fmtTraceScalar(stageCAdjointPNonlinearChange)
            << std::setw(14) << fmtTraceScalar(stageCUaMax)
            << std::setw(14) << fmtTraceScalar(stageCPaRange)
            << std::setw(14) << fmtTraceScalar(stageCAdjointPaMean)
            << std::setw(14) << fmtTraceScalar(stageCAdjointPhiaAbsSum)
            << std::setw(14) << fmtTraceScalar(stageCAdjointPhiaNet)
            << "\n";
      }

      if (stageCAdjointFlowConverged)
      {
          break;
      }
      if (stageCPendingRollback)
      {
          Ua.primitiveFieldRef() = stageCBestUaInternal;
          pa.primitiveFieldRef() = stageCBestPaInternal;
          phia.primitiveFieldRef() = stageCBestPhiaInternal;
          Ua.correctBoundaryConditions();
          pa.correctBoundaryConditions();
          fvOptions.correct(Ua);

          stageCAdjointRolledBack = true;
          stageCAdjointRollbackCorrector = stageCAdjointBestCorrector;
          if (stageCAdjointBestCorrector > 0)
          {
              stageCAdjointUNonlinearChange = stageCAdjointBestUNonlinearChange;
              stageCAdjointPNonlinearChange = stageCAdjointBestPNonlinearChange;
              stageCAdjointCombinedMetric = stageCAdjointBestCombinedMetric;
              stageCAdjointUaMax = stageCAdjointBestUaMax;
              stageCAdjointPaRange = stageCAdjointBestPaRange;
              stageCAdjointPaMean = stageCAdjointBestPaMean;
              stageCAdjointPhiaAbsSum = stageCAdjointBestPhiaAbsSum;
              stageCAdjointPhiaNet = stageCAdjointBestPhiaNet;
          }
          if (stageCGrowthDiverged)
          {
              stageCAdjointStopReason = "growth-diverged-rollback";
          }
          else if (stageCNoImproveDiverged)
          {
              stageCAdjointStopReason += "-rollback";
          }
          else if (stageCMagnitudeGuardTriggered)
          {
              stageCAdjointStopReason = "magnitude-guard-rollback";
          }
      }
      if (stageCMagnitudeGuardTriggered || stageCGrowthDiverged || stageCNoImproveDiverged)
      {
          break;
      }
  }

  Ua.storePrevIter();
  pa.storePrevIter();
  phia.storePrevIter();
}
