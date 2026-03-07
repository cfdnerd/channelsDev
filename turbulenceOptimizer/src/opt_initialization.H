Info<< "Initialization\n" << endl;
const char* optimizationHistoryFile = "optimization.hst";
const char* maxPD0SugFile = "maxPD0Sug.log";
const char* debugOptimizerJsonFile = "debugOptimizer.jsonl";
const char* debugOptimizerLogFile = "debugOptimizer.log";
const char* debugAdjointCycleLogFile = debugOptimizerLogFile;
const char* legacyDebugAdjointCycleLogFile = "debugAdjointCycle.log";
const char* solverConvergenceFile = "solverConvergences.log";
const char* optimizerDirectionEventFile = "optimizerDirectionEvents.log";
const char* firstFailureSnapshotFile = "firstFailureSnapshot.log";
if (Pstream::master())
{
  std::remove(optimizationHistoryFile);
  std::remove(maxPD0SugFile);
  std::remove(debugOptimizerJsonFile);
  std::remove(debugOptimizerLogFile);
  std::remove(legacyDebugAdjointCycleLogFile);
  std::remove(solverConvergenceFile);
  std::remove(optimizerDirectionEventFile);
  std::remove(firstFailureSnapshotFile);
  std::ofstream historyFile(optimizationHistoryFile, std::ios::out);
  historyFile
    << "# MTO_ThermalFluid optimization history\n"
    << "# Cols: Iter PDActLog PDCon PD0Sug V MeanT MaxT\n"
    << "# PDActLog=PDCon*PowerDiss0, PDCon=PowerDiss used by MMA\n"
    << "# PD0Sug=PowerDiss0*(PDCon/PowerDissMax) when PDCon>PowerDissRelax; else PowerDiss0\n"
    << std::left
    << std::setw(8) << "Iter"
    << std::setw(16) << "PDActLog"
    << std::setw(16) << "PDCon"
    << std::setw(16) << "PD0Sug"
    << std::setw(16) << "V"
    << std::setw(16) << "MeanT"
    << std::setw(16) << "MaxT"
    << '\n';

  std::ofstream maxPD0SugLogFile(maxPD0SugFile, std::ios::out);
  maxPD0SugLogFile
    << "# MTO_ThermalFluid cumulative PD0Sug maxima\n"
    << "# Cols: Iter PD0SugMax_1_to_Iter\n"
    << std::left
    << std::setw(8) << "Iter"
    << std::setw(20) << "PD0SugMax_1_to_Iter"
    << '\n';

  std::ofstream debugLogFile(debugOptimizerLogFile, std::ios::out);
  debugLogFile
    << "# MTO_ThermalFluid optimizer debug log\n"
    << "# Structured per-iteration diagnostics for optimizer development\n"
    << "# Format: sectioned human-readable report, one block per optimization iteration\n"
    << "# Adjoint-Cycle/Optimization-Cycle adjoint traces are consolidated in this file under dedicated sections\n";

  if (static_cast<bool>(optRogueDirectionSwitchEnable))
  {
    std::ofstream optimizerDirectionEventLogFile(optimizerDirectionEventFile, std::ios::out);
    optimizerDirectionEventLogFile
      << "# MTO_ThermalFluid optimizer direction-search event log\n"
      << "# Event lines for rogue-watch activation, grace countdown, and forced direction switches\n"
      << "# Format: Iter=<i> Event=<name> key=value ...\n";
  }

  std::ofstream solverConvergenceLogFile(solverConvergenceFile, std::ios::out);
  solverConvergenceLogFile
    << "# MTO_ThermalFluid solver convergence log\n"
    << "# Linear-solver statistics grouped by optimization iteration\n"
    << "# Config: endTime=" << runTime.endTime().value()
    << " heatSource=" << heatSource.value()
    << '\n';

  if (static_cast<bool>(optPersistFirstFailureSnapshot))
  {
    std::ofstream firstFailureSnapshotStream(firstFailureSnapshotFile, std::ios::out);
    firstFailureSnapshotStream
      << "# MTO_ThermalFluid first-failure snapshot bundle\n"
      << "# Stores the first optimizer health-failure event for reproducible debugging\n"
      << "# Fields are written at failure time if optPersistFirstFailureSnapshot=true\n";
  }
}

if (optimizationCycleTurbulentOptimizationMode && static_cast<bool>(optFormulationExecutionGateEnable))
{
  if
  (
      static_cast<bool>(optRequireTransportedClosureAdjoints)
   && !static_cast<bool>(optimizationCycleTransportedClosureAdjoints)
  )
  {
    FatalErrorInFunction
      << "Formulation-first execution gate requires "
      << "customTurbulenceProperties.optimizationCycleTransportedClosureAdjoints=true "
      << "for Optimization-Cycle runs."
      << exit(FatalError);
  }
  if
  (
      static_cast<bool>(optRequireTransportedClosureAdjoints)
   && optimizationCycleClosureAdjointFormulation != word("transportedResidualDraft")
   && optimizationCycleClosureAdjointFormulation != word("transportedCrossCoupledV3")
  )
  {
    FatalErrorInFunction
      << "Formulation-first execution gate requires "
      << "customTurbulenceProperties.optimizationCycleClosureAdjointFormulation "
      << "to be transportedResidualDraft or transportedCrossCoupledV3. Found "
      << optimizationCycleClosureAdjointFormulation
      << exit(FatalError);
  }
  if
  (
      static_cast<bool>(optRequireDualConsistentClosureBC)
   && !static_cast<bool>(optimizationCycleClosureAdjointDualConsistentBC)
  )
  {
    FatalErrorInFunction
      << "Formulation-first execution gate requires "
      << "customTurbulenceProperties.optimizationCycleClosureAdjointDualConsistentBC=true "
      << "for Optimization-Cycle runs."
      << exit(FatalError);
  }
  if (static_cast<bool>(smoothTurbulenceClipping))
  {
    FatalErrorInFunction
      << "Rule violation: smooth turbulence clipping must remain disabled. "
      << "Use first-principles/formulation fixes and fail-fast diagnostics."
      << exit(FatalError);
  }
  if
  (
      static_cast<bool>(optRequireGradientChecks)
   && !static_cast<bool>(optGradientCheckEnable)
  )
  {
    FatalErrorInFunction
      << "Formulation-first execution gate requires "
      << "tuneOptParameters.convergenceControl.optGradientCheckEnable=true "
      << "for Optimization-Cycle runs."
      << exit(FatalError);
  }
  if
  (
      static_cast<bool>(optRequireClosureDerivativeChecks)
   && !static_cast<bool>(optClosureDerivativeFDStageCheckEnable)
  )
  {
    FatalErrorInFunction
      << "Formulation-first execution gate requires "
      << "tuneOptParameters.convergenceControl.optClosureDerivativeFDStageCheckEnable=true "
      << "for Optimization-Cycle runs."
      << exit(FatalError);
  }
  if (!static_cast<bool>(optSensitivityQualityGateEnable))
  {
    FatalErrorInFunction
      << "Formulation-first execution gate requires "
      << "tuneOptParameters.convergenceControl.optSensitivityQualityGateEnable=true "
      << "for Optimization-Cycle runs."
      << exit(FatalError);
  }
  const dictionary& solverControls = mesh.solutionDict().subDict("solvers");
  if (!solverControls.found("ka"))
  {
    FatalErrorInFunction
      << "Optimization-Cycle transported closure adjoints require a dedicated 'ka' block "
      << "in fvSolution.solvers."
      << exit(FatalError);
  }
  if (!solverControls.found("omegaa"))
  {
    FatalErrorInFunction
      << "Optimization-Cycle transported closure adjoints require a dedicated 'omegaa' "
      << "block in fvSolution.solvers."
      << exit(FatalError);
  }
  if
  (
      static_cast<bool>(optDualConsistentStateMapEnforce)
   && (
         !static_cast<bool>(optimizationCycleTransportedClosureAdjoints)
      || !static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
      || optimizationCycleClosureAdjointFormulation == word("proxyNormalized")
      )
  )
  {
    FatalErrorInFunction
      << "Dual-consistent state-map enforcement requires "
      << "transported closure adjoints, bounded closure state mapping, "
      << "and non-proxy closure-adjoint formulation in Optimization-Cycle."
      << exit(FatalError);
  }
}

if (optRunMode != word("diagnostic") && optRunMode != word("production"))
{
  FatalErrorInFunction
    << "tuneOptParameters.convergenceControl.optRunMode must be "
    << "diagnostic or production. Found " << optRunMode
    << exit(FatalError);
}

struct LinearSolveStats
{
  label callCount;
  scalar initialResidualSum;
  scalar initialResidualMax;
  scalar finalResidualSum;
  scalar finalResidualMax;
  label iterationSum;
  label iterationMax;
  scalar lastInitialResidual;
  scalar lastFinalResidual;
  label lastIterations;

  LinearSolveStats()
  :
    callCount(0),
    initialResidualSum(0),
    initialResidualMax(0),
    finalResidualSum(0),
    finalResidualMax(0),
    iterationSum(0),
    iterationMax(0),
    lastInitialResidual(0),
    lastFinalResidual(0),
    lastIterations(0)
  {}
};

auto resetLinearSolveStats =
[](LinearSolveStats& stats)
{
  stats.callCount = 0;
  stats.initialResidualSum = 0;
  stats.initialResidualMax = 0;
  stats.finalResidualSum = 0;
  stats.finalResidualMax = 0;
  stats.iterationSum = 0;
  stats.iterationMax = 0;
  stats.lastInitialResidual = 0;
  stats.lastFinalResidual = 0;
  stats.lastIterations = 0;
};

auto recordLinearSolve =
[](LinearSolveStats& stats, const scalar initialResidual, const scalar finalResidual, const label iterations)
{
  stats.callCount++;
  stats.initialResidualSum += initialResidual;
  stats.initialResidualMax = Foam::max(stats.initialResidualMax, initialResidual);
  stats.finalResidualSum += finalResidual;
  stats.finalResidualMax = Foam::max(stats.finalResidualMax, finalResidual);
  stats.iterationSum += iterations;
  stats.iterationMax = Foam::max(stats.iterationMax, iterations);
  stats.lastInitialResidual = initialResidual;
  stats.lastFinalResidual = finalResidual;
  stats.lastIterations = iterations;
};

const scalar linearSolveObsResidual = 1e-6;
const scalar linearSolveWarnResidual = 1e-4;

auto classifyLinearSolveStatus =
[&](const LinearSolveStats& stats) -> const char*
{
  if (stats.callCount == 0)
  {
    return "SKIP";
  }
  if (stats.finalResidualMax > linearSolveWarnResidual)
  {
    return "WARN";
  }
  if (stats.finalResidualMax > linearSolveObsResidual)
  {
    return "OBS";
  }
  return "OK";
};

auto linearSolveSeverity =
[](const char* status) -> label
{
  if (std::strcmp(status, "WARN") == 0)
  {
    return 2;
  }
  if (std::strcmp(status, "OBS") == 0)
  {
    return 1;
  }
  return 0;
};

auto linearSolveSeverityFromStats =
[&](const LinearSolveStats& stats) -> label
{
  return linearSolveSeverity(classifyLinearSolveStatus(stats));
};

auto computeRelativeFieldChange =
[](const volScalarField& field, const scalarField& previousValues, const scalar scaleFloor) -> scalar
{
  const scalarField& currentValues = field.internalField();
  scalar localDeltaSqSum = 0.0;
  scalar localScaleSqSum = 0.0;

  forAll(currentValues, celli)
  {
    const scalar referenceScale =
      Foam::max
      (
        Foam::max(Foam::mag(previousValues[celli]), Foam::mag(currentValues[celli])),
        scaleFloor
      );
    const scalar deltaValue = currentValues[celli] - previousValues[celli];

    localDeltaSqSum += sqr(deltaValue);
    localScaleSqSum += sqr(referenceScale);
  }

  reduce(localDeltaSqSum, sumOp<scalar>());
  reduce(localScaleSqSum, sumOp<scalar>());

  return Foam::sqrt(localDeltaSqSum/Foam::max(localScaleSqSum, SMALL));
};

LinearSolveStats ULinearSolveStats;
LinearSolveStats pLinearSolveStats;
LinearSolveStats kLinearSolveStats;
LinearSolveStats omegaLinearSolveStats;
LinearSolveStats kaLinearSolveStats;
LinearSolveStats omegaaLinearSolveStats;
LinearSolveStats TLinearSolveStats;
LinearSolveStats TaLinearSolveStats;
LinearSolveStats TbLinearSolveStats;
LinearSolveStats UbLinearSolveStats;
LinearSolveStats pbLinearSolveStats;
LinearSolveStats UaLinearSolveStats;
LinearSolveStats paLinearSolveStats;
LinearSolveStats xpFilterLinearSolveStats;
LinearSolveStats fsensFilterLinearSolveStats;
LinearSolveStats gsensPowerFilterLinearSolveStats;
LinearSolveStats gsensVolFilterLinearSolveStats;

int i=0,opt=0;
double del=0;
scalar n(0),N(0),eta0(0),eta5(0),eta1(0),y0(0),y5(0);
n=x.size();//n:length of x in the i-th core
N=x.size();//N:design freedom
reduce(N, sumOp<scalar>()); 

scalar V(0),area=gSum(mesh.V()),MeanT(0),PowerDiss(0),PowerDissActual(0),PDActLog(0),PD0Sug(0),PD0SugMaxCum(-GREAT),PowerDissActualOverBase(0),maxT(0),sumT(0);
scalar thermalExcessObjective(0.0);
scalar thermalTmaxKSConstraint(0.0);
scalar thermalTminKSConstraint(0.0);
scalar previousObjectiveValue(0);
bool hasPreviousObjectiveValue(false);
scalar previousPowerConstraintMargin(0);
bool hasPreviousPowerConstraintMargin(false);
scalar previousPowerConstraintFinalMargin(0);
bool hasPreviousPowerConstraintFinalMargin(false);
scalar previousVolumeConstraintMargin(0);
bool hasPreviousVolumeConstraintMargin(false);
scalar previousThermalTmaxConstraintMargin(0);
bool hasPreviousThermalTmaxConstraintMargin(false);
scalar previousThermalTminConstraintMargin(0);
bool hasPreviousThermalTminConstraintMargin(false);
scalar previousMaxConstraintViolation(0);
bool hasPreviousMaxConstraintViolation(false);
scalar previousGuardConstraintViolation(0);
bool hasPreviousGuardConstraintViolation(false);
scalar guardConstraintViolationAtCountStart(0);
scalar guardConstraintViolationMinSinceCountStart(0);
bool hasGuardConstraintViolationAtCountStart(false);
label guardConstraintViolationCountStartIter(-1);
label powerConstraintWorsenStreak(0);
label maxViolationWorsenStreak(0);
label objectiveWorsenStreak(0);
label clipGuardStreak(0);
label clipGuardNoProgressStreak(0);
bool stopOptimization(false);
label optimizerConvergedStreak(0);
word optimizerStopReason("running");
bool constraintRecoveryActive(false);
scalar continuationRampScaleActive(1.0);
scalar delRampStepActive(0.2);
scalar optimizerMoveLimitActive(movlim);
label continuationIteration(0);
label continuationScheduleIterationActive(1);
label previousContinuationScheduleIterationActive(1);
bool previousConstraintWorsenScheduleReady(false);
bool hasPreviousConstraintWorsenScheduleReady(false);
bool continuationFrozenByReject(false);
bool continuationFrozenByLowTrust(false);
bool continuationFrozenByClipHold(false);
bool continuationAdvanceOverrideOnReject(false);
label continuationAdvanceCatchUpSteps(0);
label continuationDeadlockRecoveryCount(0);
bool meritStepRejected(false);
label meritRejectStreak(0);
label meritRejectStopEligibleStreak(0);
scalar meritTrustMoveScale(1.0);
bool meritRecoveryRequest(false);
bool rogueDirectionWatchActive(false);
label rogueDirectionWatchGraceRemaining(0);
label rogueDirectionWatchStartIteration(-1);
label rogueDirectionSwitchCount(0);
label rogueDirectionLastSwitchIteration(-1000000);
scalar previousAcceptedMerit(0);
bool hasPreviousAcceptedMerit(false);
label acceptedMeritScheduleIteration(-1);
label meritScheduleRebaseCount(0);
scalar alPenalty(Foam::max(optALInitialPenalty, scalar(1e-8)));
scalar alLambdaPower(0.0);
scalar alLambdaVolume(0.0);
std::vector<scalar> filterConstraintArchive;
std::vector<scalar> filterObjectiveArchive;
label filterAcceptedArchiveCount(0);
label gradientCheckFailStreak(0);
scalar closureDerivativeFDStageCheckMeanRatio(0);
scalar closureDerivativeFDStageCheckMaxRatio(0);
scalar closureDerivativeFDStageCheckMaxRatioK(0);
scalar closureDerivativeFDStageCheckMaxRatioOmega(0);
scalar closureDerivativeFDStageCheckMaxStoredMismatch(0);
label closureDerivativeFDStageCheckSamplesUsed(0);
label closureDerivativeFDStageCheckFailStreak(0);
scalar previousMaxFluidClipFraction(0.0);
scalar previousMaxClipFraction(0.0);
scalar previousGuardConstraintViolationDelta(0.0);
bool hasPreviousClipDiagnostics(false);
bool previousHealthGateFailed(false);
scalar previousThermalEnvelopeViolation(0.0);
scalar previousMaxLinearFinalResidual(0.0);
bool hasPreviousMaxLinearFinalResidual(false);
scalar sensitivityQualityScore(1.0);
bool sensitivityQualityWarn(false);
bool sensitivityQualityFail(false);
bool sensitivityQualityHardFail(false);
bool sensitivityQualityRequireConservativeStep(false);
bool sensitivityQualityReuseAcceptedGradient(false);
label sensitivityQualityFailStreak(0);
bool hasAcceptedGradientSnapshot(false);
bool hasAcceptedStateSnapshot(false);
bool firstFailureSnapshotWritten(false);
scalar primalThermalUnderRelaxActive(optThermalUnderRelaxMax);
scalar primalTurbulenceUnderRelaxActive(optTurbulenceUnderRelaxMax);
scalar primalPseudoTimeCoeffActive(0.0);
scalar thermalSourceRampActive(1.0);
scalar adaptiveMmaAsymIncActive(mma_inc);
scalar adaptiveMmaAsymDecActive(mma_dec);
scalar adaptiveMmaBaseMoveActive(movlim);
scalar turbulenceKNonlinearChange(0);
scalar turbulenceOmegaNonlinearChange(0);
scalar turbulenceNutNonlinearChange(0);
scalar turbulenceInletKActual(0);
scalar turbulenceInletOmegaActual(0);
scalar turbulenceInletKRelativeError(0);
scalar turbulenceInletOmegaRelativeError(0);
bool turbulenceNonlinearConverged(false);
bool turbulenceInletConsistent(true);
label turbulenceCorrectorsUsed(0);
const char* turbulenceInletConsistencyStatus("SKIP");
scalar adjointCycleAdjointUNonlinearChange(0);
scalar adjointCycleAdjointPNonlinearChange(0);
bool adjointCycleAdjointFlowConverged(false);
label adjointCycleAdjointCorrectorsUsed(0);
scalar adjointCycleAdjointCombinedMetric(0);
scalar adjointCycleAdjointBestCombinedMetric(GREAT);
scalar adjointCycleAdjointBestUNonlinearChange(GREAT);
scalar adjointCycleAdjointBestPNonlinearChange(GREAT);
scalar adjointCycleAdjointBestUaMax(GREAT);
scalar adjointCycleAdjointBestPaRange(GREAT);
scalar adjointCycleAdjointBestPaMean(0);
scalar adjointCycleAdjointBestPhiaAbsSum(0);
scalar adjointCycleAdjointBestPhiaNet(0);
label adjointCycleAdjointBestCorrector(0);
scalar adjointCycleAdjointUaMax(0);
scalar adjointCycleAdjointPaRange(0);
scalar adjointCycleAdjointPaMean(0);
scalar adjointCycleAdjointPhiaAbsSum(0);
scalar adjointCycleAdjointPhiaNet(0);
bool adjointCycleAdjointFlowDiverged(false);
label adjointCycleAdjointDivergedCorrector(0);
bool adjointCycleAdjointRolledBack(false);
label adjointCycleAdjointRollbackCorrector(0);
word adjointCycleAdjointStopReason("not-run");
double len=std::pow(area/N, 1.0/3);
b.value()=1.0/(filterR*len/3.464)/(filterR*len/3.464);
//***********************************//MMA
const label mmaPowerConstraintIndex = 0;
const label mmaVolumeConstraintIndex = 1;
const bool mmaUseThermalEnvelopeConstraints =
  static_cast<bool>(optThermalEnvelopeConstraintEnable)
  && static_cast<bool>(optThermalEnvelopeMMAConstraintEnable);
const label mmaThermalTmaxConstraintIndex = 2;
const label mmaThermalTminConstraintIndex = 3;
int m = mmaUseThermalEnvelopeConstraints ? 4 : 2;//number of constraints
MMA mma(n,m);
mma.asyminit = mma_init;
mma.asymdec = mma_dec;
mma.asyminc = mma_inc;
mma.raa0 = raa0;
mma.movelimit = movlim;
std::vector<double> xmma(n);
std::vector<double> dfdx(n);
std::vector<double> gx(m);
std::vector<std::vector<double>> dgdx(m);
std::vector<double> acceptedDfdx(n, 0.0);
std::vector<std::vector<double>> acceptedDgdx(m);
for (int ci = 0; ci < m; ci++)
{
  acceptedDgdx[ci].resize(n, 0.0);
}
for (int i = 0; i < m; i++)
{
  dgdx[i].resize(n);
}
