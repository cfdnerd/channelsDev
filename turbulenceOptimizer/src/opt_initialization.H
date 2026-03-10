Info<< "Initialization\n" << endl;
const char* optimizationHistoryFile = "optimization.hst";
const char* healthMonitorFile = "healthMonitor.hst";
const char* suggestedPowerDiss0File = "suggestedPowerDiss0.log";
const char* debugOptimizerJsonFile = "debugOptimizer.jsonl";
const char* debugOptimizerLogFile = "debugOptimizer.log";
const char* debugAdjointCycleLogFile = debugOptimizerLogFile;
const char* legacyDebugAdjointCycleLogFile = "debugAdjointCycle.log";
const char* solverConvergenceFile = "solverConvergences.log";
const char* firstFailureSnapshotFile = "firstFailureSnapshot.log";
if (Pstream::master())
{
  std::remove(optimizationHistoryFile);
  std::remove(healthMonitorFile);
  std::remove(suggestedPowerDiss0File);
  std::remove(debugOptimizerJsonFile);
  std::remove(debugOptimizerLogFile);
  std::remove(legacyDebugAdjointCycleLogFile);
  std::remove(solverConvergenceFile);
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

  std::ofstream healthFile(healthMonitorFile, std::ios::out);
  healthFile
    << "# MTO_ThermalFluid optimizer health monitor\n"
    << "# Top-10 prioritized health metrics (P1..P10) for early stability tracking\n"
    << "# HealthScore=100-PenTotal, clamped to [0,100]\n"
    << "# PenTotal=RollPen+constPen+LinPen+gradPen+TrustPen\n"
    << "# RollPen: rollback/divergence/magnitude-guard penalties (cap 45)\n"
    << "# LinPen: 6*log10(1 + MaxLinRes/linearSolveWarnResidual), +5 if residual growth violated, cap 20\n"
    << "# constPen: 5*log10(1 + GuardCon/optConstraintTol), cap 25\n"
    << "# gradPen: gradient-consistency penalty for opt>10 (cap 10)\n"
    << "# TrustPen: trust-ratio deviation penalty when trust ratio is defined (cap 10)\n"
    << "# FlowRB/HeatRB: binary rollback indicators (0 or 1)\n"
    << "# Priority: P1=HealthScore P2=PenTotal P3=RollPen P4=LinPen P5=constPen P6=gradPen P7=TrustPen P8=FlowRB P9=HeatRB\n"
    << "# Cols: Iter HealthScore PenTotal RollPen LinPen constPen gradPen TrustPen FlowRB HeatRB\n"
    << std::left
    << std::setw(8) << "Iter"
    << std::setw(14) << "HealthScore"
    << std::setw(14) << "PenTotal"
    << std::setw(14) << "RollPen"
    << std::setw(14) << "LinPen"
    << std::setw(14) << "constPen"
    << std::setw(14) << "gradPen"
    << std::setw(14) << "TrustPen"
    << std::setw(10) << "FlowRB"
    << std::setw(10) << "HeatRB"
    << '\n';

  std::ofstream debugLogFile(debugOptimizerLogFile, std::ios::out);
  debugLogFile
    << "# MTO_ThermalFluid optimizer debug log\n"
    << "# Structured per-iteration diagnostics for optimizer development\n"
    << "# Format: sectioned human-readable report, one block per optimization iteration\n"
    << "# Adjoint-Cycle/Optimization-Cycle adjoint traces are consolidated in this file under dedicated sections\n";

  std::ofstream solverConvergenceLogFile(solverConvergenceFile, std::ios::out);
  solverConvergenceLogFile
    << "# MTO_ThermalFluid solver convergence log\n"
    << "# Linear-solver statistics grouped by optimization iteration\n"
    << "# Config: endTime=" << runTime.endTime().value()
    << " heatSource=" << heatSource.value()
    << '\n';

  if (static_cast<bool>(optFirstFailureSnapshotEnable))
  {
    std::ofstream firstFailureSnapshotStream(firstFailureSnapshotFile, std::ios::out);
    firstFailureSnapshotStream
      << "# MTO_ThermalFluid first-failure snapshot bundle\n"
      << "# Stores the first optimizer health-failure event for reproducible debugging\n"
      << "# Fields are written at failure time if optFirstFailureSnapshotEnable=true\n";
  }
}

if (optimizationCycleTurbulentOptimizationMode)
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

auto auditDivisionDenominator =
[&](const volScalarField& denominator, const char* context, const scalar minAllowed)
{
  if (!static_cast<bool>(optDivisionAuditEnable))
  {
    return;
  }

  const scalar denominatorMin = gMin(denominator);
  const scalar denominatorMax = gMax(denominator);

  label localInvalidCount = 0;
  label localFirstInvalidCell = -1;
  label localFirstInvalidPatch = -1;
  label localFirstInvalidFace = -1;
  scalar localFirstInvalidValue = 0.0;

  const scalarField& denominatorInternal = denominator.internalField();
  forAll(denominatorInternal, celli)
  {
    const scalar denominatorValue = denominatorInternal[celli];
    const bool invalidDenominator =
         !std::isfinite(static_cast<double>(denominatorValue))
      || denominatorValue <= minAllowed;
    if (!invalidDenominator)
    {
      continue;
    }

    localInvalidCount++;
    if (localFirstInvalidCell < 0)
    {
      localFirstInvalidCell = celli;
      localFirstInvalidValue = denominatorValue;
    }
  }

  const auto& denominatorBoundary = denominator.boundaryField();
  forAll(denominatorBoundary, patchi)
  {
    const fvPatchScalarField& denominatorPatch = denominatorBoundary[patchi];
    forAll(denominatorPatch, facei)
    {
      const scalar denominatorValue = denominatorPatch[facei];
      const bool invalidDenominator =
           !std::isfinite(static_cast<double>(denominatorValue))
        || denominatorValue <= minAllowed;
      if (!invalidDenominator)
      {
        continue;
      }

      localInvalidCount++;
      if (localFirstInvalidPatch < 0)
      {
        localFirstInvalidPatch = patchi;
        localFirstInvalidFace = facei;
        localFirstInvalidValue = denominatorValue;
      }
    }
  }

  label globalInvalidCount = localInvalidCount;
  reduce(globalInvalidCount, sumOp<label>());

  const bool denominatorRangeValid =
       std::isfinite(static_cast<double>(denominatorMin))
    && std::isfinite(static_cast<double>(denominatorMax))
    && denominatorMin > minAllowed;

  if (denominatorRangeValid && globalInvalidCount == 0)
  {
    return;
  }

  FatalErrorInFunction
    << "Division denominator audit failed."
    << " Context=" << context
    << " minAllowed=" << minAllowed
    << " globalMin=" << denominatorMin
    << " globalMax=" << denominatorMax
    << " globalInvalidCount=" << globalInvalidCount
    << " localInvalidCount=" << localInvalidCount
    << " localFirstInvalidCell=" << localFirstInvalidCell
    << " localFirstInvalidPatch=" << localFirstInvalidPatch
    << " localFirstInvalidFace=" << localFirstInvalidFace
    << " localFirstInvalidValue=" << localFirstInvalidValue
    << " proc=" << Pstream::myProcNo()
    << exit(FatalError);
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
scalar PD0SugSumCum(0);
label PD0SugCountCum(0);
bool PD0SugLogInitialized(false);
scalar thermalExcessObjective(0.0);
scalar thermalTmaxKSConstraint(0.0);
scalar thermalTminKSConstraint(0.0);
scalar grayPenaltyObjective(0.0);
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
bool stopOptimization(false);
label optimizerConvergedStreak(0);
word optimizerStopReason("running");
scalar optimizerMoveLimitActive(movlim);
scalar previousTrustRatio(1.0);
bool hasPreviousTrustRatio(false);
scalar previousConstraintPenaltyForGovernor(0.0);
scalar previousLinearPenaltyForGovernor(0.0);
bool previousTrustRatioDefinedForGovernor(false);
bool hasPreviousGovernorMetrics(false);
scalar turbulenceSuppressionRampScaleActive(1.0);
bool continuationGovernorStressActive(false);
scalar projectionBetaScaleActive(optProjectionBetaScaleStart);
scalar projectionDelStepActive(optProjectionDelStepStart);
scalar filterRadiusScaleActive(1.0);
scalar materialInterpolationExponentActive(optMaterialInterpolationExponentStart);
scalar grayPenaltyWeightActive(optGrayPenaltyWeightStart);
label optimizerTinyStepStreak(0);
bool optimizerMMAStagnationRecoveryActive(false);
bool optimizerResetMMAAsymptotesNextIter(false);
scalar optimizerGradientSanityScaleActive(1.0);
label continuationIteration(1);
label continuationScheduleIterationActive(1);
scalar previousMaxLinearFinalResidual(0.0);
bool hasPreviousMaxLinearFinalResidual(false);
bool firstFailureSnapshotWritten(false);
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
bool adjointCycleAdjointFlowControlledStop(false);
label adjointCycleAdjointRollbackCorrector(0);
word adjointCycleAdjointStopReason("not-run");
bool optimizationCycleHeatAdjointRolledBack(false);
bool optimizationCycleHeatAdjointControlledStop(false);
label optimizationCycleHeatAdjointRollbackCorrector(0);
word optimizationCycleHeatAdjointStopReason("not-run");
label proModeAdjointRollbackFreeStreak(0);
const label proModeRollbackDiagnosticWindow(10);
label proModeFlowRollbackTotal(0);
label proModeHeatRollbackTotal(0);
label proModeFlowMagnitudeGuardTotal(0);
label proModeHeatMagnitudeGuardTotal(0);
label proModeFlowNoImproveTotal(0);
label proModeHeatNoImproveTotal(0);
label proModeFlowGrowthTotal(0);
label proModeFlowRollbackWindowCount(0);
label proModeHeatRollbackWindowCount(0);
label proModeFlowMagnitudeGuardWindowCount(0);
label proModeHeatMagnitudeGuardWindowCount(0);
label proModeFlowNoImproveWindowCount(0);
label proModeHeatNoImproveWindowCount(0);
label proModeFlowGrowthWindowCount(0);
bool proModeForceAdjointResetNextIter(false);
volScalarField proModeLaggedKa
(
  IOobject
  (
    "proModeLaggedKa",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  ka
);
volScalarField proModeLaggedOmegaa
(
  IOobject
  (
    "proModeLaggedOmegaa",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  omegaa
);
bool proModeLaggedClosureInitialized(false);
double len=std::pow(area/N, 1.0/3);
const scalar filterBCoeffBaseValue = 1.0/(filterR*len/3.464)/(filterR*len/3.464);
b.value() = filterBCoeffBaseValue;
//***********************************//MMA
const label mmaPowerConstraintIndex = 0;
const label mmaVolumeConstraintIndex = 1;
const bool mmaUseThermalEnvelopeConstraints =
  !turbulenceProModeActive
  &&
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
for (int i = 0; i < m; i++)
{
  dgdx[i].resize(n);
}
