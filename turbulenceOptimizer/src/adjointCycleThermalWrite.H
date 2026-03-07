const scalar TaMinAdjointCycle = gMin(Ta);
const scalar TaMaxAdjointCycle = gMax(Ta);
const volScalarField UMagAdjointCycle(mag(U));
scalar UMeanAdjointCycleLocal = 0.0;
scalar TaMeanAdjointCycleLocal = 0.0;
forAll(UMagAdjointCycle, celli)
{
  UMeanAdjointCycleLocal += UMagAdjointCycle[celli]*mesh.V()[celli];
  TaMeanAdjointCycleLocal += Ta[celli]*mesh.V()[celli];
}
reduce(UMeanAdjointCycleLocal, sumOp<scalar>());
reduce(TaMeanAdjointCycleLocal, sumOp<scalar>());
const scalar UMeanAdjointCycle = UMeanAdjointCycleLocal/area;
const scalar TaMeanAdjointCycle = TaMeanAdjointCycleLocal/area;
const scalar UMaxAdjointCycle = gMax(UMagAdjointCycle);
const scalar dPAdjointCycle = gMax(p) - gMin(p);
const scalar TMaxAdjointCycle = gMax(T);
const scalar TMinAdjointCycle = gMin(T);
const scalar alphaEffMaxAdjointCycle = gMax(thermalDiffusivityBlock);
const scalar PkMaxAdjointCycle = gMax(PkField);
const scalar omegaProdMaxAdjointCycle = gMax(omegaProductionBlock);
const char* TaSolveStatus = classifyLinearSolveStatus(TaLinearSolveStats);
const bool adjointCycleFullMode = customTurbulenceActive && adjointCycleFullAdjointVerification;
const bool adjointCycleFlowAdjointStalled =
    adjointCycleFullMode
 && !adjointCycleAdjointFlowConverged
 && !adjointCycleAdjointFlowDiverged
 && adjointCycleAdjointCorrectorsUsed >= adjointCycleAdjointCorrectors;
const bool adjointCycleUaTolMet = adjointCycleAdjointUNonlinearChange <= adjointCycleAdjointUNonlinearTol;
const bool adjointCyclePaTolMet = adjointCycleAdjointPNonlinearChange <= adjointCycleAdjointPNonlinearTol;
const scalar adjointCycleUaTargetRatio =
    adjointCycleAdjointUNonlinearChange/Foam::max(adjointCycleAdjointUNonlinearTol, SMALL);
const scalar adjointCyclePaTargetRatio =
    adjointCycleAdjointPNonlinearChange/Foam::max(adjointCycleAdjointPNonlinearTol, SMALL);
const char* adjointCycleFlowAdjointBlocker =
    (adjointCycleUaTolMet && adjointCyclePaTolMet)
  ? "none"
  : (adjointCycleUaTargetRatio > 1.1*adjointCyclePaTargetRatio)
    ? "Ua"
    : (adjointCyclePaTargetRatio > 1.1*adjointCycleUaTargetRatio)
      ? "pa"
      : "Ua+pa";
const char* adjointCycleFlowViscosityMode = adjointCycleUseNuEff ? "nuEff" : "nu";
const bool adjointCycleNoImproveStop =
    adjointCycleAdjointStopReason.find("no-improve") != std::string::npos;
const char* UaSolveStatus = classifyLinearSolveStatus(UaLinearSolveStats);
const char* paSolveStatus = classifyLinearSolveStatus(paLinearSolveStats);
if (adjointCycleFlowAdjointStalled || adjointCycleAdjointFlowDiverged)
{
  if (std::strcmp(UaSolveStatus, "WARN") != 0)
  {
    UaSolveStatus = "WARN";
  }
  if (std::strcmp(paSolveStatus, "WARN") != 0)
  {
    paSolveStatus = "WARN";
  }
}

scalar UaMaxAdjointCycle = 0.0;
scalar paRangeAdjointCycle = 0.0;
scalar kaMinAdjointCycle = 0.0;
scalar kaMaxAdjointCycle = 0.0;
scalar omegaaMinAdjointCycle = 0.0;
scalar omegaaMaxAdjointCycle = 0.0;
if (adjointCycleFullMode)
{
  UaMaxAdjointCycle = adjointCycleAdjointUaMax;
  paRangeAdjointCycle = adjointCycleAdjointPaRange;
  kaMinAdjointCycle = gMin(ka);
  kaMaxAdjointCycle = gMax(ka);
  omegaaMinAdjointCycle = gMin(omegaa);
  omegaaMaxAdjointCycle = gMax(omegaa);
}

if (Pstream::master())
{
  std::ofstream debugAdjointCycleLogFileStream(debugAdjointCycleLogFile, std::ios::app);
  const int itemWidth = 24;
  const int valueWidth = 36;

  auto fmtScalar =
  [&](const scalar value)
  {
    std::ostringstream os;
    os << std::scientific << std::setprecision(6) << value;
    return os.str();
  };

  auto fmtCount =
  [&](const label value)
  {
    std::ostringstream os;
    os << value;
    return os.str();
  };

  auto writeSectionHeader =
  [&](const char* title)
  {
    debugAdjointCycleLogFileStream
      << "\n[" << title << "]\n"
      << "  "
      << std::left
      << std::setw(itemWidth) << "Item"
      << std::setw(valueWidth) << "Value 1"
      << std::setw(valueWidth) << "Value 2"
      << std::setw(valueWidth) << "Value 3"
      << std::setw(valueWidth) << "Value 4"
      << '\n'
      << "  "
      << std::string(itemWidth + 4*valueWidth, '-')
      << '\n';
  };

  auto writeRow =
  [&](const std::string& item, const std::string& v1, const std::string& v2 = "", const std::string& v3 = "", const std::string& v4 = "")
  {
    debugAdjointCycleLogFileStream
      << "  "
      << std::left
      << std::setw(itemWidth) << item
      << std::setw(valueWidth) << v1
      << std::setw(valueWidth) << v2
      << std::setw(valueWidth) << v3
      << std::setw(valueWidth) << v4
      << '\n';
  };

  debugAdjointCycleLogFileStream
    << "\n=========================================================================================================\n"
    << " Adjoint-Cycle Adjoint Verification Report\n"
    << "---------------------------------------------------------------------------------------------------------\n"
    << "  Iteration        : " << opt << "\n"
    << "  Mode             : "
    << (adjointCycleFullMode ? "full-adjoint-verification" : "thermal-adjoint-only")
    << "\n"
    << "  Ta Solve Status  : " << TaSolveStatus << "\n"
    << "=========================================================================================================\n";

  writeSectionHeader("Primal State");
  writeRow("Velocity", "UMax="+fmtScalar(UMaxAdjointCycle), "UMean="+fmtScalar(UMeanAdjointCycle), "", "");
  writeRow("Pressure", "dP="+fmtScalar(dPAdjointCycle), "", "", "");
  writeRow("Temperature", "TMax="+fmtScalar(TMaxAdjointCycle), "TMin="+fmtScalar(TMinAdjointCycle), "", "");

  writeSectionHeader("Thermal Adjoint");
  writeRow("Ta field", "TaMin="+fmtScalar(TaMinAdjointCycle), "TaMax="+fmtScalar(TaMaxAdjointCycle), "TaMean="+fmtScalar(TaMeanAdjointCycle), "");
  writeRow("Ta solve", "final="+fmtScalar(TaLinearSolveStats.lastFinalResidual), "max="+fmtScalar(TaLinearSolveStats.finalResidualMax), "iters="+fmtCount(TaLinearSolveStats.iterationMax), "calls="+fmtCount(TaLinearSolveStats.callCount));

  if (adjointCycleFullMode)
  {
    writeSectionHeader("Flow Adjoint");
    writeRow("Ua field", "|Ua|max="+fmtScalar(UaMaxAdjointCycle), "paRange="+fmtScalar(paRangeAdjointCycle), "", "");
    writeRow("Ua solve", "status="+std::string(UaSolveStatus), "final="+fmtScalar(UaLinearSolveStats.lastFinalResidual), "max="+fmtScalar(UaLinearSolveStats.finalResidualMax), "iters="+fmtCount(UaLinearSolveStats.iterationMax));
    writeRow("pa solve", "status="+std::string(paSolveStatus), "final="+fmtScalar(paLinearSolveStats.lastFinalResidual), "max="+fmtScalar(paLinearSolveStats.finalResidualMax), "iters="+fmtCount(paLinearSolveStats.iterationMax));
    writeRow("Ua loop", "used="+fmtCount(adjointCycleAdjointCorrectorsUsed), "converged="+std::string(adjointCycleAdjointFlowConverged ? "true" : "false"), "dU="+fmtScalar(adjointCycleAdjointUNonlinearChange), "dP="+fmtScalar(adjointCycleAdjointPNonlinearChange));
    writeRow("Ua target", "uTol="+fmtScalar(adjointCycleAdjointUNonlinearTol), "uOK="+std::string(adjointCycleUaTolMet ? "true" : "false"), "pTol="+fmtScalar(adjointCycleAdjointPNonlinearTol), "pOK="+std::string(adjointCyclePaTolMet ? "true" : "false"));
    writeRow("Ua gate", "uRatio="+fmtScalar(adjointCycleUaTargetRatio), "pRatio="+fmtScalar(adjointCyclePaTargetRatio), "blocker="+std::string(adjointCycleFlowAdjointBlocker), "");
    writeRow("Ua stop", "reason="+std::string(adjointCycleAdjointStopReason.c_str()), "bestPass="+fmtCount(adjointCycleAdjointBestCorrector), std::string(adjointCycleNoImproveStop ? "stagnated=" : "diverged=")+(adjointCycleAdjointFlowDiverged ? "true" : "false"), "at="+fmtCount(adjointCycleAdjointDivergedCorrector));
    writeRow("Ua rollback", "used="+std::string(adjointCycleAdjointRolledBack ? "true" : "false"), "toPass="+fmtCount(adjointCycleAdjointRollbackCorrector), "", "");
    writeRow("Ua best", "worstRatio="+fmtScalar(adjointCycleAdjointBestCombinedMetric), "bestDU="+fmtScalar(adjointCycleAdjointBestUNonlinearChange), "bestDP="+fmtScalar(adjointCycleAdjointBestPNonlinearChange), "");
    writeRow("Ua config", "reset="+std::string(adjointCycleResetAdjointState ? "true" : "false"), "visc="+std::string(adjointCycleFlowViscosityMode), "TaSrc="+std::string(adjointCycleIncludeThermalMomentumSource ? "true" : "false"), "legacyPRelax="+std::string(adjointCycleUseLegacyPressureRelax ? "true" : "false"));
    writeRow("pa diag", "paMean="+fmtScalar(adjointCycleAdjointPaMean), "|phia|sum="+fmtScalar(adjointCycleAdjointPhiaAbsSum), "phiaNet="+fmtScalar(adjointCycleAdjointPhiaNet), "pRef="+fmtCount(paRefCell)+":"+fmtScalar(paRefValue));

    writeSectionHeader("Turbulence Adjoint Proxies");
    writeRow("ka proxy", "kaMin="+fmtScalar(kaMinAdjointCycle), "kaMax="+fmtScalar(kaMaxAdjointCycle), "mode=lagged-explicit", "");
    writeRow("omegaa proxy", "wAMin="+fmtScalar(omegaaMinAdjointCycle), "wAMax="+fmtScalar(omegaaMaxAdjointCycle), "mode=lagged-explicit", "");
  }

  writeSectionHeader("Coupling");
  writeRow("Transport", "alphaEffMax="+fmtScalar(alphaEffMaxAdjointCycle), "PkMax="+fmtScalar(PkMaxAdjointCycle), "wProdMax="+fmtScalar(omegaProdMaxAdjointCycle), "");
  writeRow("Inlet", "status="+std::string(turbulenceInletConsistencyStatus), "kErr="+fmtScalar(turbulenceInletKRelativeError), "wErr="+fmtScalar(turbulenceInletOmegaRelativeError), "");
}
