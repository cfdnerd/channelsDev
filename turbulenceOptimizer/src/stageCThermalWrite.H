const scalar TaMinStageC = gMin(Ta);
const scalar TaMaxStageC = gMax(Ta);
const volScalarField UMagStageC(mag(U));
scalar UMeanStageCLocal = 0.0;
scalar TaMeanStageCLocal = 0.0;
forAll(UMagStageC, celli)
{
  UMeanStageCLocal += UMagStageC[celli]*mesh.V()[celli];
  TaMeanStageCLocal += Ta[celli]*mesh.V()[celli];
}
reduce(UMeanStageCLocal, sumOp<scalar>());
reduce(TaMeanStageCLocal, sumOp<scalar>());
const scalar UMeanStageC = UMeanStageCLocal/area;
const scalar TaMeanStageC = TaMeanStageCLocal/area;
const scalar UMaxStageC = gMax(UMagStageC);
const scalar dPStageC = gMax(p) - gMin(p);
const scalar TMaxStageC = gMax(T);
const scalar TMinStageC = gMin(T);
const scalar alphaEffMaxStageC = gMax(thermalDiffusivityBlock);
const scalar PkMaxStageC = gMax(PkField);
const scalar omegaProdMaxStageC = gMax(omegaProductionBlock);
const char* TaSolveStatus = classifyLinearSolveStatus(TaLinearSolveStats);
const bool stageCFullMode = customTurbulenceActive && stageCFullAdjointVerification;
const bool stageCFlowAdjointStalled =
    stageCFullMode
 && !stageCAdjointFlowConverged
 && !stageCAdjointFlowDiverged
 && stageCAdjointCorrectorsUsed >= stageCAdjointCorrectors;
const bool stageCUaTolMet = stageCAdjointUNonlinearChange <= stageCAdjointUNonlinearTol;
const bool stageCPaTolMet = stageCAdjointPNonlinearChange <= stageCAdjointPNonlinearTol;
const scalar stageCUaTargetRatio =
    stageCAdjointUNonlinearChange/Foam::max(stageCAdjointUNonlinearTol, SMALL);
const scalar stageCPaTargetRatio =
    stageCAdjointPNonlinearChange/Foam::max(stageCAdjointPNonlinearTol, SMALL);
const char* stageCFlowAdjointBlocker =
    (stageCUaTolMet && stageCPaTolMet)
  ? "none"
  : (stageCUaTargetRatio > 1.1*stageCPaTargetRatio)
    ? "Ua"
    : (stageCPaTargetRatio > 1.1*stageCUaTargetRatio)
      ? "pa"
      : "Ua+pa";
const char* stageCFlowViscosityMode = stageCUseNuEff ? "nuEff" : "nu";
const bool stageCNoImproveStop =
    stageCAdjointStopReason.find("no-improve") != std::string::npos;
const char* UaSolveStatus = classifyLinearSolveStatus(UaLinearSolveStats);
const char* paSolveStatus = classifyLinearSolveStatus(paLinearSolveStats);
if (stageCFlowAdjointStalled || stageCAdjointFlowDiverged)
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

scalar UaMaxStageC = 0.0;
scalar paRangeStageC = 0.0;
scalar kaMinStageC = 0.0;
scalar kaMaxStageC = 0.0;
scalar omegaaMinStageC = 0.0;
scalar omegaaMaxStageC = 0.0;
if (stageCFullMode)
{
  UaMaxStageC = stageCAdjointUaMax;
  paRangeStageC = stageCAdjointPaRange;
  kaMinStageC = gMin(ka);
  kaMaxStageC = gMax(ka);
  omegaaMinStageC = gMin(omegaa);
  omegaaMaxStageC = gMax(omegaa);
}

if (Pstream::master())
{
  std::ofstream debugStageCLogFileStream(debugStageCLogFile, std::ios::app);
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
    debugStageCLogFileStream
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
    debugStageCLogFileStream
      << "  "
      << std::left
      << std::setw(itemWidth) << item
      << std::setw(valueWidth) << v1
      << std::setw(valueWidth) << v2
      << std::setw(valueWidth) << v3
      << std::setw(valueWidth) << v4
      << '\n';
  };

  debugStageCLogFileStream
    << "\n=========================================================================================================\n"
    << " Stage-C Adjoint Verification Report\n"
    << "---------------------------------------------------------------------------------------------------------\n"
    << "  Iteration        : " << opt << "\n"
    << "  Mode             : "
    << (stageCFullMode ? "full-adjoint-verification" : "thermal-adjoint-only")
    << "\n"
    << "  Ta Solve Status  : " << TaSolveStatus << "\n"
    << "=========================================================================================================\n";

  writeSectionHeader("Primal State");
  writeRow("Velocity", "UMax="+fmtScalar(UMaxStageC), "UMean="+fmtScalar(UMeanStageC), "", "");
  writeRow("Pressure", "dP="+fmtScalar(dPStageC), "", "", "");
  writeRow("Temperature", "TMax="+fmtScalar(TMaxStageC), "TMin="+fmtScalar(TMinStageC), "", "");

  writeSectionHeader("Thermal Adjoint");
  writeRow("Ta field", "TaMin="+fmtScalar(TaMinStageC), "TaMax="+fmtScalar(TaMaxStageC), "TaMean="+fmtScalar(TaMeanStageC), "");
  writeRow("Ta solve", "final="+fmtScalar(TaLinearSolveStats.lastFinalResidual), "max="+fmtScalar(TaLinearSolveStats.finalResidualMax), "iters="+fmtCount(TaLinearSolveStats.iterationMax), "calls="+fmtCount(TaLinearSolveStats.callCount));

  if (stageCFullMode)
  {
    writeSectionHeader("Flow Adjoint");
    writeRow("Ua field", "|Ua|max="+fmtScalar(UaMaxStageC), "paRange="+fmtScalar(paRangeStageC), "", "");
    writeRow("Ua solve", "status="+std::string(UaSolveStatus), "final="+fmtScalar(UaLinearSolveStats.lastFinalResidual), "max="+fmtScalar(UaLinearSolveStats.finalResidualMax), "iters="+fmtCount(UaLinearSolveStats.iterationMax));
    writeRow("pa solve", "status="+std::string(paSolveStatus), "final="+fmtScalar(paLinearSolveStats.lastFinalResidual), "max="+fmtScalar(paLinearSolveStats.finalResidualMax), "iters="+fmtCount(paLinearSolveStats.iterationMax));
    writeRow("Ua loop", "used="+fmtCount(stageCAdjointCorrectorsUsed), "converged="+std::string(stageCAdjointFlowConverged ? "true" : "false"), "dU="+fmtScalar(stageCAdjointUNonlinearChange), "dP="+fmtScalar(stageCAdjointPNonlinearChange));
    writeRow("Ua target", "uTol="+fmtScalar(stageCAdjointUNonlinearTol), "uOK="+std::string(stageCUaTolMet ? "true" : "false"), "pTol="+fmtScalar(stageCAdjointPNonlinearTol), "pOK="+std::string(stageCPaTolMet ? "true" : "false"));
    writeRow("Ua gate", "uRatio="+fmtScalar(stageCUaTargetRatio), "pRatio="+fmtScalar(stageCPaTargetRatio), "blocker="+std::string(stageCFlowAdjointBlocker), "");
    writeRow("Ua stop", "reason="+std::string(stageCAdjointStopReason.c_str()), "bestPass="+fmtCount(stageCAdjointBestCorrector), std::string(stageCNoImproveStop ? "stagnated=" : "diverged=")+(stageCAdjointFlowDiverged ? "true" : "false"), "at="+fmtCount(stageCAdjointDivergedCorrector));
    writeRow("Ua rollback", "used="+std::string(stageCAdjointRolledBack ? "true" : "false"), "toPass="+fmtCount(stageCAdjointRollbackCorrector), "", "");
    writeRow("Ua best", "worstRatio="+fmtScalar(stageCAdjointBestCombinedMetric), "bestDU="+fmtScalar(stageCAdjointBestUNonlinearChange), "bestDP="+fmtScalar(stageCAdjointBestPNonlinearChange), "");
    writeRow("Ua config", "reset="+std::string(stageCResetAdjointState ? "true" : "false"), "visc="+std::string(stageCFlowViscosityMode), "TaSrc="+std::string(stageCIncludeThermalMomentumSource ? "true" : "false"), "legacyPRelax="+std::string(stageCUseLegacyPressureRelax ? "true" : "false"));
    writeRow("pa diag", "paMean="+fmtScalar(stageCAdjointPaMean), "|phia|sum="+fmtScalar(stageCAdjointPhiaAbsSum), "phiaNet="+fmtScalar(stageCAdjointPhiaNet), "pRef="+fmtCount(paRefCell)+":"+fmtScalar(paRefValue));

    writeSectionHeader("Turbulence Adjoint Proxies");
    writeRow("ka proxy", "kaMin="+fmtScalar(kaMinStageC), "kaMax="+fmtScalar(kaMaxStageC), "mode=lagged-explicit", "");
    writeRow("omegaa proxy", "wAMin="+fmtScalar(omegaaMinStageC), "wAMax="+fmtScalar(omegaaMaxStageC), "mode=lagged-explicit", "");
  }

  writeSectionHeader("Coupling");
  writeRow("Transport", "alphaEffMax="+fmtScalar(alphaEffMaxStageC), "PkMax="+fmtScalar(PkMaxStageC), "wProdMax="+fmtScalar(omegaProdMaxStageC), "");
  writeRow("Inlet", "status="+std::string(turbulenceInletConsistencyStatus), "kErr="+fmtScalar(turbulenceInletKRelativeError), "wErr="+fmtScalar(turbulenceInletOmegaRelativeError), "");
}
