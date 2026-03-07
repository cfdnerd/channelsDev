if (Pstream::master())
{
  std::ofstream solverConvergenceLogFile(solverConvergenceFile, std::ios::app);

  auto writeTableHeader =
  [&]()
  {
    solverConvergenceLogFile
      << "  "
      << std::left
      << std::setw(22) << "Solve"
      << std::setw(6) << "Stat"
      << std::right
      << std::setw(7) << "Calls"
      << std::setw(14) << "AvgInit"
      << std::setw(14) << "MaxInit"
      << std::setw(14) << "AvgFinal"
      << std::setw(14) << "MaxFinal"
      << std::setw(10) << "AvgIter"
      << std::setw(10) << "MaxIter"
      << std::setw(14) << "LastFinal"
      << std::setw(10) << "LastItr"
      << "\n"
      << "  "
      << std::string(22 + 6 + 7 + 14 + 14 + 14 + 14 + 10 + 10 + 14 + 10, '-')
      << "\n";
  };

  auto writeSolveStats =
  [&](const char* solveName, const LinearSolveStats& stats)
  {
    const scalar avgInitialResidual =
      stats.callCount > 0 ? stats.initialResidualSum/static_cast<scalar>(stats.callCount) : 0;
    const scalar avgFinalResidual =
      stats.callCount > 0 ? stats.finalResidualSum/static_cast<scalar>(stats.callCount) : 0;
    const scalar avgIterations =
      stats.callCount > 0 ? static_cast<scalar>(stats.iterationSum)/static_cast<scalar>(stats.callCount) : 0;
    const bool stageCFlowAdjointStalled =
        customTurbulenceActive
     && stageCFullAdjointVerification
     &&
     (
          stageCAdjointFlowDiverged
       || (!stageCAdjointFlowConverged
        && stageCAdjointCorrectorsUsed >= stageCAdjointCorrectors)
     );
    const bool stageCFlowSolve =
        std::strcmp(solveName, "Ua adjoint momentum") == 0
     || std::strcmp(solveName, "pa adjoint pressure") == 0;
    const char* status = classifyLinearSolveStatus(stats);
    if
    (
        stageCFlowAdjointStalled
     && stageCFlowSolve
     && std::strcmp(status, "WARN") != 0
    )
    {
      status = "WARN";
    }

    solverConvergenceLogFile
      << "  "
      << std::left << std::setw(22) << solveName
      << std::setw(6) << status
      << std::right
      << std::setw(7) << stats.callCount
      << std::scientific << std::setprecision(3)
      << std::setw(14) << avgInitialResidual
      << std::setw(14) << stats.initialResidualMax
      << std::setw(14) << avgFinalResidual
      << std::setw(14) << stats.finalResidualMax
      << std::fixed << std::setprecision(2)
      << std::setw(10) << avgIterations
      << std::setw(10) << stats.iterationMax
      << std::scientific << std::setprecision(3)
      << std::setw(14) << stats.lastFinalResidual
      << std::fixed << std::setprecision(0)
      << std::setw(10) << stats.lastIterations
      << "\n";
  };

  solverConvergenceLogFile
    << "\n"
    << "=========================================================================================================\n"
    << " Solver Convergence Report\n"
    << "---------------------------------------------------------------------------------------------------------\n"
    << "  Optimization Iteration : " << opt << "\n"
    << "  Status Legend         : OK=healthy  OBS=watch  WARN=likely issue  SKIP=not called\n"
    << "=========================================================================================================\n"
    << "[Primary Equations]\n";
  writeTableHeader();
  writeSolveStats("U momentum", ULinearSolveStats);
  writeSolveStats("p pressure", pLinearSolveStats);
  writeSolveStats("k turbulence", kLinearSolveStats);
  writeSolveStats("omega turbulence", omegaLinearSolveStats);
  writeSolveStats("T thermal", TLinearSolveStats);
  writeSolveStats("Ta adjoint thermal", TaLinearSolveStats);
  writeSolveStats("Ub adjoint momentum", UbLinearSolveStats);
  writeSolveStats("pb adjoint pressure", pbLinearSolveStats);
  writeSolveStats("Ua adjoint momentum", UaLinearSolveStats);
  writeSolveStats("pa adjoint pressure", paLinearSolveStats);

  solverConvergenceLogFile
    << "\n"
    << "[Filter And Projection Equations]\n";
  writeTableHeader();
  writeSolveStats("xp filter", xpFilterLinearSolveStats);
  writeSolveStats("fsens filter", fsensFilterLinearSolveStats);
  writeSolveStats("gsensPower filter", gsensPowerFilterLinearSolveStats);
  writeSolveStats("gsensVol filter", gsensVolFilterLinearSolveStats);

  if (customTurbulenceActive && stageCFullAdjointVerification)
  {
    const bool stageCFlowAdjointStalled =
        !stageCAdjointFlowConverged
     && !stageCAdjointFlowDiverged
     && stageCAdjointCorrectorsUsed >= stageCAdjointCorrectors;
    const bool stageCUaTolMet =
        stageCAdjointUNonlinearChange <= stageCAdjointUNonlinearTol;
    const bool stageCPaTolMet =
        stageCAdjointPNonlinearChange <= stageCAdjointPNonlinearTol;
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
    const bool stageCNoImproveStop =
        stageCAdjointStopReason.find("no-improve") != std::string::npos;
    solverConvergenceLogFile
      << "\n"
      << "[Stage-C Flow Adjoint Nonlinear Check]\n"
      << "  Status               : "
      << ((stageCFlowAdjointStalled || stageCAdjointFlowDiverged) ? "WARN" : "OK") << "\n"
      << "  Correctors Used      : " << stageCAdjointCorrectorsUsed
      << " / " << stageCAdjointCorrectors << "\n"
      << "  RMS Ua Change        : " << std::scientific << std::setprecision(3)
      << stageCAdjointUNonlinearChange << "\n"
      << "  Ua Target / Met      : " << std::scientific << std::setprecision(3)
      << stageCAdjointUNonlinearTol << " / "
      << (stageCUaTolMet ? "true" : "false") << "\n"
      << "  RMS pa Change        : " << std::scientific << std::setprecision(3)
      << stageCAdjointPNonlinearChange << "\n"
      << "  pa Target / Met      : " << std::scientific << std::setprecision(3)
      << stageCAdjointPNonlinearTol << " / "
      << (stageCPaTolMet ? "true" : "false") << "\n"
      << "  Dominant Blocker     : " << stageCFlowAdjointBlocker << "\n"
      << "  Stop Reason          : " << stageCAdjointStopReason << "\n"
      << "  Best Corrector       : " << stageCAdjointBestCorrector << "\n"
      << "  Best Worst-Ratio     : " << std::scientific << std::setprecision(3)
      << stageCAdjointBestCombinedMetric << "\n"
      << "  Stop Flag / At       : "
      << (stageCNoImproveStop ? "stagnated" : "diverged")
      << "=" << (stageCAdjointFlowDiverged ? "true" : "false")
      << " / " << stageCAdjointDivergedCorrector << "\n"
      << "  Rolled Back / To     : "
      << (stageCAdjointRolledBack ? "true" : "false")
      << " / " << stageCAdjointRollbackCorrector << "\n"
      << "  |Ua|max / paRange    : " << std::scientific << std::setprecision(3)
      << stageCAdjointUaMax << " / " << stageCAdjointPaRange << "\n"
      << "  paMean               : " << std::scientific << std::setprecision(3)
      << stageCAdjointPaMean << "\n"
      << "  |phia|sum / phiaNet  : " << std::scientific << std::setprecision(3)
      << stageCAdjointPhiaAbsSum << " / " << stageCAdjointPhiaNet << "\n"
      << "  pRefCell / Value     : " << paRefCell << " / "
      << std::scientific << std::setprecision(3) << paRefValue << "\n"
      << "  Config               : reset="
      << (stageCResetAdjointState ? "true" : "false")
      << " viscosity=" << (stageCUseNuEff ? "nuEff" : "nu")
      << " TaSource="
      << (stageCIncludeThermalMomentumSource ? "true" : "false")
      << " legacyPaRelax="
      << (stageCUseLegacyPressureRelax ? "true" : "false")
      << " divergenceGuard="
      << (stageCDivergenceCheck ? "true" : "false")
      << " stagnationGuard="
      << (stageCStagnationCheck ? "true" : "false")
      << " growth=" << stageCDivergenceGrowthFactor
      << " minPass=" << stageCDivergenceMinCorrectors
      << " patience=" << stageCStagnationPatience << "\n";
  }
}
