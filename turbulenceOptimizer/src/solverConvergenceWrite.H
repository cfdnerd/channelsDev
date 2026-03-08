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
    const bool adjointCycleFlowAdjointStalled =
        customTurbulenceActive
     && adjointCycleFullAdjointVerification
     &&
     (
          adjointCycleAdjointFlowDiverged
       || (!adjointCycleAdjointFlowConverged
        && adjointCycleAdjointCorrectorsUsed >= adjointCycleAdjointCorrectors)
     );
    const bool adjointCycleFlowSolve =
        std::strcmp(solveName, "Ua adjoint momentum") == 0
     || std::strcmp(solveName, "pa adjoint pressure") == 0;
    const char* status = classifyLinearSolveStatus(stats);
    if
    (
        adjointCycleFlowAdjointStalled
     && adjointCycleFlowSolve
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
  writeSolveStats("ka adjoint closure", kaLinearSolveStats);
  writeSolveStats("omegaa adjoint closure", omegaaLinearSolveStats);
  writeSolveStats("T thermal", TLinearSolveStats);
  writeSolveStats("Ta adjoint thermal", TaLinearSolveStats);
  writeSolveStats("Tb adjoint thermal", TbLinearSolveStats);
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

  if (customTurbulenceActive && adjointCycleFullAdjointVerification)
  {
    const bool adjointCycleFlowAdjointStalled =
        !adjointCycleAdjointFlowConverged
     && !adjointCycleAdjointFlowDiverged
     && adjointCycleAdjointCorrectorsUsed >= adjointCycleAdjointCorrectors;
    const bool adjointCycleUaTolMet =
        adjointCycleAdjointUNonlinearChange <= adjointCycleAdjointUNonlinearTol;
    const bool adjointCyclePaTolMet =
        adjointCycleAdjointPNonlinearChange <= adjointCycleAdjointPNonlinearTol;
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
    const bool adjointCycleNoImproveStop =
        adjointCycleAdjointStopReason.find("no-improve") != std::string::npos;
    solverConvergenceLogFile
      << "\n"
      << "[Adjoint-Cycle Flow Adjoint Nonlinear Check]\n"
      << "  Status               : "
      << ((adjointCycleFlowAdjointStalled || adjointCycleAdjointFlowDiverged) ? "WARN" : "OK") << "\n"
      << "  Correctors Used      : " << adjointCycleAdjointCorrectorsUsed
      << " / " << adjointCycleAdjointCorrectors << "\n"
      << "  RMS Ua Change        : " << std::scientific << std::setprecision(3)
      << adjointCycleAdjointUNonlinearChange << "\n"
      << "  Ua Target / Met      : " << std::scientific << std::setprecision(3)
      << adjointCycleAdjointUNonlinearTol << " / "
      << (adjointCycleUaTolMet ? "true" : "false") << "\n"
      << "  RMS pa Change        : " << std::scientific << std::setprecision(3)
      << adjointCycleAdjointPNonlinearChange << "\n"
      << "  pa Target / Met      : " << std::scientific << std::setprecision(3)
      << adjointCycleAdjointPNonlinearTol << " / "
      << (adjointCyclePaTolMet ? "true" : "false") << "\n"
      << "  Dominant Blocker     : " << adjointCycleFlowAdjointBlocker << "\n"
      << "  Stop Reason          : " << adjointCycleAdjointStopReason << "\n"
      << "  Best Corrector       : " << adjointCycleAdjointBestCorrector << "\n"
      << "  Best Worst-Ratio     : " << std::scientific << std::setprecision(3)
      << adjointCycleAdjointBestCombinedMetric << "\n"
      << "  Stop Flag / At       : "
      << (adjointCycleNoImproveStop ? "stagnated" : "diverged")
      << "=" << (adjointCycleAdjointFlowDiverged ? "true" : "false")
      << " / " << adjointCycleAdjointDivergedCorrector << "\n"
      << "  Rolled Back / To     : "
      << (adjointCycleAdjointRolledBack ? "true" : "false")
      << " / " << adjointCycleAdjointRollbackCorrector << "\n"
      << "  |Ua|max / paRange    : " << std::scientific << std::setprecision(3)
      << adjointCycleAdjointUaMax << " / " << adjointCycleAdjointPaRange << "\n"
      << "  paMean               : " << std::scientific << std::setprecision(3)
      << adjointCycleAdjointPaMean << "\n"
      << "  |phia|sum / phiaNet  : " << std::scientific << std::setprecision(3)
      << adjointCycleAdjointPhiaAbsSum << " / " << adjointCycleAdjointPhiaNet << "\n"
      << "  pRefCell / Value     : " << paRefCell << " / "
      << std::scientific << std::setprecision(3) << paRefValue << "\n"
      << "  Config               : reset="
      << (adjointCycleResetAdjointState ? "true" : "false")
      << " viscosity=" << (adjointCycleUseNuEff ? "nuEff" : "nu")
      << " TaSource="
      << (adjointCycleIncludeThermalMomentumSource ? "true" : "false")
      << " legacyPaRelax="
      << (adjointCycleUseLegacyPressureRelax ? "true" : "false")
      << " divergenceGuard="
      << (adjointCycleDivergenceCheck ? "true" : "false")
      << " stagnationGuard="
      << (adjointCycleStagnationCheck ? "true" : "false")
      << " growth=" << adjointCycleDivergenceGrowthFactor
      << " minPass=" << adjointCycleDivergenceMinCorrectors
      << " patience=" << adjointCycleStagnationPatience << "\n";
  }
}
