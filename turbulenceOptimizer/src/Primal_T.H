T=T2;
fvScalarMatrix TEqn
(
    fvm::div(phi,T)
   -fvm::laplacian(alphaEff, T)
);
if(DDHS)
{
 TEqn-=thermalSource*(1-xh);
}
else
{
 TEqn-=thermalSource;
}
if (!turbulenceProModeActive)
{
 TEqn.relax();
}
const auto TLinearSolvePerformance =
  proModeRobustScalarSolverActive
? TEqn.solve(proModeRobustScalarSolverControls)
: TEqn.solve();
recordLinearSolve
(
  TLinearSolveStats,
  TLinearSolvePerformance.initialResidual(),
  TLinearSolvePerformance.finalResidual(),
  TLinearSolvePerformance.nIterations()
);
if (!turbulenceProModeActive)
{
 T2 = T;
}
