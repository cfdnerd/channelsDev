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
TEqn.relax();
const auto TLinearSolvePerformance = TEqn.solve();
recordLinearSolve
(
  TLinearSolveStats,
  TLinearSolvePerformance.initialResidual(),
  TLinearSolvePerformance.finalResidual(),
  TLinearSolvePerformance.nIterations()
);
T2 = T;
