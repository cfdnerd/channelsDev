T=T2;
const dimensionedScalar thermalPseudoTime
(
    "thermalPseudoTime",
    dimensionSet(0,0,-1,0,0,0,0),
    primalPseudoTimeCoeffActive
);
const dimensionedScalar thermalSourceRamp
(
    "thermalSourceRamp",
    dimless,
    thermalSourceRampActive
);
fvScalarMatrix TEqn
(
    fvm::div(phi,T)
   -fvm::laplacian(alphaEff, T)
   +fvm::Sp(thermalPseudoTime, T)
);
if(DDHS)
{
 TEqn-=thermalSource*thermalSourceRamp*(1-xh);
}
else
{
 TEqn-=thermalSource*thermalSourceRamp;
}
TEqn.relax(primalThermalUnderRelaxActive);
const auto TLinearSolvePerformance = TEqn.solve();
recordLinearSolve
(
  TLinearSolveStats,
  TLinearSolvePerformance.initialResidual(),
  TLinearSolvePerformance.finalResidual(),
  TLinearSolvePerformance.nIterations()
);
T2 = T;
