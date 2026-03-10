
//***********************************//PDE filter
const auto xpFilterLinearSolvePerformance =
solve(fvm::laplacian(xp)-fvm::Sp(b,xp)+x*b);
recordLinearSolve
(
  xpFilterLinearSolveStats,
  xpFilterLinearSolvePerformance.initialResidual(),
  xpFilterLinearSolvePerformance.finalResidual(),
  xpFilterLinearSolvePerformance.nIterations()
);
//***********************************//Heaviside_filter
if (turbulenceProModeActive)
{
  del=Foam::min(projectionDelStepActive*opt,100.0);
}
else
{
  const scalar targetDel =
    Foam::min(projectionDelStepActive*continuationScheduleIterationActive, scalar(100));
  const scalar continuationDelStep = projectionDelStepActive;
  del = Foam::min(del + continuationDelStep, targetDel);
}
eta0=0.0001;
eta1=1.0;
y0=diff(xp,mesh.V(),del,eta0,n);
reduce(y0, sumOp<scalar>());
do
{
  eta5=(eta0+eta1)/2.0;
  y5=diff(xp,mesh.V(),del,eta5,n);
  reduce(y5, sumOp<scalar>());
  if(y0*y5<0)
  {
     eta1=eta5;
  }
  else
  {
     eta0=eta5;
     y0=y5;
  }
} while ((eta1-eta0)>0.0001);

if (turbulenceProModeActive)
{
  for(i=0;i<n;i++)
  {
     if(xp[i]<=eta5)
     {
       xh[i]=eta5*(Foam::exp(-del*(1-xp[i]/eta5))-(1-xp[i]/eta5)*Foam::exp(-del));
     }
     else
     {
       xh[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+(xp[i]-eta5)*Foam::exp(-del)/(1-eta5));
     }
  }
}
else
{
  const scalar projectionBeta =
    Foam::max(projectionBetaScaleActive*Foam::max(del, scalar(1.0)), scalar(1e-6));
  const scalar projectionDenominator =
    Foam::tanh(projectionBeta*eta5)
  + Foam::tanh(projectionBeta*(1.0 - eta5));
  const scalar projectionDenSafe = Foam::max(projectionDenominator, scalar(1e-6));
  for(i=0;i<n;i++)
  {
     if (static_cast<bool>(optProjectionSmoothingEnable))
     {
       xh[i] =
         (
           Foam::tanh(projectionBeta*eta5)
         + Foam::tanh(projectionBeta*(xp[i] - eta5))
         )/projectionDenSafe;
     }
     else
     {
       if(xp[i]<=eta5)
       {
         xh[i]=eta5*(Foam::exp(-del*(1-xp[i]/eta5))-(1-xp[i]/eta5)*Foam::exp(-del));
       }
       else
       {
         xh[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+(xp[i]-eta5)*Foam::exp(-del)/(1-eta5));
       }
     }
  }
}
