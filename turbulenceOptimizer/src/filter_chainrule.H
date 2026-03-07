scalar fsensNormalization = gMax(mag(fsenshMeanT.primitiveField()));
scalar gsensPowerNormalization = gMax(mag(gsenshPowerDiss.primitiveField()));
if
(
    fsensNormalization != fsensNormalization
 || Foam::mag(fsensNormalization) > GREAT
 || fsensNormalization < SMALL
)
{
  fsensNormalization = 1.0;
}
if
(
    gsensPowerNormalization != gsensPowerNormalization
 || Foam::mag(gsensPowerNormalization) > GREAT
 || gsensPowerNormalization < SMALL
)
{
  gsensPowerNormalization = 1.0;
}
fsenshMeanT = fsenshMeanT/fsensNormalization;
gsenshPowerDiss = gsenshPowerDiss/gsensPowerNormalization;
gsenshVol.primitiveFieldRef() = 1.0;
//***********************************//Heaviside_filter
if(opt>1)
{
  const scalar eta5Lower = Foam::max(eta5, scalar(1e-6));
  const scalar eta5Upper = Foam::max(1.0 - eta5, scalar(1e-6));
  for(i=0;i<n;i++)
  {
     if(xp[i]<=eta5)
     {
       drho[i]=del*Foam::exp(-del*(1-xp[i]/eta5Lower))+Foam::exp(-del); 
     }
     else
     {
       drho[i]=del*Foam::exp(-del*(xp[i]-eta5)/eta5Upper)+Foam::exp(-del); 
     } 
  }
  fsenshMeanT=fsenshMeanT*drho;
  gsenshPowerDiss=gsenshPowerDiss*drho;
  gsenshVol=gsenshVol*drho;
}
//***********************************//PDE filter
const auto fsensFilterLinearSolvePerformance =
solve(fvm::laplacian(fsensMeanT)-fvm::Sp(b,fsensMeanT)+fsenshMeanT*b);
recordLinearSolve
(
  fsensFilterLinearSolveStats,
  fsensFilterLinearSolvePerformance.initialResidual(),
  fsensFilterLinearSolvePerformance.finalResidual(),
  fsensFilterLinearSolvePerformance.nIterations()
);

const auto gsensPowerFilterLinearSolvePerformance =
solve(fvm::laplacian(gsensPowerDiss)-fvm::Sp(b,gsensPowerDiss)+gsenshPowerDiss*b);
recordLinearSolve
(
  gsensPowerFilterLinearSolveStats,
  gsensPowerFilterLinearSolvePerformance.initialResidual(),
  gsensPowerFilterLinearSolvePerformance.finalResidual(),
  gsensPowerFilterLinearSolvePerformance.nIterations()
);

const auto gsensVolFilterLinearSolvePerformance =
solve(fvm::laplacian(gsensVol)-fvm::Sp(b,gsensVol)+gsenshVol*b);
recordLinearSolve
(
  gsensVolFilterLinearSolveStats,
  gsensVolFilterLinearSolvePerformance.initialResidual(),
  gsensVolFilterLinearSolvePerformance.finalResidual(),
  gsensVolFilterLinearSolvePerformance.nIterations()
);

const scalar maxCellVolume = Foam::max(gMax(mesh.V()), SMALL);
fsensMeanT.primitiveFieldRef()=fsensMeanT.primitiveFieldRef()*mesh.V()/maxCellVolume;
gsensPowerDiss.primitiveFieldRef()=gsensPowerDiss.primitiveFieldRef()*mesh.V()/maxCellVolume;
gsensVol.primitiveFieldRef()=gsensVol.primitiveFieldRef()*mesh.V()/maxCellVolume;
