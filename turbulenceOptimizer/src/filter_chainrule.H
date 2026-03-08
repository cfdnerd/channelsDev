scalar fsensNormalization = gMax(mag(fsenshMeanT.primitiveField()));
scalar gsensPowerNormalization = gMax(mag(gsenshPowerDiss.primitiveField()));
scalar gsensThermalTmaxNormalization = 1.0;
scalar gsensThermalTminNormalization = 1.0;
if (useDedicatedThermalEnvelopeAdjoints)
{
  gsensThermalTmaxNormalization =
    gMax(mag(gsenshThermalTmax.primitiveField()));
  gsensThermalTminNormalization =
    gMax(mag(gsenshThermalTmin.primitiveField()));
}
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
if
(
    useDedicatedThermalEnvelopeAdjoints
 && (
       gsensThermalTmaxNormalization != gsensThermalTmaxNormalization
    || Foam::mag(gsensThermalTmaxNormalization) > GREAT
    || gsensThermalTmaxNormalization < SMALL
    )
)
{
  gsensThermalTmaxNormalization = 1.0;
}
if
(
    useDedicatedThermalEnvelopeAdjoints
 && (
       gsensThermalTminNormalization != gsensThermalTminNormalization
    || Foam::mag(gsensThermalTminNormalization) > GREAT
    || gsensThermalTminNormalization < SMALL
    )
)
{
  gsensThermalTminNormalization = 1.0;
}
fsenshMeanT = fsenshMeanT/fsensNormalization;
gsenshPowerDiss = gsenshPowerDiss/gsensPowerNormalization;
gsenshVol.primitiveFieldRef() = 1.0;
if (useDedicatedThermalEnvelopeAdjoints)
{
  gsenshThermalTmax = gsenshThermalTmax/gsensThermalTmaxNormalization;
  gsenshThermalTmin = gsenshThermalTmin/gsensThermalTminNormalization;
}
//***********************************//Heaviside_filter
if(opt>1)
{
  if (static_cast<bool>(optProjectionSmoothingEnable))
  {
    const scalar projectionBeta =
      Foam::max(optProjectionBetaScale*Foam::max(del, scalar(1.0)), scalar(1e-6));
    const scalar projectionDenominator =
      Foam::tanh(projectionBeta*eta5)
    + Foam::tanh(projectionBeta*(1.0 - eta5));
    const scalar projectionDenSafe =
      Foam::max(projectionDenominator, scalar(1e-6));
    for(i=0;i<n;i++)
    {
      const scalar tanhArg = Foam::tanh(projectionBeta*(xp[i] - eta5));
      drho[i] = projectionBeta*(1.0 - tanhArg*tanhArg)/projectionDenSafe;
    }
  }
  else
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
  }
  fsenshMeanT=fsenshMeanT*drho;
  gsenshPowerDiss=gsenshPowerDiss*drho;
  gsenshVol=gsenshVol*drho;
  if (useDedicatedThermalEnvelopeAdjoints)
  {
    gsenshThermalTmax = gsenshThermalTmax*drho;
    gsenshThermalTmin = gsenshThermalTmin*drho;
  }
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
if (useDedicatedThermalEnvelopeAdjoints)
{
  solve
  (
    fvm::laplacian(gsensThermalTmax)
   -fvm::Sp(b, gsensThermalTmax)
   +gsenshThermalTmax*b
  );
  solve
  (
    fvm::laplacian(gsensThermalTmin)
   -fvm::Sp(b, gsensThermalTmin)
   +gsenshThermalTmin*b
  );
}

const scalar maxCellVolume = Foam::max(gMax(mesh.V()), SMALL);
fsensMeanT.primitiveFieldRef()=fsensMeanT.primitiveFieldRef()*mesh.V()/maxCellVolume;
gsensPowerDiss.primitiveFieldRef()=gsensPowerDiss.primitiveFieldRef()*mesh.V()/maxCellVolume;
gsensVol.primitiveFieldRef()=gsensVol.primitiveFieldRef()*mesh.V()/maxCellVolume;
if (useDedicatedThermalEnvelopeAdjoints)
{
  gsensThermalTmax.primitiveFieldRef()=
    gsensThermalTmax.primitiveFieldRef()*mesh.V()/maxCellVolume;
  gsensThermalTmin.primitiveFieldRef()=
    gsensThermalTmin.primitiveFieldRef()*mesh.V()/maxCellVolume;
}
