
Info<<"sensitivity analysis"<<endl;

std::vector<scalar> xOld(n);
std::vector<scalar> xpOld(n);
std::vector<scalar> xhOld(n);
const vectorField UOldInternal(U.primitiveField());
const scalarField pOldInternal(p.primitiveField());
const scalarField TOldInternal(T.primitiveField());
const scalarField kOldInternal(k.primitiveField());
const scalarField omegaOldInternal(omega.primitiveField());
const scalarField nutOldInternal(nut.primitiveField());
const scalarField alphatOldInternal(alphat.primitiveField());
for(i=0;i<n;i++)
{
 xOld[i]=x[i];
 xpOld[i]=xp[i];
 xhOld[i]=xh[i];
}

const volScalarField dAlphaDxh
(
  -alphaMax*(1 + qu)*qu/((qu + xh)*(qu + xh))
);
const volScalarField dAlphaMolDxh
(
  (kf - ks)*(1 + qu)*qu/((qu + xh)*(qu + xh))/rhoCp
);
const bool optimizationCycleTurbulentSensitivityMode =
    customTurbulenceActive && optimizationCycleTurbulentOptimizationMode;
const volScalarField& objectiveThermalAdjoint =
    optimizationCycleTurbulentSensitivityMode ? optimizationCycleThermalAdjointField : Tb;

fsenshMeanT = dAlphaDxh*(U&Ub)
           + b1*dAlphaMolDxh*(fvc::grad(T)&fvc::grad(objectiveThermalAdjoint));
if(DDHS)//design dependent heat source
{
 fsenshMeanT += b1*objectiveThermalAdjoint*thermalSource;
}
gsenshPowerDiss = dAlphaDxh*(U&Ua);

if (optimizationCycleTurbulentSensitivityMode && optimizationCycleUseTaThermalAdjoint)
{
 #include "TurbulentSensitivity.H"
}


#include "filter_chainrule.H"

const scalar weightedPowerConstraint =
  optPowerConstraintWeight*
  (
    PowerDiss/Foam::max(PowerDissRelax - scalar(0.2)*opt, PowerDissMax)
   - 1.0
  );
const scalar weightedVolumeConstraint = optVolumeConstraintWeight*V;
if (turbulenceProModeActive)
{
 gx[0]=weightedPowerConstraint;
 gx[1]=weightedVolumeConstraint;

 for(i=0;i<n;i++)
 {
  xmma[i] = x[i];
  dfdx[i]=fsensMeanT[i]/N;//sensitivities of objective function
  dgdx[0][i]=optPowerConstraintWeight*(gsensPowerDiss[i]/N);//sensitivities of constraint function
  dgdx[1][i]=optVolumeConstraintWeight*(gsensVol[i]/N);//sensitivities of constraint function
 }

 mma.asyminc = mma_inc;
 mma.asymdec = mma_dec;
 mma.movelimit = movlim;
}
else
{
 gx[mmaPowerConstraintIndex]=weightedPowerConstraint;
 gx[mmaVolumeConstraintIndex]=weightedVolumeConstraint;
 if (mmaUseThermalEnvelopeConstraints && m > mmaThermalTminConstraintIndex)
 {
   gx[mmaThermalTmaxConstraintIndex] = thermalTmaxKSConstraint;
   gx[mmaThermalTminConstraintIndex] = thermalTminKSConstraint;
 }

 for(i=0;i<n;i++)
 {
  xmma[i] = x[i];
  dfdx[i]=fsensMeanT[i]/N;//sensitivities of objective function
  // Build raw constraint gradients first; apply user weights after L2 normalization.
  dgdx[mmaPowerConstraintIndex][i]=gsensPowerDiss[i]/N;//sensitivities of constraint function
  dgdx[mmaVolumeConstraintIndex][i]=gsensVol[i]/N;//sensitivities of constraint function
  if (mmaUseThermalEnvelopeConstraints && m > mmaThermalTminConstraintIndex)
  {
    dgdx[mmaThermalTmaxConstraintIndex][i]=fsensMeanT[i]/N;
    dgdx[mmaThermalTminConstraintIndex][i]=-fsensMeanT[i]/N;
  }
 }

 // Gradient nondimensionalization keeps objective and constraints on comparable
 // scales before MMA asymptote updates.
 scalar dfdxL2Local = 0.0;
 scalar dgdxPowerL2Local = 0.0;
 scalar dgdxVolumeL2Local = 0.0;
 for (i = 0; i < n; i++)
 {
   const scalar dObj = static_cast<scalar>(dfdx[i]);
   const scalar dPow = static_cast<scalar>(dgdx[mmaPowerConstraintIndex][i]);
   const scalar dVol = static_cast<scalar>(dgdx[mmaVolumeConstraintIndex][i]);
   dfdxL2Local += dObj*dObj;
   dgdxPowerL2Local += dPow*dPow;
   dgdxVolumeL2Local += dVol*dVol;
 }
 reduce(dfdxL2Local, sumOp<scalar>());
 reduce(dgdxPowerL2Local, sumOp<scalar>());
 reduce(dgdxVolumeL2Local, sumOp<scalar>());

 const scalar dfdxNorm =
   Foam::max(Foam::sqrt(dfdxL2Local), optSensitivityGradientL2Eps);
 const scalar dgdxPowerNorm =
   Foam::max(Foam::sqrt(dgdxPowerL2Local), optSensitivityGradientL2Eps);
 const scalar dgdxVolumeNorm =
   Foam::max(Foam::sqrt(dgdxVolumeL2Local), optSensitivityGradientL2Eps);
 scalar dgdxThermalTmaxNorm = optSensitivityGradientL2Eps;
 scalar dgdxThermalTminNorm = optSensitivityGradientL2Eps;
 if (mmaUseThermalEnvelopeConstraints && m > mmaThermalTminConstraintIndex)
 {
   scalar dgdxThermalTmaxL2Local = 0.0;
   scalar dgdxThermalTminL2Local = 0.0;
   for (i = 0; i < n; i++)
   {
     const scalar dTmax = static_cast<scalar>(dgdx[mmaThermalTmaxConstraintIndex][i]);
     const scalar dTmin = static_cast<scalar>(dgdx[mmaThermalTminConstraintIndex][i]);
     dgdxThermalTmaxL2Local += dTmax*dTmax;
     dgdxThermalTminL2Local += dTmin*dTmin;
   }
   reduce(dgdxThermalTmaxL2Local, sumOp<scalar>());
   reduce(dgdxThermalTminL2Local, sumOp<scalar>());
   dgdxThermalTmaxNorm =
     Foam::max(Foam::sqrt(dgdxThermalTmaxL2Local), optSensitivityGradientL2Eps);
   dgdxThermalTminNorm =
     Foam::max(Foam::sqrt(dgdxThermalTminL2Local), optSensitivityGradientL2Eps);
 }
 for (i = 0; i < n; i++)
 {
   dfdx[i] /= dfdxNorm;
   dgdx[mmaPowerConstraintIndex][i] =
     optPowerConstraintWeight*(dgdx[mmaPowerConstraintIndex][i]/dgdxPowerNorm);
   dgdx[mmaVolumeConstraintIndex][i] =
     optVolumeConstraintWeight*(dgdx[mmaVolumeConstraintIndex][i]/dgdxVolumeNorm);
   if (mmaUseThermalEnvelopeConstraints && m > mmaThermalTminConstraintIndex)
   {
     dgdx[mmaThermalTmaxConstraintIndex][i] /= dgdxThermalTmaxNorm;
     dgdx[mmaThermalTminConstraintIndex][i] /= dgdxThermalTminNorm;
   }
 }

 mma.asyminc = mma_inc;
 mma.asymdec = mma_dec;
 mma.movelimit = optimizerMoveLimitActive;
}

Info<< "\n run MMA \n" << endl;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
for(i=0;i<n;i++)
{
   x[i]=xmma[i];
}

#include "filter_x.H"
#include "debugOptimizer.H"
#include "solverConvergenceWrite.H"
