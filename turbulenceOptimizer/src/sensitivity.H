
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
const volScalarField xhPositive(Foam::max(xh, scalar(0.0)));
const volScalarField xhMaterial
(
  Foam::pow(xhPositive, materialInterpolationExponentActive)
);
const volScalarField xhMaterialDerivative
(
  materialInterpolationExponentActive
 *Foam::pow(xhPositive, materialInterpolationExponentActive - scalar(1.0))
);

const volScalarField dAlphaDxh
(
  -alphaMax
  *(1 + qu)
  *qu
  /(
      (qu + xhMaterial)
     *(qu + xhMaterial)
    )
  *xhMaterialDerivative
);
const volScalarField dAlphaMolDxh
(
  (kf - ks)
  *(1 + qu)
  *qu
  /(
      (qu + xhMaterial)
     *(qu + xhMaterial)
    )
  /rhoCp
  *xhMaterialDerivative
);
const bool optimizationCycleTurbulentSensitivityMode =
    customTurbulenceActive && optimizationCycleTurbulentOptimizationMode;
const volScalarField& objectiveThermalAdjoint =
    optimizationCycleTurbulentSensitivityMode ? optimizationCycleThermalAdjointField : Tb;

fsenshMeanT = dAlphaDxh*(U&Ub)
           + b1*dAlphaMolDxh*(fvc::grad(T)&fvc::grad(objectiveThermalAdjoint));
if (static_cast<bool>(optGrayPenaltyEnable) && grayPenaltyWeightActive > SMALL)
{
 const dimensionedScalar grayPenaltyWeight
 (
   "grayPenaltyWeight",
   fsenshMeanT.dimensions(),
   grayPenaltyWeightActive
 );
 fsenshMeanT += grayPenaltyWeight*(scalar(1.0) - scalar(2.0)*xh);
}
if(DDHS)//design dependent heat source
{
 fsenshMeanT += b1*objectiveThermalAdjoint*thermalSource;
}
gsenshPowerDiss = dAlphaDxh*(U&Ua);

if (optimizationCycleTurbulentSensitivityMode)
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
std::vector<scalar> constraintGradientWeights(m, scalar(1.0));
constraintGradientWeights[mmaPowerConstraintIndex] = optPowerConstraintWeight;
constraintGradientWeights[mmaVolumeConstraintIndex] = optVolumeConstraintWeight;
if (turbulenceProModeActive)
{
 gx[0]=weightedPowerConstraint;
 gx[1]=weightedVolumeConstraint;

 for(i=0;i<n;i++)
 {
  xmma[i] = x[i];
  dfdx[i]=fsensMeanT[i]/N;//sensitivities of objective function
  dgdx[0][i]=gsensPowerDiss[i]/N;//sensitivities of constraint function
  dgdx[1][i]=gsensVol[i]/N;//sensitivities of constraint function
 }

 mma.asyminc = mma_inc;
 mma.asymdec = mma_dec;
 mma.movelimit = optimizerMoveLimitActive;
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
  dgdx[mmaPowerConstraintIndex][i]=gsensPowerDiss[i]/N;//sensitivities of constraint function
  dgdx[mmaVolumeConstraintIndex][i]=gsensVol[i]/N;//sensitivities of constraint function
  if (mmaUseThermalEnvelopeConstraints && m > mmaThermalTminConstraintIndex)
  {
    dgdx[mmaThermalTmaxConstraintIndex][i]=fsensMeanT[i]/N;
     dgdx[mmaThermalTminConstraintIndex][i]=-fsensMeanT[i]/N;
  }
 }
 mma.asyminc = mma_inc;
 mma.asymdec = mma_dec;
 mma.movelimit = optimizerMoveLimitActive;
}

// Gradient nondimensionalization keeps objective and constraints on comparable
// scales before MMA asymptote updates.
scalar dfdxL2Local = 0.0;
std::vector<scalar> dgdxL2Local(m, scalar(0.0));
for (i = 0; i < n; i++)
{
  const scalar dObj = static_cast<scalar>(dfdx[i]);
  dfdxL2Local += dObj*dObj;
  for (label constrainti = 0; constrainti < m; ++constrainti)
  {
    const scalar dConstraint = static_cast<scalar>(dgdx[constrainti][i]);
    dgdxL2Local[constrainti] += dConstraint*dConstraint;
  }
}
reduce(dfdxL2Local, sumOp<scalar>());
for (label constrainti = 0; constrainti < m; ++constrainti)
{
  reduce(dgdxL2Local[constrainti], sumOp<scalar>());
}
const scalar dfdxNorm =
  Foam::max(Foam::sqrt(dfdxL2Local), optSensitivityGradientL2Eps);
std::vector<scalar> dgdxNorm(m, optSensitivityGradientL2Eps);
for (label constrainti = 0; constrainti < m; ++constrainti)
{
  dgdxNorm[constrainti] =
    Foam::max(Foam::sqrt(dgdxL2Local[constrainti]), optSensitivityGradientL2Eps);
}
const scalar gradientSanityScale =
  static_cast<bool>(optGradientSanityRescaleEnable)
  ? optimizerGradientSanityScaleActive
  : scalar(1.0);
for (i = 0; i < n; i++)
{
  dfdx[i] = gradientSanityScale*(dfdx[i]/dfdxNorm);
  for (label constrainti = 0; constrainti < m; ++constrainti)
  {
    dgdx[constrainti][i] =
      gradientSanityScale
     *constraintGradientWeights[constrainti]
     *(dgdx[constrainti][i]/dgdxNorm[constrainti]);
  }
}
if (optimizerResetMMAAsymptotesNextIter)
{
  mma.iter = 1;
  optimizerResetMMAAsymptotesNextIter = false;
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
