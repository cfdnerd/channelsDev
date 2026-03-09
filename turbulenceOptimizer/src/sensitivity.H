
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
const bool useDedicatedThermalEnvelopeAdjoints =
    mmaUseThermalEnvelopeConstraints
 && m > mmaThermalTminConstraintIndex
 && static_cast<bool>(optDedicatedThermalEnvelopeAdjointsEnable);
volScalarField gsenshThermalTmax
(
  IOobject
  (
    "gsenshThermalTmax",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("zeroGsenshThermalTmax", fsenshMeanT.dimensions(), 0.0)
);
volScalarField gsenshThermalTmin
(
  IOobject
  (
    "gsenshThermalTmin",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
  ),
  mesh,
  dimensionedScalar("zeroGsenshThermalTmin", fsenshMeanT.dimensions(), 0.0)
);
volScalarField gsensThermalTmax(gsenshThermalTmax);
volScalarField gsensThermalTmin(gsenshThermalTmin);

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

if (useDedicatedThermalEnvelopeAdjoints)
{
  auto auditFiniteScalarField =
  [&](const volScalarField& field, const char* context)
  {
    const scalar fieldMin = gMin(field);
    const scalar fieldMax = gMax(field);

    label localInvalidCount = 0;
    label localFirstInvalidCell = -1;
    scalar localFirstInvalidValue = 0.0;
    const scalarField& fieldInternal = field.internalField();
    forAll(fieldInternal, celli)
    {
      const scalar value = fieldInternal[celli];
      const bool invalidValue =
           !std::isfinite(static_cast<double>(value))
        || Foam::mag(value) > GREAT;
      if (!invalidValue)
      {
        continue;
      }

      localInvalidCount++;
      if (localFirstInvalidCell < 0)
      {
        localFirstInvalidCell = celli;
        localFirstInvalidValue = value;
      }
    }

    label globalInvalidCount = localInvalidCount;
    reduce(globalInvalidCount, sumOp<label>());
    if
    (
         std::isfinite(static_cast<double>(fieldMin))
      && std::isfinite(static_cast<double>(fieldMax))
      && globalInvalidCount == 0
    )
    {
      return;
    }

    FatalErrorInFunction
      << "Non-finite field detected during dedicated thermal-envelope adjoint assembly."
      << " Context=" << context
      << " min=" << fieldMin
      << " max=" << fieldMax
      << " globalInvalidCount=" << globalInvalidCount
      << " localInvalidCount=" << localInvalidCount
      << " localFirstInvalidCell=" << localFirstInvalidCell
      << " localFirstInvalidValue=" << localFirstInvalidValue
      << " proc=" << Pstream::myProcNo()
      << exit(FatalError);
  };

  scalar thetaMaxPeakLocal = -GREAT;
  scalar thetaMinPeakLocal = -GREAT;
  for (i = 0; i < n; i++)
  {
    const scalar Ti = T[i];
    const scalar thetaMax = (Ti - optThermalTmaxLimit)/Foam::max(optThermalExcessSoftplusScale, SMALL);
    const scalar thetaMin = (optThermalTminLimit - Ti)/Foam::max(optThermalExcessSoftplusScale, SMALL);
    thetaMaxPeakLocal = Foam::max(thetaMaxPeakLocal, thetaMax);
    thetaMinPeakLocal = Foam::max(thetaMinPeakLocal, thetaMin);
  }
  reduce(thetaMaxPeakLocal, maxOp<scalar>());
  reduce(thetaMinPeakLocal, maxOp<scalar>());

  scalar ksMaxExpSumLocal = 0.0;
  scalar ksMinExpSumLocal = 0.0;
  for (i = 0; i < n; i++)
  {
    const scalar Ti = T[i];
    const scalar cellV = mesh.V()[i];
    const scalar thetaMax = (Ti - optThermalTmaxLimit)/Foam::max(optThermalExcessSoftplusScale, SMALL);
    const scalar thetaMin = (optThermalTminLimit - Ti)/Foam::max(optThermalExcessSoftplusScale, SMALL);
    ksMaxExpSumLocal += std::exp(thermalKSBeta*(thetaMax - thetaMaxPeakLocal))*cellV;
    ksMinExpSumLocal += std::exp(thermalKSBeta*(thetaMin - thetaMinPeakLocal))*cellV;
  }
  reduce(ksMaxExpSumLocal, sumOp<scalar>());
  reduce(ksMinExpSumLocal, sumOp<scalar>());
  const scalar ksMaxExpDen = Foam::max(ksMaxExpSumLocal, SMALL);
  const scalar ksMinExpDen = Foam::max(ksMinExpSumLocal, SMALL);
  const scalar thermalScale = Foam::max(optThermalExcessSoftplusScale, SMALL);

  volScalarField thermalTmaxWeight
  (
    IOobject
    (
      "thermalTmaxConstraintWeight",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroThermalTmaxWeight", dimless, 0.0)
  );
  volScalarField thermalTminWeight
  (
    IOobject
    (
      "thermalTminConstraintWeight",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroThermalTminWeight", dimless, 0.0)
  );
  scalarField& thermalTmaxWeightInternal = thermalTmaxWeight.primitiveFieldRef();
  scalarField& thermalTminWeightInternal = thermalTminWeight.primitiveFieldRef();
  const scalarField& TInternal = T.internalField();
  const scalarField& VInternal = mesh.V();
  forAll(thermalTmaxWeightInternal, celli)
  {
    const scalar Ti = TInternal[celli];
    const scalar cellV = VInternal[celli];
    const scalar thetaMax = (Ti - optThermalTmaxLimit)/thermalScale;
    const scalar thetaMin = (optThermalTminLimit - Ti)/thermalScale;
    const scalar expMax = std::exp(thermalKSBeta*(thetaMax - thetaMaxPeakLocal));
    const scalar expMin = std::exp(thermalKSBeta*(thetaMin - thetaMinPeakLocal));
    thermalTmaxWeightInternal[celli] = (expMax*cellV/ksMaxExpDen)/thermalScale;
    thermalTminWeightInternal[celli] = -(expMin*cellV/ksMinExpDen)/thermalScale;
  }
  thermalTmaxWeight.correctBoundaryConditions();
  thermalTminWeight.correctBoundaryConditions();
  auditFiniteScalarField(thermalTmaxWeight, "thermalTmaxWeight");
  auditFiniteScalarField(thermalTminWeight, "thermalTminWeight");
  auditFiniteScalarField(DT, "DT");

  volScalarField TbThermalTmax
  (
    IOobject
    (
      "TbThermalTmax",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    objectiveThermalAdjoint
  );
  TbThermalTmax =
    dimensionedScalar("zeroTbThermalTmax", T.dimensions(), 0.0);
  volScalarField TbThermalTmin
  (
    IOobject
    (
      "TbThermalTmin",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    objectiveThermalAdjoint
  );
  TbThermalTmin =
    dimensionedScalar("zeroTbThermalTmin", T.dimensions(), 0.0);
  const dimensionedScalar thermalAdjointDiagReg
  (
    "thermalAdjointDiagReg",
    dimensionSet(0,0,-1,0,0,0,0),
    1e-12
  );

  fvScalarMatrix TbThermalTmaxEqn
  (
     -fvm::div(-phi, TbThermalTmax)
     +fvm::laplacian(DT, TbThermalTmax)
     +fvm::Sp(thermalAdjointDiagReg, TbThermalTmax)
  );
  TbThermalTmaxEqn -= unitThermalSource*thermalTmaxWeight;
  Info<< "Solving dedicated thermal-envelope adjoint TbThermalTmax\n" << endl;
  TbThermalTmaxEqn.solve();

  fvScalarMatrix TbThermalTminEqn
  (
     -fvm::div(-phi, TbThermalTmin)
     +fvm::laplacian(DT, TbThermalTmin)
     +fvm::Sp(thermalAdjointDiagReg, TbThermalTmin)
  );
  TbThermalTminEqn -= unitThermalSource*thermalTminWeight;
  Info<< "Solving dedicated thermal-envelope adjoint TbThermalTmin\n" << endl;
  TbThermalTminEqn.solve();

  gsenshThermalTmax =
      dAlphaDxh*(U&Ub)
    + b1*dAlphaMolDxh*(fvc::grad(T)&fvc::grad(TbThermalTmax));
  gsenshThermalTmin =
      dAlphaDxh*(U&Ub)
    + b1*dAlphaMolDxh*(fvc::grad(T)&fvc::grad(TbThermalTmin));
  if (DDHS)
  {
    gsenshThermalTmax += b1*TbThermalTmax*thermalSource;
    gsenshThermalTmin += b1*TbThermalTmin*thermalSource;
  }
}


#include "filter_chainrule.H"

const scalar weightedPowerConstraint =
  optPowerConstraintWeight*
  (
    PowerDiss/Foam::max(PowerDissRelax - scalar(0.2)*opt, PowerDissMax)
   - 1.0
  );
const scalar weightedVolumeConstraint = optVolumeConstraintWeight*V;
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
   if (useDedicatedThermalEnvelopeAdjoints)
   {
     dgdx[mmaThermalTmaxConstraintIndex][i]=gsensThermalTmax[i]/N;
     dgdx[mmaThermalTminConstraintIndex][i]=gsensThermalTmin[i]/N;
   }
   else
   {
     // Fallback surrogate map for backward compatibility.
     dgdx[mmaThermalTmaxConstraintIndex][i]=fsensMeanT[i]/N;
     dgdx[mmaThermalTminConstraintIndex][i]=-fsensMeanT[i]/N;
   }
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

Info<< "\n run MMA \n" << endl;
mma.movelimit = optimizerMoveLimitActive;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
for(i=0;i<n;i++)
{
   x[i]=xmma[i];
}

#include "filter_x.H"
#include "debugOptimizer.H"
#include "solverConvergenceWrite.H"
