//Tb=T2;
sumT=0;
maxT=gMax(T);
const scalar maxTSafe = Foam::max(Foam::mag(maxT), SMALL);
if (objFunction==2)
{
   for(i=0;i<n;i++)
   {
      sumT=sumT+std::pow(T[i]/maxTSafe,Pnorm);
   }
   reduce(sumT, sumOp<scalar>());
   sumT=std::pow(sumT,1.0/Pnorm-1);
}

fvScalarMatrix TEqnb
(
   -fvm::div(-phi,Tb) 
   +fvm::laplacian(DT, Tb)
);
if (objFunction==1)
{
   TEqnb -= unitThermalSource;
}
if (objFunction==2)
{
   TEqnb -= b2*sumT*Foam::pow(T/maxTSafe,Pnorm-1);
}
if (objFunction==3 && static_cast<bool>(optThermalBoundedObjectiveEnable))
{
   volScalarField thermalExcessWeight
   (
      IOobject
      (
         "thermalExcessWeightTb",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("zeroWeightTb", dimless, 0.0)
   );
   const scalarField& TInternal = T.internalField();
   scalarField& weightInternal = thermalExcessWeight.primitiveFieldRef();
   const scalar thermalScale = Foam::max(thermalObjectiveSoftplusScale.value(), SMALL);
   const scalar thermalRef = thermalObjectiveTref.value();
   forAll(weightInternal, celli)
   {
      const scalar z = (TInternal[celli] - thermalRef)/thermalScale;
      if (z >= 0.0)
      {
         const scalar expNeg = std::exp(-z);
         weightInternal[celli] = 1.0/(1.0 + expNeg);
      }
      else
      {
         const scalar expPos = std::exp(z);
         weightInternal[celli] = expPos/(1.0 + expPos);
      }
   }
   thermalExcessWeight.correctBoundaryConditions();
   TEqnb -= unitThermalSource*thermalExcessWeight;
}
const auto TbLinearSolvePerformance = TEqnb.solve();
recordLinearSolve
(
  TbLinearSolveStats,
  TbLinearSolvePerformance.initialResidual(),
  TbLinearSolvePerformance.finalResidual(),
  TbLinearSolvePerformance.nIterations()
);
