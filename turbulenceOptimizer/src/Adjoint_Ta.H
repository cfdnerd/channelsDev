sumT = 0;
maxT = gMax(T);
const scalar maxTAdjointScale = Foam::max(maxT, SMALL);
for (i = 0; i < n; i++)
{
   sumT = sumT + std::pow(T[i]/maxTAdjointScale, Pnorm);
}
reduce(sumT, sumOp<scalar>());
sumT = std::pow(sumT, 1.0/Pnorm - 1);

fvScalarMatrix TaEqn
(
   -fvm::div(-phi, Ta)
   + fvm::laplacian(alphaEff, Ta)
);
if (objFunction == 1)
{
   TaEqn -= unitThermalSource;
}
if (objFunction == 2)
{
   TaEqn -= b2*sumT*Foam::pow(T/maxTAdjointScale, Pnorm - 1);
}
if (objFunction == 3 && static_cast<bool>(optThermalBoundedObjectiveEnable))
{
   volScalarField thermalExcessWeight
   (
      IOobject
      (
         "thermalExcessWeightTa",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("zeroWeightTa", dimless, 0.0)
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
   TaEqn -= unitThermalSource*thermalExcessWeight;
}

// Stage-C thermal adjoint is currently a standalone verification block.
// Explicitly pin one cell to remove the null mode seen in parallel runs.
if (mesh.nCells() > 0)
{
   TaEqn.setReference(0, 0.0);
}

const auto TaLinearSolvePerformance = TaEqn.solve();
recordLinearSolve
(
  TaLinearSolveStats,
  TaLinearSolvePerformance.initialResidual(),
  TaLinearSolvePerformance.finalResidual(),
  TaLinearSolvePerformance.nIterations()
);
