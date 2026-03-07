auto computePatchAreaWeightedMean =
[&](const fvPatchScalarField& patchField) -> scalar
{
  const label patchi = patchField.patch().index();
  const scalarField& magSf = mesh.boundary()[patchi].magSf();
  scalar localWeightedSum = 0.0;
  scalar localAreaSum = 0.0;

  forAll(patchField, facei)
  {
    localWeightedSum += patchField[facei]*magSf[facei];
    localAreaSum += magSf[facei];
  }

  reduce(localWeightedSum, sumOp<scalar>());
  reduce(localAreaSum, sumOp<scalar>());

  return localWeightedSum/Foam::max(localAreaSum, SMALL);
};

auto updateTurbulenceResidualBlocks =
[&]()
{
  momentumResistanceBlock = alpha;
  momentumDiffusivityBlock = nuEff;
  thermalDiffusivityBlock = alphaEff;
  kDiffusionBlock = nu + alphaKField*nut;
  omegaDiffusionBlock = nu + alphaOmegaField*nut;
  kDestructionBlock = betaStarCoeff*omega;
  omegaProductionBlock = gammaField*PkField/Foam::max(nut, nutFloor);
  omegaDestructionBlock = betaField*omega;
  turbulenceDampingBlock = turbulenceDampingCoeff*alpha;
};

auto updateInletTurbulenceConsistency =
[&]()
{
  turbulenceInletKActual = turbulenceInletKTarget;
  turbulenceInletOmegaActual = turbulenceInletOmegaTarget;
  turbulenceInletKRelativeError = 0.0;
  turbulenceInletOmegaRelativeError = 0.0;
  turbulenceInletConsistencyStatus = "SKIP";

  if (!customTurbulenceActive || inletPatchId < 0)
  {
    turbulenceInletConsistent = true;
    return;
  }

  turbulenceInletKActual =
    computePatchAreaWeightedMean(k.boundaryField()[inletPatchId]);
  turbulenceInletOmegaActual =
    computePatchAreaWeightedMean(omega.boundaryField()[inletPatchId]);

  turbulenceInletKRelativeError =
    Foam::mag(turbulenceInletKActual - turbulenceInletKTarget)
    / Foam::max(Foam::mag(turbulenceInletKTarget), SMALL);
  turbulenceInletOmegaRelativeError =
    Foam::mag(turbulenceInletOmegaActual - turbulenceInletOmegaTarget)
    / Foam::max(Foam::mag(turbulenceInletOmegaTarget), SMALL);

  const scalar worstInletRelativeError =
    Foam::max(turbulenceInletKRelativeError, turbulenceInletOmegaRelativeError);

  if (worstInletRelativeError > 2.0*inletConsistencyRelTol)
  {
    turbulenceInletConsistencyStatus = "WARN";
    turbulenceInletConsistent = false;
  }
  else if (worstInletRelativeError > inletConsistencyRelTol)
  {
    turbulenceInletConsistencyStatus = "OBS";
    turbulenceInletConsistent = false;
  }
  else
  {
    turbulenceInletConsistencyStatus = "OK";
    turbulenceInletConsistent = true;
  }
};

auto updateTurbulenceStageBState =
[&]()
{
  updateTurbulenceResidualBlocks();
  updateInletTurbulenceConsistency();
};
