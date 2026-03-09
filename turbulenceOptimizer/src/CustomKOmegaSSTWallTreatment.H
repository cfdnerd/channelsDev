if (customTurbulenceActive)
{
  k.correctBoundaryConditions();
  omega.correctBoundaryConditions();

  const scalar cmuQuarter = Foam::pow(betaStarCoeff.value(), 0.25);
  const scalar sqrtBetaStar = Foam::sqrt(betaStarCoeff.value());
  const scalar kappa(0.41);
  const scalar yPlusSwitch(11.0);
  const scalar yPlusBlendWidth(2.0);

  forAll(mesh.boundary(), patchi)
  {
    if (mesh.boundary()[patchi].type() != "wall")
    {
      continue;
    }

    fvPatchScalarField& kPatch = k.boundaryFieldRef()[patchi];
    fvPatchScalarField& omegaPatch = omega.boundaryFieldRef()[patchi];
    const labelUList& faceCells = mesh.boundary()[patchi].faceCells();
    const scalarField& deltaCoeffs = mesh.boundary()[patchi].deltaCoeffs();

    forAll(kPatch, facei)
    {
      const scalar yFace =
        Foam::max
        (
          1.0/Foam::max(Foam::mag(deltaCoeffs[facei]), SMALL),
          yMin.value()
        );

      if (!wallTreatmentModeWallFunctionCompatible)
      {
        const scalar omegaWallValue =
          omegaWallFactor*nu.value()/(beta1Coeff.value()*sqr(yFace));
        if (!(omegaWallValue > 0.0) || omegaWallValue >= GREAT)
        {
          FatalErrorInFunction
            << "Non-physical wall omega generated on patch "
            << mesh.boundary()[patchi].name()
            << " face " << facei
            << ": omegaWall=" << omegaWallValue
            << ". Clipping is disabled by policy; inspect wall-treatment inputs."
            << exit(FatalError);
        }

        kPatch[facei] = kMin.value();
        omegaPatch[facei] = omegaWallValue;
      }
      else
      {
        const label celli = faceCells[facei];
        const scalar kWallValue = Foam::max(k[celli], kMin.value());
        const scalar uTau = cmuQuarter*Foam::sqrt(kWallValue);
        const scalar yPlus = uTau*yFace/Foam::max(nu.value(), SMALL);
        const scalar omegaViscous =
          omegaWallFactor*nu.value()/(beta1Coeff.value()*sqr(yFace));
        const scalar omegaLog =
          uTau/Foam::max(sqrtBetaStar*kappa*yFace, SMALL);
        const scalar blend =
          0.5*(1.0 + Foam::tanh((yPlus - yPlusSwitch)/yPlusBlendWidth));
        const scalar omegaWallValue =
          Foam::max((1.0 - blend)*omegaViscous + blend*omegaLog, omegaMin.value());

        if (!(omegaWallValue > 0.0) || omegaWallValue >= GREAT)
        {
          FatalErrorInFunction
            << "Non-physical wall omega generated on patch "
            << mesh.boundary()[patchi].name()
            << " face " << facei
            << ": omegaWall=" << omegaWallValue
            << ". Clipping is disabled by policy; inspect wall-treatment inputs."
            << exit(FatalError);
        }

        kPatch[facei] = kWallValue;
        omegaPatch[facei] = omegaWallValue;
      }
    }
  }
}
else
{
  k.correctBoundaryConditions();
  omega.correctBoundaryConditions();
}
