if (customTurbulenceActive)
{
  k.correctBoundaryConditions();
  omega.correctBoundaryConditions();

  forAll(mesh.boundary(), patchi)
  {
    if (mesh.boundary()[patchi].type() != "wall")
    {
      continue;
    }

    fvPatchScalarField& kPatch = k.boundaryFieldRef()[patchi];
    fvPatchScalarField& omegaPatch = omega.boundaryFieldRef()[patchi];
    const scalarField& deltaCoeffs = mesh.boundary()[patchi].deltaCoeffs();

    forAll(kPatch, facei)
    {
      const scalar yFace = Foam::max(1.0/Foam::max(deltaCoeffs[facei], SMALL), yMin.value());
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
  }
}
else
{
  k.correctBoundaryConditions();
  omega.correctBoundaryConditions();
}
