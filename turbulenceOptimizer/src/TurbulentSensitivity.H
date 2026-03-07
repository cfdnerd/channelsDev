if (customTurbulenceActive && stageDTurbulentOptimizationMode)
{
  const bool stageDUseTransportedClosureAdjoints =
    static_cast<bool>(stageDTransportedClosureAdjoints)
    && static_cast<bool>(optDualConsistentStateMapEnforce);
  auto stageDGlobalAbsMax =
  [&](const volScalarField& field) -> scalar
  {
    scalar localMax = 0.0;
    const scalarField& values = field.internalField();
    forAll(values, celli)
    {
      localMax = Foam::max(localMax, Foam::mag(values[celli]));
    }
    reduce(localMax, maxOp<scalar>());
    return localMax;
  };

  const volSymmTensorField stageDSymmGradU(symm(fvc::grad(U)));
  const volSymmTensorField stageDSymmGradUa(symm(fvc::grad(Ua)));
  const volSymmTensorField stageDSymmGradUb(symm(fvc::grad(Ub)));

  const volScalarField stageDKObjectiveThermal
  (
    dAlphaEffDkBlock*(fvc::grad(T) & fvc::grad(Ta))
  );
  const volScalarField stageDKObjectiveMomentum
  (
    2.0*dNutDkBlock*(stageDSymmGradU && stageDSymmGradUb)
  );
  const volScalarField stageDOmegaObjectiveThermal
  (
    dAlphaEffDOmegaBlock*(fvc::grad(T) & fvc::grad(Ta))
  );
  const volScalarField stageDOmegaObjectiveMomentum
  (
    2.0*dNutDOmegaBlock*(stageDSymmGradU && stageDSymmGradUb)
  );
  const volScalarField stageDKConstraintMomentum
  (
    2.0*dNutDkBlock*(stageDSymmGradU && stageDSymmGradUa)
  );
  const volScalarField stageDOmegaConstraintMomentum
  (
    2.0*dNutDOmegaBlock*(stageDSymmGradU && stageDSymmGradUa)
  );
  const volScalarField stageDKObjectiveThermalMagnitude
  (
    mag(stageDKObjectiveThermal)
  );
  const volScalarField stageDKObjectiveMomentumMagnitude
  (
    mag(stageDKObjectiveMomentum)
  );
  const volScalarField stageDOmegaObjectiveThermalMagnitude
  (
    mag(stageDOmegaObjectiveThermal)
  );
  const volScalarField stageDOmegaObjectiveMomentumMagnitude
  (
    mag(stageDOmegaObjectiveMomentum)
  );
  const volScalarField stageDKConstraintMomentumMagnitude
  (
    mag(stageDKConstraintMomentum)
  );
  const volScalarField stageDOmegaConstraintMomentumMagnitude
  (
    mag(stageDOmegaConstraintMomentum)
  );

  const scalar stageDKObjectiveThermalScaleValue =
    Foam::max(stageDGlobalAbsMax(stageDKObjectiveThermalMagnitude), SMALL);
  const scalar stageDKObjectiveMomentumScaleValue =
    Foam::max(stageDGlobalAbsMax(stageDKObjectiveMomentumMagnitude), SMALL);
  const scalar stageDOmegaObjectiveThermalScaleValue =
    Foam::max(stageDGlobalAbsMax(stageDOmegaObjectiveThermalMagnitude), SMALL);
  const scalar stageDOmegaObjectiveMomentumScaleValue =
    Foam::max(stageDGlobalAbsMax(stageDOmegaObjectiveMomentumMagnitude), SMALL);
  const scalar stageDKConstraintMomentumScaleValue =
    Foam::max(stageDGlobalAbsMax(stageDKConstraintMomentumMagnitude), SMALL);
  const scalar stageDOmegaConstraintMomentumScaleValue =
    Foam::max(stageDGlobalAbsMax(stageDOmegaConstraintMomentumMagnitude), SMALL);

  const dimensionedScalar stageDKObjectiveThermalScale
  (
    "stageDKObjectiveThermalScale",
    stageDKObjectiveThermal.dimensions(),
    stageDKObjectiveThermalScaleValue
  );
  const dimensionedScalar stageDKObjectiveMomentumScale
  (
    "stageDKObjectiveMomentumScale",
    stageDKObjectiveMomentum.dimensions(),
    stageDKObjectiveMomentumScaleValue
  );
  const dimensionedScalar stageDOmegaObjectiveThermalScale
  (
    "stageDOmegaObjectiveThermalScale",
    stageDOmegaObjectiveThermal.dimensions(),
    stageDOmegaObjectiveThermalScaleValue
  );
  const dimensionedScalar stageDOmegaObjectiveMomentumScale
  (
    "stageDOmegaObjectiveMomentumScale",
    stageDOmegaObjectiveMomentum.dimensions(),
    stageDOmegaObjectiveMomentumScaleValue
  );
  const dimensionedScalar stageDKConstraintMomentumScale
  (
    "stageDKConstraintMomentumScale",
    stageDKConstraintMomentum.dimensions(),
    stageDKConstraintMomentumScaleValue
  );
  const dimensionedScalar stageDOmegaConstraintMomentumScale
  (
    "stageDOmegaConstraintMomentumScale",
    stageDOmegaConstraintMomentum.dimensions(),
    stageDOmegaConstraintMomentumScaleValue
  );

  const volScalarField stageDKObjectiveProxy
  (
      stageDKObjectiveThermal/stageDKObjectiveThermalScale
    + stageDKObjectiveMomentum/stageDKObjectiveMomentumScale
  );
  const volScalarField stageDOmegaObjectiveProxy
  (
      stageDOmegaObjectiveThermal/stageDOmegaObjectiveThermalScale
    + stageDOmegaObjectiveMomentum/stageDOmegaObjectiveMomentumScale
  );
  const volScalarField stageDKConstraintProxy
  (
    stageDKConstraintMomentum/stageDKConstraintMomentumScale
  );
  const volScalarField stageDOmegaConstraintProxy
  (
    stageDOmegaConstraintMomentum/stageDOmegaConstraintMomentumScale
  );
  const volScalarField stageDKMagnitude
  (
    mag
    (
      static_cast<bool>(stageDClosureConsistencyUseBoundedState)
      ? kClosureState
      : k
    )
  );
  const volScalarField stageDOmegaMagnitude
  (
    mag
    (
      static_cast<bool>(stageDClosureConsistencyUseBoundedState)
      ? omegaClosureState
      : omega
    )
  );

  const volScalarField stageDDampingDxh
  (
    turbulenceDampingCoeff*dAlphaDxh
  );
  const volScalarField stageDDampingMagnitude
  (
    mag(stageDDampingDxh)
  );
  const dimensionedScalar stageDKScale
  (
    "stageDKScale",
    k.dimensions(),
    Foam::max(stageDGlobalAbsMax(stageDKMagnitude), SMALL)
  );
  const dimensionedScalar stageDOmegaScale
  (
    "stageDOmegaScale",
    omega.dimensions(),
    Foam::max(stageDGlobalAbsMax(stageDOmegaMagnitude), SMALL)
  );
  const dimensionedScalar stageDDampingScale
  (
    "stageDDampingScale",
    stageDDampingDxh.dimensions(),
    Foam::max(stageDGlobalAbsMax(stageDDampingMagnitude), SMALL)
  );
  const volScalarField stageDDampingNormalized
  (
    stageDDampingDxh/stageDDampingScale
  );
  const volScalarField stageDKObjectiveNormalized
  (
    (
      (
        static_cast<bool>(stageDClosureConsistencyUseBoundedState)
        ? kClosureState
        : k
      )/stageDKScale
    )*stageDKObjectiveProxy
  );
  const volScalarField stageDOmegaObjectiveNormalized
  (
    (
      (
        static_cast<bool>(stageDClosureConsistencyUseBoundedState)
        ? omegaClosureState
        : omega
      )/stageDOmegaScale
    )*stageDOmegaObjectiveProxy
  );
  const volScalarField stageDKConstraintNormalized
  (
    (
      (
        static_cast<bool>(stageDClosureConsistencyUseBoundedState)
        ? kClosureState
        : k
      )/stageDKScale
    )*stageDKConstraintProxy
  );
  const volScalarField stageDOmegaConstraintNormalized
  (
    (
      (
        static_cast<bool>(stageDClosureConsistencyUseBoundedState)
        ? omegaClosureState
        : omega
      )/stageDOmegaScale
    )*stageDOmegaConstraintProxy
  );
  const volScalarField stageDObjectiveNormalizedMix
  (
    stageDKObjectiveNormalized + stageDOmegaObjectiveNormalized
  );
  const volScalarField stageDConstraintNormalizedMix
  (
    stageDKConstraintNormalized + stageDOmegaConstraintNormalized
  );
  const volScalarField stageDTransportedClosureMix
  (
    ka + omegaa
  );
  const volScalarField stageDLogClosureMix
  (
    logKClosureState + logOmegaClosureState
  );
  const volScalarField stageDLogClosureMixMagnitude
  (
    mag(stageDLogClosureMix)
  );
  const dimensionedScalar stageDLogClosureScale
  (
    "stageDLogClosureScale",
    stageDLogClosureMix.dimensions(),
    Foam::max(stageDGlobalAbsMax(stageDLogClosureMixMagnitude), SMALL)
  );
  const volScalarField stageDTransportedClosureMixMagnitude
  (
    mag(stageDTransportedClosureMix)
  );
  const dimensionedScalar stageDTransportedClosureScale
  (
    "stageDTransportedClosureScale",
    stageDTransportedClosureMix.dimensions(),
    Foam::max(stageDGlobalAbsMax(stageDTransportedClosureMixMagnitude), SMALL)
  );
  const volScalarField stageDTransportedClosureNormalized
  (
    stageDTransportedClosureMix/stageDTransportedClosureScale
  );
  volScalarField stageDObjectiveClosureMix
  (
    stageDUseTransportedClosureAdjoints
    ? stageDTransportedClosureNormalized
    : stageDObjectiveNormalizedMix
  );
  if
  (
      stageDUseTransportedClosureAdjoints
   && static_cast<bool>(stageDUseLogClosureState)
  )
  {
    stageDObjectiveClosureMix +=
      0.05*stageDLogClosureMix/stageDLogClosureScale;
  }
  const volScalarField stageDConstraintClosureMix
  (
    stageDUseTransportedClosureAdjoints
    ? stageDTransportedClosureNormalized
    : stageDConstraintNormalizedMix
  );
  const volScalarField stageDObjectiveClosureContribution
  (
    -fsenshMeanT*stageDDampingNormalized*stageDObjectiveClosureMix
  );
  const volScalarField stageDConstraintClosureContribution
  (
    -gsenshPowerDiss*stageDDampingNormalized*stageDConstraintClosureMix
  );

  fsenshMeanT += stageDObjectiveClosureContribution;
  gsenshPowerDiss += stageDConstraintClosureContribution;
  const volScalarField stageDObjectiveClosureMagnitude
  (
    mag(stageDObjectiveClosureContribution)
  );
  const volScalarField stageDConstraintClosureMagnitude
  (
    mag(stageDConstraintClosureContribution)
  );
  const volScalarField stageDKObjectiveProxyMagnitude
  (
    mag(stageDKObjectiveProxy)
  );
  const volScalarField stageDOmegaObjectiveProxyMagnitude
  (
    mag(stageDOmegaObjectiveProxy)
  );
  const volScalarField stageDKConstraintProxyMagnitude
  (
    mag(stageDKConstraintProxy)
  );
  const volScalarField stageDOmegaConstraintProxyMagnitude
  (
    mag(stageDOmegaConstraintProxy)
  );

  const scalar stageDObjectiveClosureMax =
    stageDGlobalAbsMax(stageDObjectiveClosureMagnitude);
  const scalar stageDConstraintClosureMax =
    stageDGlobalAbsMax(stageDConstraintClosureMagnitude);
  const scalar stageDKObjectiveProxyMax =
    stageDGlobalAbsMax(stageDKObjectiveProxyMagnitude);
  const scalar stageDOmegaObjectiveProxyMax =
    stageDGlobalAbsMax(stageDOmegaObjectiveProxyMagnitude);
  const scalar stageDKConstraintProxyMax =
    stageDGlobalAbsMax(stageDKConstraintProxyMagnitude);
  const scalar stageDOmegaConstraintProxyMax =
    stageDGlobalAbsMax(stageDOmegaConstraintProxyMagnitude);
  const scalar stageDTransportedClosureMixMax =
    stageDGlobalAbsMax(stageDTransportedClosureMixMagnitude);

  if (Pstream::master())
  {
    std::ofstream debugLogFile(debugOptimizerLogFile, std::ios::app);
    debugLogFile
      << "\n[Stage-D Turbulent Sensitivity]\n"
      << "  Obj closure |max|       : " << stageDObjectiveClosureMax << "\n"
      << "  Con closure |max|       : " << stageDConstraintClosureMax << "\n"
      << "  Obj k-proxy |max|       : " << stageDKObjectiveProxyMax << "\n"
      << "  Obj omega-proxy |max|   : " << stageDOmegaObjectiveProxyMax << "\n"
      << "  Con k-proxy |max|       : " << stageDKConstraintProxyMax << "\n"
      << "  Con omega-proxy |max|   : " << stageDOmegaConstraintProxyMax << "\n"
      << "  Closure mix mode        : "
      << (stageDUseTransportedClosureAdjoints ? "transported" : "proxy")
      << "\n"
      << "  Transported mix |max|   : " << stageDTransportedClosureMixMax << "\n";
  }
}
