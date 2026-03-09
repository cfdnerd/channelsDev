if (customTurbulenceActive && optimizationCycleTurbulentOptimizationMode)
{
  const bool optimizationCycleUseTransportedClosureAdjoints =
    static_cast<bool>(optimizationCycleTransportedClosureAdjoints)
    && static_cast<bool>(optDualConsistentStateMapEnforce);
  auto optimizationCycleGlobalAbsMax =
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

  const volSymmTensorField optimizationCycleSymmGradU(symm(fvc::grad(U)));
  const volSymmTensorField optimizationCycleSymmGradUa(symm(fvc::grad(Ua)));
  const volSymmTensorField optimizationCycleSymmGradUb(symm(fvc::grad(Ub)));

  const volScalarField optimizationCycleKObjectiveThermal
  (
    dAlphaEffDkBlock*(fvc::grad(T) & fvc::grad(optimizationCycleThermalAdjointField))
  );
  const volScalarField optimizationCycleKObjectiveMomentum
  (
    2.0*dNutDkBlock*(optimizationCycleSymmGradU && optimizationCycleSymmGradUb)
  );
  const volScalarField optimizationCycleOmegaObjectiveThermal
  (
    dAlphaEffDOmegaBlock*(fvc::grad(T) & fvc::grad(optimizationCycleThermalAdjointField))
  );
  const volScalarField optimizationCycleOmegaObjectiveMomentum
  (
    2.0*dNutDOmegaBlock*(optimizationCycleSymmGradU && optimizationCycleSymmGradUb)
  );
  const volScalarField optimizationCycleKConstraintMomentum
  (
    2.0*dNutDkBlock*(optimizationCycleSymmGradU && optimizationCycleSymmGradUa)
  );
  const volScalarField optimizationCycleOmegaConstraintMomentum
  (
    2.0*dNutDOmegaBlock*(optimizationCycleSymmGradU && optimizationCycleSymmGradUa)
  );
  const volScalarField optimizationCycleKObjectiveThermalMagnitude
  (
    mag(optimizationCycleKObjectiveThermal)
  );
  const volScalarField optimizationCycleKObjectiveMomentumMagnitude
  (
    mag(optimizationCycleKObjectiveMomentum)
  );
  const volScalarField optimizationCycleOmegaObjectiveThermalMagnitude
  (
    mag(optimizationCycleOmegaObjectiveThermal)
  );
  const volScalarField optimizationCycleOmegaObjectiveMomentumMagnitude
  (
    mag(optimizationCycleOmegaObjectiveMomentum)
  );
  const volScalarField optimizationCycleKConstraintMomentumMagnitude
  (
    mag(optimizationCycleKConstraintMomentum)
  );
  const volScalarField optimizationCycleOmegaConstraintMomentumMagnitude
  (
    mag(optimizationCycleOmegaConstraintMomentum)
  );

  const scalar optimizationCycleKObjectiveThermalScaleValue =
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleKObjectiveThermalMagnitude), SMALL);
  const scalar optimizationCycleKObjectiveMomentumScaleValue =
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleKObjectiveMomentumMagnitude), SMALL);
  const scalar optimizationCycleOmegaObjectiveThermalScaleValue =
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleOmegaObjectiveThermalMagnitude), SMALL);
  const scalar optimizationCycleOmegaObjectiveMomentumScaleValue =
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleOmegaObjectiveMomentumMagnitude), SMALL);
  const scalar optimizationCycleKConstraintMomentumScaleValue =
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleKConstraintMomentumMagnitude), SMALL);
  const scalar optimizationCycleOmegaConstraintMomentumScaleValue =
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleOmegaConstraintMomentumMagnitude), SMALL);

  const dimensionedScalar optimizationCycleKObjectiveThermalScale
  (
    "optimizationCycleKObjectiveThermalScale",
    optimizationCycleKObjectiveThermal.dimensions(),
    optimizationCycleKObjectiveThermalScaleValue
  );
  const dimensionedScalar optimizationCycleKObjectiveMomentumScale
  (
    "optimizationCycleKObjectiveMomentumScale",
    optimizationCycleKObjectiveMomentum.dimensions(),
    optimizationCycleKObjectiveMomentumScaleValue
  );
  const dimensionedScalar optimizationCycleOmegaObjectiveThermalScale
  (
    "optimizationCycleOmegaObjectiveThermalScale",
    optimizationCycleOmegaObjectiveThermal.dimensions(),
    optimizationCycleOmegaObjectiveThermalScaleValue
  );
  const dimensionedScalar optimizationCycleOmegaObjectiveMomentumScale
  (
    "optimizationCycleOmegaObjectiveMomentumScale",
    optimizationCycleOmegaObjectiveMomentum.dimensions(),
    optimizationCycleOmegaObjectiveMomentumScaleValue
  );
  const dimensionedScalar optimizationCycleKConstraintMomentumScale
  (
    "optimizationCycleKConstraintMomentumScale",
    optimizationCycleKConstraintMomentum.dimensions(),
    optimizationCycleKConstraintMomentumScaleValue
  );
  const dimensionedScalar optimizationCycleOmegaConstraintMomentumScale
  (
    "optimizationCycleOmegaConstraintMomentumScale",
    optimizationCycleOmegaConstraintMomentum.dimensions(),
    optimizationCycleOmegaConstraintMomentumScaleValue
  );

  const volScalarField optimizationCycleKObjectiveProxy
  (
      optimizationCycleKObjectiveThermal/optimizationCycleKObjectiveThermalScale
    + optimizationCycleKObjectiveMomentum/optimizationCycleKObjectiveMomentumScale
  );
  const volScalarField optimizationCycleOmegaObjectiveProxy
  (
      optimizationCycleOmegaObjectiveThermal/optimizationCycleOmegaObjectiveThermalScale
    + optimizationCycleOmegaObjectiveMomentum/optimizationCycleOmegaObjectiveMomentumScale
  );
  const volScalarField optimizationCycleKConstraintProxy
  (
    optimizationCycleKConstraintMomentum/optimizationCycleKConstraintMomentumScale
  );
  const volScalarField optimizationCycleOmegaConstraintProxy
  (
    optimizationCycleOmegaConstraintMomentum/optimizationCycleOmegaConstraintMomentumScale
  );
  const volScalarField optimizationCycleKMagnitude
  (
    mag
    (
      static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
      ? kClosureState
      : k
    )
  );
  const volScalarField optimizationCycleOmegaMagnitude
  (
    mag
    (
      static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
      ? omegaClosureState
      : omega
    )
  );

  const volScalarField optimizationCycleDampingDxh
  (
    turbulenceDampingCoeff*dAlphaDxh
  );
  const volScalarField optimizationCycleDampingMagnitude
  (
    mag(optimizationCycleDampingDxh)
  );
  const dimensionedScalar optimizationCycleKScale
  (
    "optimizationCycleKScale",
    k.dimensions(),
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleKMagnitude), SMALL)
  );
  const dimensionedScalar optimizationCycleOmegaScale
  (
    "optimizationCycleOmegaScale",
    omega.dimensions(),
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleOmegaMagnitude), SMALL)
  );
  const dimensionedScalar optimizationCycleDampingScale
  (
    "optimizationCycleDampingScale",
    optimizationCycleDampingDxh.dimensions(),
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleDampingMagnitude), SMALL)
  );
  const volScalarField optimizationCycleDampingNormalized
  (
    optimizationCycleDampingDxh/optimizationCycleDampingScale
  );
  const volScalarField optimizationCycleKObjectiveNormalized
  (
    (
      (
        static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
        ? kClosureState
        : k
      )/optimizationCycleKScale
    )*optimizationCycleKObjectiveProxy
  );
  const volScalarField optimizationCycleOmegaObjectiveNormalized
  (
    (
      (
        static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
        ? omegaClosureState
        : omega
      )/optimizationCycleOmegaScale
    )*optimizationCycleOmegaObjectiveProxy
  );
  const volScalarField optimizationCycleKConstraintNormalized
  (
    (
      (
        static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
        ? kClosureState
        : k
      )/optimizationCycleKScale
    )*optimizationCycleKConstraintProxy
  );
  const volScalarField optimizationCycleOmegaConstraintNormalized
  (
    (
      (
        static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
        ? omegaClosureState
        : omega
      )/optimizationCycleOmegaScale
    )*optimizationCycleOmegaConstraintProxy
  );
  const volScalarField optimizationCycleObjectiveNormalizedMix
  (
    optimizationCycleKObjectiveNormalized + optimizationCycleOmegaObjectiveNormalized
  );
  const volScalarField optimizationCycleConstraintNormalizedMix
  (
    optimizationCycleKConstraintNormalized + optimizationCycleOmegaConstraintNormalized
  );
  const bool optimizationCycleUseSAClosureAdjoint =
    turbulenceProModeActive
    && proModeClosureAdjointSolvesActive
    && proModeSAAdjointClosureActive;
  volScalarField optimizationCycleTransportedClosureMix
  (
    IOobject
    (
      "optimizationCycleTransportedClosureMix",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
      "zeroOptimizationCycleTransportedClosureMix",
      ka.dimensions(),
      0.0
    )
  );
  optimizationCycleTransportedClosureMix = ka;
  if (!optimizationCycleUseSAClosureAdjoint)
  {
    optimizationCycleTransportedClosureMix += omegaa;
  }
  const volScalarField optimizationCycleLogClosureMix
  (
    logKClosureState + logOmegaClosureState
  );
  const volScalarField optimizationCycleLogClosureMixMagnitude
  (
    mag(optimizationCycleLogClosureMix)
  );
  const dimensionedScalar optimizationCycleLogClosureScale
  (
    "optimizationCycleLogClosureScale",
    optimizationCycleLogClosureMix.dimensions(),
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleLogClosureMixMagnitude), SMALL)
  );
  const volScalarField optimizationCycleTransportedClosureMixMagnitude
  (
    mag(optimizationCycleTransportedClosureMix)
  );
  const dimensionedScalar optimizationCycleTransportedClosureScale
  (
    "optimizationCycleTransportedClosureScale",
    optimizationCycleTransportedClosureMix.dimensions(),
    Foam::max(optimizationCycleGlobalAbsMax(optimizationCycleTransportedClosureMixMagnitude), SMALL)
  );
  const volScalarField optimizationCycleTransportedClosureNormalized
  (
    optimizationCycleTransportedClosureMix/optimizationCycleTransportedClosureScale
  );
  volScalarField optimizationCycleObjectiveClosureMix
  (
    optimizationCycleUseTransportedClosureAdjoints
    ? optimizationCycleTransportedClosureNormalized
    : optimizationCycleObjectiveNormalizedMix
  );
  if
  (
      optimizationCycleUseTransportedClosureAdjoints
   && static_cast<bool>(optimizationCycleUseLogClosureState)
  )
  {
    optimizationCycleObjectiveClosureMix +=
      0.05*optimizationCycleLogClosureMix/optimizationCycleLogClosureScale;
  }
  const volScalarField optimizationCycleConstraintClosureMix
  (
    optimizationCycleUseTransportedClosureAdjoints
    ? optimizationCycleTransportedClosureNormalized
    : optimizationCycleConstraintNormalizedMix
  );
  const volScalarField optimizationCycleObjectiveClosureContribution
  (
    -fsenshMeanT*optimizationCycleDampingNormalized*optimizationCycleObjectiveClosureMix
  );
  const volScalarField optimizationCycleConstraintClosureContribution
  (
    -gsenshPowerDiss*optimizationCycleDampingNormalized*optimizationCycleConstraintClosureMix
  );

  fsenshMeanT += optimizationCycleObjectiveClosureContribution;
  gsenshPowerDiss += optimizationCycleConstraintClosureContribution;
  const volScalarField optimizationCycleObjectiveClosureMagnitude
  (
    mag(optimizationCycleObjectiveClosureContribution)
  );
  const volScalarField optimizationCycleConstraintClosureMagnitude
  (
    mag(optimizationCycleConstraintClosureContribution)
  );
  const volScalarField optimizationCycleKObjectiveProxyMagnitude
  (
    mag(optimizationCycleKObjectiveProxy)
  );
  const volScalarField optimizationCycleOmegaObjectiveProxyMagnitude
  (
    mag(optimizationCycleOmegaObjectiveProxy)
  );
  const volScalarField optimizationCycleKConstraintProxyMagnitude
  (
    mag(optimizationCycleKConstraintProxy)
  );
  const volScalarField optimizationCycleOmegaConstraintProxyMagnitude
  (
    mag(optimizationCycleOmegaConstraintProxy)
  );

  const scalar optimizationCycleObjectiveClosureMax =
    optimizationCycleGlobalAbsMax(optimizationCycleObjectiveClosureMagnitude);
  const scalar optimizationCycleConstraintClosureMax =
    optimizationCycleGlobalAbsMax(optimizationCycleConstraintClosureMagnitude);
  const scalar optimizationCycleKObjectiveProxyMax =
    optimizationCycleGlobalAbsMax(optimizationCycleKObjectiveProxyMagnitude);
  const scalar optimizationCycleOmegaObjectiveProxyMax =
    optimizationCycleGlobalAbsMax(optimizationCycleOmegaObjectiveProxyMagnitude);
  const scalar optimizationCycleKConstraintProxyMax =
    optimizationCycleGlobalAbsMax(optimizationCycleKConstraintProxyMagnitude);
  const scalar optimizationCycleOmegaConstraintProxyMax =
    optimizationCycleGlobalAbsMax(optimizationCycleOmegaConstraintProxyMagnitude);
  const scalar optimizationCycleTransportedClosureMixMax =
    optimizationCycleGlobalAbsMax(optimizationCycleTransportedClosureMixMagnitude);

  if (Pstream::master())
  {
    std::ofstream debugLogFile(debugOptimizerLogFile, std::ios::app);
    debugLogFile
      << "\n[Optimization-Cycle Turbulent Sensitivity]\n"
      << "  Obj closure |max|       : " << optimizationCycleObjectiveClosureMax << "\n"
      << "  Con closure |max|       : " << optimizationCycleConstraintClosureMax << "\n"
      << "  Obj k-proxy |max|       : " << optimizationCycleKObjectiveProxyMax << "\n"
      << "  Obj omega-proxy |max|   : " << optimizationCycleOmegaObjectiveProxyMax << "\n"
      << "  Con k-proxy |max|       : " << optimizationCycleKConstraintProxyMax << "\n"
      << "  Con omega-proxy |max|   : " << optimizationCycleOmegaConstraintProxyMax << "\n"
      << "  Closure mix mode        : "
      << (optimizationCycleUseTransportedClosureAdjoints ? "transported" : "proxy")
      << "\n"
      << "  Transported mix |max|   : " << optimizationCycleTransportedClosureMixMax << "\n";
  }
}
