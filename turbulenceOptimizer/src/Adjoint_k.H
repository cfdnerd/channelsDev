if (customTurbulenceActive)
{
  const bool optimizationCycleUseTransportedClosureAdjoints =
    turbulenceAdjointClosureMode
    && static_cast<bool>(optimizationCycleTransportedClosureAdjoints)
    && static_cast<bool>(optDualConsistentStateMapEnforce);
  const bool optimizationCycleUseResidualDraftClosureAdjoints =
    optimizationCycleUseTransportedClosureAdjoints
    && optimizationCycleClosureAdjointFormulation == word("transportedResidualDraft");
  const bool optimizationCycleUseCrossCoupledClosureAdjoints =
    optimizationCycleUseTransportedClosureAdjoints
    && optimizationCycleClosureAdjointFormulation == word("transportedCrossCoupledV3");
  auto applyDualConsistentKaBC =
  [&]()
  {
    if (optimizationCycleUseTransportedClosureAdjoints && static_cast<bool>(optimizationCycleClosureAdjointDualConsistentBC))
    {
      forAll(ka.boundaryFieldRef(), patchi)
      {
        if (k.boundaryField()[patchi].fixesValue())
        {
          ka.boundaryFieldRef()[patchi] == scalar(0.0);
        }
      }
    }
  };

  const volSymmTensorField adjointCycleSymmGradU(symm(fvc::grad(U)));
  const volSymmTensorField adjointCycleSymmGradUa(symm(fvc::grad(Ua)));
  const volScalarField& optimizationCycleOmegaStateRef =
    static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
    ? omegaClosureState
    : omega;
  const volScalarField adjointCycleKThermalCoupling
  (
    dAlphaEffDkBlock*(fvc::grad(T) & fvc::grad(optimizationCycleThermalAdjointField))
  );
  const volScalarField adjointCycleKMomentumCoupling
  (
    2.0*dNutDkBlock*(adjointCycleSymmGradU && adjointCycleSymmGradUa)
  );
  const volScalarField adjointCycleKThermalMagnitude(mag(adjointCycleKThermalCoupling));
  const volScalarField adjointCycleKMomentumMagnitude(mag(adjointCycleKMomentumCoupling));
  const scalar adjointCycleKThermalScaleValue =
    Foam::max(gMax(adjointCycleKThermalMagnitude), SMALL);
  const scalar adjointCycleKMomentumScaleValue =
    Foam::max(gMax(adjointCycleKMomentumMagnitude), SMALL);
  const dimensionedScalar adjointCycleKThermalScale
  (
    "adjointCycleKThermalScale",
    adjointCycleKThermalCoupling.dimensions(),
    adjointCycleKThermalScaleValue
  );
  const dimensionedScalar adjointCycleKMomentumScale
  (
    "adjointCycleKMomentumScale",
    adjointCycleKMomentumCoupling.dimensions(),
    adjointCycleKMomentumScaleValue
  );
  volScalarField adjointCycleKSourceProxy
  (
    adjointCycleKThermalCoupling/adjointCycleKThermalScale
  + adjointCycleKMomentumCoupling/adjointCycleKMomentumScale
  );
  if (static_cast<bool>(optimizationCycleUseLogClosureState))
  {
    const volScalarField logKMagnitude(mag(logKClosureState));
    const dimensionedScalar logKScale
    (
      "logKScale",
      logKClosureState.dimensions(),
      Foam::max(gMax(logKMagnitude), SMALL)
    );
    adjointCycleKSourceProxy += 0.05*logKClosureState/logKScale;
  }
  if (optimizationCycleUseTransportedClosureAdjoints)
  {
    applyDualConsistentKaBC();

    const dimensionedScalar optimizationCycleClosureInvTau
    (
      "optimizationCycleClosureInvTau",
      dimless/dimTime,
      optimizationCycleClosureAdjointInvTau
    );
    volScalarField optimizationCycleKReactionCoeff
    (
      betaStarCoeff
     *optimizationCycleOmegaStateRef
    + turbulenceDampingCoeff*alpha
    );
    if (!optimizationCycleUseResidualDraftClosureAdjoints)
    {
      optimizationCycleKReactionCoeff =
        dimensionedScalar
        (
          "zeroOptimizationCycleKReactionCoeff",
          optimizationCycleKReactionCoeff.dimensions(),
          0.0
        );
    }
    volScalarField optimizationCycleKCrossSource
    (
      IOobject
      (
        "optimizationCycleKCrossSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("zeroOptimizationCycleKCrossSource", dimless/dimTime, 0.0)
    );
    if (optimizationCycleUseCrossCoupledClosureAdjoints)
    {
      const volScalarField omegaaMagnitude(mag(omegaa));
      const dimensionedScalar omegaaScale
      (
        "omegaaScale",
        omegaa.dimensions(),
        Foam::max(gMax(omegaaMagnitude), SMALL)
      );
      optimizationCycleKCrossSource =
        optimizationCycleClosureInvTau
       *optimizationCycleClosureCrossCoupling
       *(omegaa/omegaaScale);
    }
    const volScalarField kaTransportDiffusivity
    (
      turbulenceAdjointTransportDiffusionScale*(nu + alphaKField*nut)
    );
    fvScalarMatrix kaEqn
    (
        fvm::div(phi, ka)
      - fvm::laplacian(kaTransportDiffusivity, ka)
      + fvm::Sp(optimizationCycleClosureInvTau, ka)
      + fvm::Sp(optimizationCycleKReactionCoeff, ka)
     ==
        optimizationCycleClosureInvTau*adjointCycleKSourceProxy
      + optimizationCycleKCrossSource
    );

    kaEqn.relax();

    const auto kaLinearSolvePerformance = kaEqn.solve();
    recordLinearSolve
    (
      kaLinearSolveStats,
      kaLinearSolvePerformance.initialResidual(),
      kaLinearSolvePerformance.finalResidual(),
      kaLinearSolvePerformance.nIterations()
    );

    applyDualConsistentKaBC();
  }
  else
  {
    ka = adjointCycleKSourceProxy;
  }

  const volScalarField kaMagnitudeCheck(mag(ka));
  const scalar kaAbsMax = gMax(kaMagnitudeCheck);
  if (!(kaAbsMax < optimizationCycleAdjointMaxMagnitude))
  {
    FatalErrorInFunction
      << "Non-physical ka magnitude detected: max(|ka|)=" << kaAbsMax
      << ". Clipping is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  ka.correctBoundaryConditions();
  applyDualConsistentKaBC();
}
