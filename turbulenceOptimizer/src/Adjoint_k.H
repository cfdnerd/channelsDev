if (customTurbulenceActive)
{
  const bool stageDUseTransportedClosureAdjoints =
    stageDTurbulentOptimizationMode
    && static_cast<bool>(stageDTransportedClosureAdjoints)
    && static_cast<bool>(optDualConsistentStateMapEnforce);
  const bool stageDUseResidualDraftClosureAdjoints =
    stageDUseTransportedClosureAdjoints
    && stageDClosureAdjointFormulation == word("transportedResidualDraft");
  const bool stageDUseCrossCoupledClosureAdjoints =
    stageDUseTransportedClosureAdjoints
    && stageDClosureAdjointFormulation == word("transportedCrossCoupledV3");
  auto applyDualConsistentKaBC =
  [&]()
  {
    if (stageDUseTransportedClosureAdjoints && static_cast<bool>(stageDClosureAdjointDualConsistentBC))
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

  const volSymmTensorField stageCSymmGradU(symm(fvc::grad(U)));
  const volSymmTensorField stageCSymmGradUa(symm(fvc::grad(Ua)));
  const volScalarField& stageDOmegaStateRef =
    static_cast<bool>(stageDClosureConsistencyUseBoundedState)
    ? omegaClosureState
    : omega;
  const volScalarField stageCKThermalCoupling
  (
    dAlphaEffDkBlock*(fvc::grad(T) & fvc::grad(Ta))
  );
  const volScalarField stageCKMomentumCoupling
  (
    2.0*dNutDkBlock*(stageCSymmGradU && stageCSymmGradUa)
  );
  const volScalarField stageCKThermalMagnitude(mag(stageCKThermalCoupling));
  const volScalarField stageCKMomentumMagnitude(mag(stageCKMomentumCoupling));
  const scalar stageCKThermalScaleValue =
    Foam::max(gMax(stageCKThermalMagnitude), SMALL);
  const scalar stageCKMomentumScaleValue =
    Foam::max(gMax(stageCKMomentumMagnitude), SMALL);
  const dimensionedScalar stageCKThermalScale
  (
    "stageCKThermalScale",
    stageCKThermalCoupling.dimensions(),
    stageCKThermalScaleValue
  );
  const dimensionedScalar stageCKMomentumScale
  (
    "stageCKMomentumScale",
    stageCKMomentumCoupling.dimensions(),
    stageCKMomentumScaleValue
  );
  volScalarField stageCKSourceProxy
  (
    stageCKThermalCoupling/stageCKThermalScale
  + stageCKMomentumCoupling/stageCKMomentumScale
  );
  if (static_cast<bool>(stageDUseLogClosureState))
  {
    const volScalarField logKMagnitude(mag(logKClosureState));
    const dimensionedScalar logKScale
    (
      "logKScale",
      logKClosureState.dimensions(),
      Foam::max(gMax(logKMagnitude), SMALL)
    );
    stageCKSourceProxy += 0.05*logKClosureState/logKScale;
  }
  if (stageDUseTransportedClosureAdjoints)
  {
    applyDualConsistentKaBC();

    const dimensionedScalar stageDClosureInvTau
    (
      "stageDClosureInvTau",
      dimless/dimTime,
      stageDClosureAdjointInvTau
    );
    volScalarField stageDKReactionCoeff
    (
      betaStarCoeff
     *stageDOmegaStateRef
    + turbulenceDampingCoeff*alpha
    );
    if (!stageDUseResidualDraftClosureAdjoints)
    {
      stageDKReactionCoeff =
        dimensionedScalar
        (
          "zeroStageDKReactionCoeff",
          stageDKReactionCoeff.dimensions(),
          0.0
        );
    }
    volScalarField stageDKCrossSource
    (
      IOobject
      (
        "stageDKCrossSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("zeroStageDKCrossSource", dimless/dimTime, 0.0)
    );
    if (stageDUseCrossCoupledClosureAdjoints)
    {
      const volScalarField omegaaMagnitude(mag(omegaa));
      const dimensionedScalar omegaaScale
      (
        "omegaaScale",
        omegaa.dimensions(),
        Foam::max(gMax(omegaaMagnitude), SMALL)
      );
      stageDKCrossSource =
        stageDClosureInvTau
       *stageDClosureCrossCoupling
       *(omegaa/omegaaScale);
    }
    fvScalarMatrix kaEqn
    (
        fvm::div(phi, ka)
      - fvm::laplacian(nu + alphaKField*nut, ka)
      + fvm::Sp(stageDClosureInvTau, ka)
      + fvm::Sp(stageDKReactionCoeff, ka)
     ==
        stageDClosureInvTau*stageCKSourceProxy
      + stageDKCrossSource
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
    ka = stageCKSourceProxy;
  }

  const volScalarField kaMagnitudeCheck(mag(ka));
  const scalar kaAbsMax = gMax(kaMagnitudeCheck);
  if (!(kaAbsMax < stageDAdjointMaxMagnitude))
  {
    FatalErrorInFunction
      << "Non-physical ka magnitude detected: max(|ka|)=" << kaAbsMax
      << ". Clipping is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  ka.correctBoundaryConditions();
  applyDualConsistentKaBC();
}
