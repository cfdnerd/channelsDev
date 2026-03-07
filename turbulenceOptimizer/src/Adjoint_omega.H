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
  auto applyDualConsistentOmegaaBC =
  [&]()
  {
    if (stageDUseTransportedClosureAdjoints && static_cast<bool>(stageDClosureAdjointDualConsistentBC))
    {
      forAll(omegaa.boundaryFieldRef(), patchi)
      {
        if (omega.boundaryField()[patchi].fixesValue())
        {
          omegaa.boundaryFieldRef()[patchi] == scalar(0.0);
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
  const volScalarField stageCOmegaThermalCoupling
  (
    dAlphaEffDOmegaBlock*(fvc::grad(T) & fvc::grad(Ta))
  );
  const volScalarField stageCOmegaMomentumCoupling
  (
    2.0*dNutDOmegaBlock*(stageCSymmGradU && stageCSymmGradUa)
  );
  const volScalarField stageCOmegaThermalMagnitude(mag(stageCOmegaThermalCoupling));
  const volScalarField stageCOmegaMomentumMagnitude(mag(stageCOmegaMomentumCoupling));
  const scalar stageCOmegaThermalScaleValue =
    Foam::max(gMax(stageCOmegaThermalMagnitude), SMALL);
  const scalar stageCOmegaMomentumScaleValue =
    Foam::max(gMax(stageCOmegaMomentumMagnitude), SMALL);
  const dimensionedScalar stageCOmegaThermalScale
  (
    "stageCOmegaThermalScale",
    stageCOmegaThermalCoupling.dimensions(),
    stageCOmegaThermalScaleValue
  );
  const dimensionedScalar stageCOmegaMomentumScale
  (
    "stageCOmegaMomentumScale",
    stageCOmegaMomentumCoupling.dimensions(),
    stageCOmegaMomentumScaleValue
  );
  volScalarField stageCOmegaSourceProxy
  (
    stageCOmegaThermalCoupling/stageCOmegaThermalScale
  + stageCOmegaMomentumCoupling/stageCOmegaMomentumScale
  );
  if (static_cast<bool>(stageDUseLogClosureState))
  {
    const volScalarField logOmegaMagnitude(mag(logOmegaClosureState));
    const dimensionedScalar logOmegaScale
    (
      "logOmegaScale",
      logOmegaClosureState.dimensions(),
      Foam::max(gMax(logOmegaMagnitude), SMALL)
    );
    stageCOmegaSourceProxy += 0.05*logOmegaClosureState/logOmegaScale;
  }
  if (stageDUseTransportedClosureAdjoints)
  {
    applyDualConsistentOmegaaBC();

    const dimensionedScalar stageDClosureInvTau
    (
      "stageDClosureInvTau",
      dimless/dimTime,
      stageDClosureAdjointInvTau
    );
    volScalarField stageDOmegaReactionCoeff
    (
      scalar(2.0)
     *betaField
     *stageDOmegaStateRef
    + turbulenceDampingCoeff*alpha
    );
    if (!stageDUseResidualDraftClosureAdjoints)
    {
      stageDOmegaReactionCoeff =
        dimensionedScalar
        (
          "zeroStageDOmegaReactionCoeff",
          stageDOmegaReactionCoeff.dimensions(),
          0.0
        );
    }
    volScalarField stageDOmegaCrossSource
    (
      IOobject
      (
        "stageDOmegaCrossSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("zeroStageDOmegaCrossSource", dimless/dimTime, 0.0)
    );
    if (stageDUseCrossCoupledClosureAdjoints)
    {
      const volScalarField kaMagnitude(mag(ka));
      const dimensionedScalar kaScale
      (
        "kaScale",
        ka.dimensions(),
        Foam::max(gMax(kaMagnitude), SMALL)
      );
      stageDOmegaCrossSource =
        stageDClosureInvTau
       *stageDClosureCrossCoupling
       *(ka/kaScale);
    }
    fvScalarMatrix omegaaEqn
    (
        fvm::div(phi, omegaa)
      - fvm::laplacian(nu + alphaOmegaField*nut, omegaa)
      + fvm::Sp(stageDClosureInvTau, omegaa)
      + fvm::Sp(stageDOmegaReactionCoeff, omegaa)
     ==
        stageDClosureInvTau*stageCOmegaSourceProxy
      + stageDOmegaCrossSource
    );

    omegaaEqn.relax();

    const auto omegaaLinearSolvePerformance = omegaaEqn.solve();
    recordLinearSolve
    (
      omegaaLinearSolveStats,
      omegaaLinearSolvePerformance.initialResidual(),
      omegaaLinearSolvePerformance.finalResidual(),
      omegaaLinearSolvePerformance.nIterations()
    );

    applyDualConsistentOmegaaBC();
  }
  else
  {
    omegaa = stageCOmegaSourceProxy;
  }

  const volScalarField omegaaMagnitudeCheck(mag(omegaa));
  const scalar omegaaAbsMax = gMax(omegaaMagnitudeCheck);
  if (!(omegaaAbsMax < stageDAdjointMaxMagnitude))
  {
    FatalErrorInFunction
      << "Non-physical omegaa magnitude detected: max(|omegaa|)=" << omegaaAbsMax
      << ". Clipping is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  omegaa.correctBoundaryConditions();
  applyDualConsistentOmegaaBC();
}
