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
  auto applyDualConsistentOmegaaBC =
  [&]()
  {
    if (optimizationCycleUseTransportedClosureAdjoints && static_cast<bool>(optimizationCycleClosureAdjointDualConsistentBC))
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

  const volSymmTensorField adjointCycleSymmGradU(symm(fvc::grad(U)));
  const volSymmTensorField adjointCycleSymmGradUa(symm(fvc::grad(Ua)));
  const volScalarField& optimizationCycleOmegaStateRef =
    static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
    ? omegaClosureState
    : omega;
  const volScalarField adjointCycleOmegaThermalCoupling
  (
    dAlphaEffDOmegaBlock*(fvc::grad(T) & fvc::grad(optimizationCycleThermalAdjointField))
  );
  const volScalarField adjointCycleOmegaMomentumCoupling
  (
    2.0*dNutDOmegaBlock*(adjointCycleSymmGradU && adjointCycleSymmGradUa)
  );
  const volScalarField adjointCycleOmegaThermalMagnitude(mag(adjointCycleOmegaThermalCoupling));
  const volScalarField adjointCycleOmegaMomentumMagnitude(mag(adjointCycleOmegaMomentumCoupling));
  const scalar adjointCycleOmegaThermalScaleValue =
    Foam::max(gMax(adjointCycleOmegaThermalMagnitude), SMALL);
  const scalar adjointCycleOmegaMomentumScaleValue =
    Foam::max(gMax(adjointCycleOmegaMomentumMagnitude), SMALL);
  const dimensionedScalar adjointCycleOmegaThermalScale
  (
    "adjointCycleOmegaThermalScale",
    adjointCycleOmegaThermalCoupling.dimensions(),
    adjointCycleOmegaThermalScaleValue
  );
  const dimensionedScalar adjointCycleOmegaMomentumScale
  (
    "adjointCycleOmegaMomentumScale",
    adjointCycleOmegaMomentumCoupling.dimensions(),
    adjointCycleOmegaMomentumScaleValue
  );
  volScalarField adjointCycleOmegaSourceProxy
  (
    adjointCycleOmegaThermalCoupling/adjointCycleOmegaThermalScale
  + adjointCycleOmegaMomentumCoupling/adjointCycleOmegaMomentumScale
  );
  if (static_cast<bool>(optimizationCycleUseLogClosureState))
  {
    const volScalarField logOmegaMagnitude(mag(logOmegaClosureState));
    const dimensionedScalar logOmegaScale
    (
      "logOmegaScale",
      logOmegaClosureState.dimensions(),
      Foam::max(gMax(logOmegaMagnitude), SMALL)
    );
    adjointCycleOmegaSourceProxy += 0.05*logOmegaClosureState/logOmegaScale;
  }
  if (optimizationCycleUseTransportedClosureAdjoints)
  {
    applyDualConsistentOmegaaBC();

    const dimensionedScalar optimizationCycleClosureInvTau
    (
      "optimizationCycleClosureInvTau",
      dimless/dimTime,
      optimizationCycleClosureAdjointInvTau
    );
    volScalarField optimizationCycleOmegaReactionCoeff
    (
      scalar(2.0)
     *betaField
     *optimizationCycleOmegaStateRef
    + turbulenceDampingCoeff*alpha
    );
    if (!optimizationCycleUseResidualDraftClosureAdjoints)
    {
      optimizationCycleOmegaReactionCoeff =
        dimensionedScalar
        (
          "zeroOptimizationCycleOmegaReactionCoeff",
          optimizationCycleOmegaReactionCoeff.dimensions(),
          0.0
        );
    }
    volScalarField optimizationCycleOmegaCrossSource
    (
      IOobject
      (
        "optimizationCycleOmegaCrossSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("zeroOptimizationCycleOmegaCrossSource", dimless/dimTime, 0.0)
    );
    if (optimizationCycleUseCrossCoupledClosureAdjoints)
    {
      const volScalarField kaMagnitude(mag(ka));
      const dimensionedScalar kaScale
      (
        "kaScale",
        ka.dimensions(),
        Foam::max(gMax(kaMagnitude), SMALL)
      );
      optimizationCycleOmegaCrossSource =
        optimizationCycleClosureInvTau
       *optimizationCycleClosureCrossCoupling
       *(ka/kaScale);
    }
    fvScalarMatrix omegaaEqn
    (
        fvm::div(phi, omegaa)
      - fvm::laplacian(nu + alphaOmegaField*nut, omegaa)
      + fvm::Sp(optimizationCycleClosureInvTau, omegaa)
      + fvm::Sp(optimizationCycleOmegaReactionCoeff, omegaa)
     ==
        optimizationCycleClosureInvTau*adjointCycleOmegaSourceProxy
      + optimizationCycleOmegaCrossSource
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
    omegaa = adjointCycleOmegaSourceProxy;
  }

  const volScalarField omegaaMagnitudeCheck(mag(omegaa));
  const scalar omegaaAbsMax = gMax(omegaaMagnitudeCheck);
  if (!(omegaaAbsMax < optimizationCycleAdjointMaxMagnitude))
  {
    FatalErrorInFunction
      << "Non-physical omegaa magnitude detected: max(|omegaa|)=" << omegaaAbsMax
      << ". Clipping is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  omegaa.correctBoundaryConditions();
  applyDualConsistentOmegaaBC();
}
