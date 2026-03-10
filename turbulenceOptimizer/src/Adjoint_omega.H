if (customTurbulenceActive)
{
  const bool optimizationCycleUseTransportedClosureAdjoints =
    turbulenceAdjointClosureMode
    && closureAdjointTransportedActive
    && static_cast<bool>(optDualConsistentStateMapEnforce);
  const bool optimizationCycleUseResidualDraftClosureAdjoints =
    optimizationCycleUseTransportedClosureAdjoints
    && closureAdjointFormulationActive == word("transportedResidualDraft");
  const bool optimizationCycleUseCrossCoupledClosureAdjoints =
    optimizationCycleUseTransportedClosureAdjoints
    && closureAdjointFormulationActive == word("transportedCrossCoupledV3");
  const bool optimizationCycleUseSAClosureAdjoint =
    turbulenceProModeActive
    && proModeClosureAdjointSolvesActive
    && proModeSAAdjointClosureActive;
  const bool proModeClosureAdjointMode =
      turbulenceProModeActive
   && proModeAdjointStabilityPackActive
   && proModeClosureAdjointSolvesActive;
  const scalar proModeWarmupProgressClosure =
      proModeClosureAdjointMode
    ? Foam::min
      (
        scalar(1.0),
        Foam::max
        (
          scalar(0.0),
          scalar(opt - 1)
         /Foam::max(scalar(proModeAdjointWarmupIterations - 1), scalar(1.0))
        )
      )
    : scalar(1.0);
  const scalar proModeRollbackFreeProgressClosure =
      proModeClosureAdjointMode
    ? Foam::min
      (
        scalar(1.0),
        scalar(proModeAdjointRollbackFreeStreak)
       /Foam::max(scalar(proModeAdjointRampRollbackFreeWindow), scalar(1.0))
      )
    : scalar(1.0);
  const scalar proModeClosureRamp =
      proModeClosureAdjointMode
    ? Foam::min
      (
        scalar(1.0),
        scalar(0.5)*(proModeWarmupProgressClosure + proModeRollbackFreeProgressClosure)
      )
    : scalar(1.0);
  const scalar proModeClosureEqnRelaxFactor =
      proModeClosureAdjointMode
    ? proModeClosureEqnRelaxMin
    + (proModeClosureEqnRelaxMax - proModeClosureEqnRelaxMin)*proModeClosureRamp
    : scalar(1.0);
  const scalar proModeClosurePseudoDecay =
      proModeClosureAdjointMode
    ? Foam::max
      (
        proModePseudoInvTauFloorRatio,
        scalar(1.0)
      - (scalar(1.0) - proModePseudoInvTauFloorRatio)*proModeWarmupProgressClosure
      )
    : scalar(1.0);
  const scalar omegaaPseudoInvTauScale =
      proModeClosureAdjointMode
    ? proModeClosurePseudoInvTauScale*proModeClosurePseudoDecay
    : scalar(0.75);
  const bool useProModeRobustClosureLinearSolvers =
      proModeClosureAdjointMode
   && proModeRobustLinearSolverActive;
  auto applyDualConsistentOmegaaBC =
  [&]()
  {
    if (optimizationCycleUseTransportedClosureAdjoints && closureAdjointDualConsistentBCActive)
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
  auto applyZeroOmegaaBC =
  [&]()
  {
    forAll(omegaa.boundaryFieldRef(), patchi)
    {
      omegaa.boundaryFieldRef()[patchi] == scalar(0.0);
    }
  };

  if (optimizationCycleUseSAClosureAdjoint)
  {
    omegaa = dimensionedScalar("zeroOmegaaAdjoint", omegaa.dimensions(), 0.0);
    applyZeroOmegaaBC();
    omegaa.correctBoundaryConditions();
    applyZeroOmegaaBC();
  }
  else
  {
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
        closureAdjointInvTauActive
      );
      const dimensionedScalar omegaaPseudoInvTau
      (
        "omegaaPseudoInvTau",
        dimless/dimTime,
        omegaaPseudoInvTauScale*closureAdjointInvTauActive
      );
      const volScalarField omegaaOldState(omegaa);
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
        const volScalarField& kaCrossReference =
          (proModeClosureAdjointMode && proModeLaggedClosureInitialized)
        ? proModeLaggedKa
        : ka;
        const volScalarField kaMagnitude(mag(kaCrossReference));
        const dimensionedScalar kaScale
        (
          "kaScale",
          kaCrossReference.dimensions(),
          Foam::max(gMax(kaMagnitude), SMALL)
        );
        optimizationCycleOmegaCrossSource =
          optimizationCycleClosureInvTau
         *closureAdjointCrossCouplingActive
         *(kaCrossReference/kaScale);
      }
      const volScalarField omegaaTransportDiffusivity
      (
        turbulenceAdjointTransportDiffusionScaleActive*(nu + alphaOmegaField*nut)
      );
      fvScalarMatrix omegaaEqn
      (
          fvm::div(phi, omegaa)
        - fvm::laplacian(omegaaTransportDiffusivity, omegaa)
        + fvm::Sp(optimizationCycleClosureInvTau, omegaa)
        + fvm::Sp(optimizationCycleOmegaReactionCoeff, omegaa)
        + fvm::Sp(omegaaPseudoInvTau, omegaa)
       ==
          optimizationCycleClosureInvTau*adjointCycleOmegaSourceProxy
        + optimizationCycleOmegaCrossSource
        + omegaaPseudoInvTau*omegaaOldState
      );

      if (proModeClosureAdjointMode)
      {
        omegaaEqn.relax(proModeClosureEqnRelaxFactor);
      }
      else
      {
        omegaaEqn.relax();
      }

      const auto omegaaLinearSolvePerformance =
        useProModeRobustClosureLinearSolvers
      ? omegaaEqn.solve(proModeRobustNonsymmetricSolverControls)
      : omegaaEqn.solve();
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
  }

  const volScalarField omegaaMagnitudeCheck(mag(omegaa));
  const scalar omegaaAbsMax = gMax(omegaaMagnitudeCheck);
  if (!(omegaaAbsMax < optimizationCycleAdjointMaxMagnitude))
  {
    FatalErrorInFunction
      << "Non-physical omegaa magnitude detected: max(|omegaa|)=" << omegaaAbsMax
      << ". Hard limiting is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  omegaa.correctBoundaryConditions();
  if (optimizationCycleUseSAClosureAdjoint)
  {
    applyZeroOmegaaBC();
  }
  else
  {
    applyDualConsistentOmegaaBC();
  }
  if (proModeClosureAdjointMode)
  {
    proModeLaggedKa = ka;
    proModeLaggedOmegaa = omegaa;
    proModeLaggedKa.correctBoundaryConditions();
    proModeLaggedOmegaa.correctBoundaryConditions();
    proModeLaggedClosureInitialized = true;
  }
}
