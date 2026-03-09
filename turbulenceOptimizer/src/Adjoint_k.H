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
  const scalar kaPseudoInvTauScale =
      proModeClosureAdjointMode
    ? proModeClosurePseudoInvTauScale*proModeClosurePseudoDecay
    : scalar(0.75);
  const bool useProModeRobustClosureLinearSolvers =
      proModeClosureAdjointMode
   && proModeRobustLinearSolverActive;
  auto applyDualConsistentKaBC =
  [&]()
  {
    if (optimizationCycleUseTransportedClosureAdjoints && closureAdjointDualConsistentBCActive)
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

  if (optimizationCycleUseSAClosureAdjoint)
  {
    const scalar saSigma = 2.0/3.0;
    const scalar saCb1 = 0.1355;
    const scalar saCb2 = 0.622;
    const scalar saKappa = 0.41;
    const scalar saCw2 = 0.3;
    const scalar saCw3 = 2.0;
    const scalar saCv1 = 7.1;
    const scalar saCt3 = 1.2;
    const scalar saCt4 = 0.5;
    const scalar saCw1 = saCb1/(saKappa*saKappa) + (1.0 + saCb2)/saSigma;
    const scalar saCv1Cubed = saCv1*saCv1*saCv1;
    const scalar saCw3Pow6 = Foam::pow(saCw3, 6.0);
    const scalar saInvKappa2 = 1.0/(saKappa*saKappa);

    const volTensorField adjointCycleGradU(fvc::grad(U));
    const volTensorField adjointCycleVorticityTensor(skew(adjointCycleGradU));
    const volScalarField adjointCycleVorticityMagnitude
    (
      Foam::sqrt(scalar(2.0))*mag(adjointCycleVorticityTensor)
    );

    volScalarField saNuTildaState
    (
      IOobject
      (
        "saNuTildaState",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("zeroSaNuTildaState", nu.dimensions(), 0.0)
    );
    volScalarField dNutDNuTilda
    (
      IOobject
      (
        "dNutDNuTilda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
      ),
      mesh,
      dimensionedScalar("zeroDNutDNuTilda", dimless, 0.0)
    );
    const scalarField& nutInternal = nut.internalField();
    const scalar nuUniform = Foam::max(nu.value(), SMALL);
    scalarField& saNuTildaInternal = saNuTildaState.primitiveFieldRef();
    scalarField& dNutDNuTildaInternal = dNutDNuTilda.primitiveFieldRef();
    forAll(saNuTildaInternal, celli)
    {
      const scalar nuCell = nuUniform;
      const scalar nutCell = Foam::max(nutInternal[celli], scalar(0.0));
      scalar nuTildaCell = Foam::max(nutCell, scalar(1e-6)*nuCell);
      for (label saNewtonIter = 0; saNewtonIter < 6; ++saNewtonIter)
      {
        const scalar chi = nuTildaCell/nuCell;
        const scalar chi2 = chi*chi;
        const scalar chi3 = chi2*chi;
        const scalar denom = Foam::max(chi3 + saCv1Cubed, SMALL);
        const scalar fv1 = chi3/denom;
        const scalar dfv1dChi =
          3.0*chi2*saCv1Cubed/Foam::max(denom*denom, SMALL);
        const scalar dfv1dNuTilda = dfv1dChi/nuCell;
        const scalar residual = nuTildaCell*fv1 - nutCell;
        const scalar jacobian = Foam::max(fv1 + nuTildaCell*dfv1dNuTilda, SMALL);
        nuTildaCell = Foam::max(nuTildaCell - residual/jacobian, scalar(0.0));
      }
      const scalar chi = nuTildaCell/nuCell;
      const scalar chi2 = chi*chi;
      const scalar chi3 = chi2*chi;
      const scalar denom = Foam::max(chi3 + saCv1Cubed, SMALL);
      const scalar fv1 = chi3/denom;
      const scalar dfv1dChi =
        3.0*chi2*saCv1Cubed/Foam::max(denom*denom, SMALL);
      const scalar dfv1dNuTilda = dfv1dChi/nuCell;
      saNuTildaInternal[celli] = nuTildaCell;
      dNutDNuTildaInternal[celli] =
        Foam::max(fv1 + nuTildaCell*dfv1dNuTilda, scalar(0.0));
    }
    saNuTildaState.correctBoundaryConditions();
    dNutDNuTilda.correctBoundaryConditions();

    const volScalarField dAlphaEffDNuTilda(dNutDNuTilda/Prt);
    const volSymmTensorField adjointCycleSymmGradU(symm(fvc::grad(U)));
    const volSymmTensorField adjointCycleSymmGradUa(symm(fvc::grad(Ua)));
    const volScalarField adjointCycleSAThermalCoupling
    (
      dAlphaEffDNuTilda*(fvc::grad(T) & fvc::grad(optimizationCycleThermalAdjointField))
    );
    const volScalarField adjointCycleSAMomentumCoupling
    (
      2.0*dNutDNuTilda*(adjointCycleSymmGradU && adjointCycleSymmGradUa)
    );
    const volScalarField adjointCycleSAThermalMagnitude(mag(adjointCycleSAThermalCoupling));
    const volScalarField adjointCycleSAMomentumMagnitude(mag(adjointCycleSAMomentumCoupling));
    const scalar adjointCycleSAThermalScaleValue =
      Foam::max(gMax(adjointCycleSAThermalMagnitude), SMALL);
    const scalar adjointCycleSAMomentumScaleValue =
      Foam::max(gMax(adjointCycleSAMomentumMagnitude), SMALL);
    const dimensionedScalar adjointCycleSAThermalScale
    (
      "adjointCycleSAThermalScale",
      adjointCycleSAThermalCoupling.dimensions(),
      adjointCycleSAThermalScaleValue
    );
    const dimensionedScalar adjointCycleSAMomentumScale
    (
      "adjointCycleSAMomentumScale",
      adjointCycleSAMomentumCoupling.dimensions(),
      adjointCycleSAMomentumScaleValue
    );
    volScalarField adjointCycleSASourceProxy
    (
      adjointCycleSAThermalCoupling/adjointCycleSAThermalScale
    + adjointCycleSAMomentumCoupling/adjointCycleSAMomentumScale
    );
    if (static_cast<bool>(optimizationCycleUseLogClosureState))
    {
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
        Foam::max(gMax(optimizationCycleLogClosureMixMagnitude), SMALL)
      );
      adjointCycleSASourceProxy +=
        0.05*optimizationCycleLogClosureMix/optimizationCycleLogClosureScale;
    }

    if (optimizationCycleUseTransportedClosureAdjoints)
    {
      applyDualConsistentKaBC();
      const volScalarField kaOldState(ka);

      const dimensionedScalar optimizationCycleClosureInvTau
      (
        "optimizationCycleClosureInvTau",
        dimless/dimTime,
        closureAdjointInvTauActive
      );
      const dimensionedScalar kaPseudoInvTau
      (
        "kaPseudoInvTau",
        dimless/dimTime,
        kaPseudoInvTauScale*closureAdjointInvTauActive
      );
      volScalarField saTransportDiffusivity
      (
        IOobject
        (
          "saTransportDiffusivity",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroSaTransportDiffusivity", nu.dimensions(), 0.0)
      );
      volScalarField saReactionCoeff
      (
        IOobject
        (
          "saReactionCoeff",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroSaReactionCoeff", dimless/dimTime, 0.0)
      );
      scalarField& saTransportDiffusivityInternal =
        saTransportDiffusivity.primitiveFieldRef();
      scalarField& saReactionCoeffInternal = saReactionCoeff.primitiveFieldRef();
      const scalarField& yWallInternal = yWall.internalField();
      const scalarField& alphaInternal = alpha.internalField();
      const scalarField& vorticityMagInternal =
        adjointCycleVorticityMagnitude.internalField();
      forAll(saReactionCoeffInternal, celli)
      {
        const scalar nuCell = nuUniform;
        const scalar nuTildaCell = Foam::max(saNuTildaInternal[celli], scalar(0.0));
        const scalar yCell = Foam::max(yWallInternal[celli], yMin.value());
        const scalar chi = nuTildaCell/nuCell;
        const scalar chi2 = chi*chi;
        const scalar chi3 = chi2*chi;
        const scalar denom = Foam::max(chi3 + saCv1Cubed, SMALL);
        const scalar fv1 = chi3/denom;
        const scalar fv2 = 1.0 - chi/Foam::max(1.0 + chi*fv1, SMALL);
        const scalar ft2 = saCt3*Foam::exp(-saCt4*chi2);
        const scalar omegaTilde =
          Foam::max
          (
            vorticityMagInternal[celli]
          + nuTildaCell*saInvKappa2*fv2/Foam::max(yCell*yCell, SMALL),
            SMALL
          );
        const scalar r =
          Foam::min
          (
            nuTildaCell
           /Foam::max(omegaTilde*(saKappa*saKappa)*yCell*yCell, SMALL),
            scalar(10.0)
          );
        const scalar r2 = r*r;
        const scalar r6 = r2*r2*r2;
        const scalar g = r + saCw2*(r6 - r);
        const scalar g2 = g*g;
        const scalar g6 = g2*g2*g2;
        const scalar fw =
          g
         *Foam::pow
          (
            (1.0 + saCw3Pow6)/Foam::max(g6 + saCw3Pow6, SMALL),
            1.0/6.0
          );

        saTransportDiffusivityInternal[celli] =
          turbulenceAdjointTransportDiffusionScaleActive
         *Foam::max((nuCell + nuTildaCell)/saSigma, SMALL);
        saReactionCoeffInternal[celli] =
          saCb1*(1.0 - ft2)*omegaTilde
        + saCw1*fw*nuTildaCell/Foam::max(yCell*yCell, SMALL)
        + turbulenceDampingCoeff*alphaInternal[celli];
      }
      saTransportDiffusivity.correctBoundaryConditions();
      saReactionCoeff.correctBoundaryConditions();
      const volScalarField kaPatankarRef
      (
        Foam::max
        (
          mag(kaOldState),
          dimensionedScalar("kaPatankarRefFloor", ka.dimensions(), 1.0)
        )
      );
      const volScalarField kaSourceRaw
      (
        optimizationCycleClosureInvTau*adjointCycleSASourceProxy
      );
      const dimensionedScalar kaSourceZero
      (
        "kaSourceZero",
        kaSourceRaw.dimensions(),
        0.0
      );
      const volScalarField kaSourcePositive
      (
        Foam::max(kaSourceRaw, kaSourceZero)
      );
      const volScalarField kaSourceNegative
      (
        Foam::max(-kaSourceRaw, kaSourceZero)
      );
      const volScalarField kaSourceNegativeCoeff
      (
        kaSourceNegative/kaPatankarRef
      );
      fvScalarMatrix kaEqn
      (
          fvm::div(phi, ka)
        - fvm::laplacian(saTransportDiffusivity, ka)
        + fvm::Sp(optimizationCycleClosureInvTau, ka)
        + fvm::Sp(saReactionCoeff, ka)
        + fvm::Sp(kaPseudoInvTau, ka)
        + fvm::Sp(kaSourceNegativeCoeff, ka)
       ==
          kaSourcePositive
        + kaPseudoInvTau*kaOldState
      );

      if (proModeClosureAdjointMode)
      {
        kaEqn.relax(proModeClosureEqnRelaxFactor);
      }
      else
      {
        kaEqn.relax();
      }

      const auto kaLinearSolvePerformance =
        useProModeRobustClosureLinearSolvers
      ? kaEqn.solve(proModeRobustNonsymmetricSolverControls)
      : kaEqn.solve();
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
      ka = adjointCycleSASourceProxy;
    }
  }
  else
  {
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

    if (optimizationCycleUseTransportedClosureAdjoints)
    {
      applyDualConsistentKaBC();
      const volScalarField kaOldState(ka);

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
      const dimensionedScalar optimizationCycleClosureInvTau
      (
        "optimizationCycleClosureInvTau",
        dimless/dimTime,
        closureAdjointInvTauActive
      );
      const dimensionedScalar kaPseudoInvTau
      (
        "kaPseudoInvTau",
        dimless/dimTime,
        kaPseudoInvTauScale*closureAdjointInvTauActive
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
        const volScalarField& omegaaCrossReference =
          (proModeClosureAdjointMode && proModeLaggedClosureInitialized)
        ? proModeLaggedOmegaa
        : omegaa;
        const volScalarField omegaaMagnitude(mag(omegaaCrossReference));
        const dimensionedScalar omegaaScale
        (
          "omegaaScale",
          omegaaCrossReference.dimensions(),
          Foam::max(gMax(omegaaMagnitude), SMALL)
        );
        optimizationCycleKCrossSource =
          optimizationCycleClosureInvTau
         *closureAdjointCrossCouplingActive
         *(omegaaCrossReference/omegaaScale);
      }
      const volScalarField kaPatankarRef
      (
        Foam::max
        (
          mag(kaOldState),
          dimensionedScalar("kaPatankarRefFloor", ka.dimensions(), 1.0)
        )
      );
      const volScalarField kaSourceRaw
      (
          optimizationCycleClosureInvTau*adjointCycleKSourceProxy
        + optimizationCycleKCrossSource
      );
      const dimensionedScalar kaSourceZero
      (
        "kaSourceZero",
        kaSourceRaw.dimensions(),
        0.0
      );
      const volScalarField kaSourcePositive
      (
        Foam::max(kaSourceRaw, kaSourceZero)
      );
      const volScalarField kaSourceNegative
      (
        Foam::max(-kaSourceRaw, kaSourceZero)
      );
      const volScalarField kaSourceNegativeCoeff
      (
        kaSourceNegative/kaPatankarRef
      );
      const volScalarField kaTransportDiffusivity
      (
        turbulenceAdjointTransportDiffusionScaleActive*(nu + alphaKField*nut)
      );
      fvScalarMatrix kaEqn
      (
          fvm::div(phi, ka)
        - fvm::laplacian(kaTransportDiffusivity, ka)
        + fvm::Sp(optimizationCycleClosureInvTau, ka)
        + fvm::Sp(optimizationCycleKReactionCoeff, ka)
        + fvm::Sp(kaPseudoInvTau, ka)
        + fvm::Sp(kaSourceNegativeCoeff, ka)
       ==
          kaSourcePositive
        + kaPseudoInvTau*kaOldState
      );

      if (proModeClosureAdjointMode)
      {
        kaEqn.relax(proModeClosureEqnRelaxFactor);
      }
      else
      {
        kaEqn.relax();
      }

      const auto kaLinearSolvePerformance =
        useProModeRobustClosureLinearSolvers
      ? kaEqn.solve(proModeRobustNonsymmetricSolverControls)
      : kaEqn.solve();
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
      ka = adjointCycleKSourceProxy;
    }
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
