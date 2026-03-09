if (customTurbulenceActive)
{
  const scalar kMinBefore = gMin(k);
  const scalar omegaMinBefore = gMin(omega);
  if (!(omegaMinBefore > 0.0))
  {
    FatalErrorInFunction
      << "Detected non-physical closure state before k solve: "
      << "min(k)=" << kMinBefore << ", min(omega)=" << omegaMinBefore
      << ". Clipping is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  if (!(kMinBefore > 0.0))
  {
    FatalErrorInFunction
      << "Detected non-physical k state before k solve: "
      << "min(k)=" << kMinBefore
      << ". Positivity regularization is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }

  volScalarField omegaDestructionState(omega);
  if (static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState))
  {
    omegaDestructionState = omegaClosureState;
  }

  volScalarField turbulenceProductionMask
  (
    Foam::min
    (
      Foam::max
      (
        Foam::pow(xh, porousTurbulenceProductionExponent),
        scalar(porousTurbulenceProductionFloor)
      ),
      scalar(1.0)
    )
  );
  volScalarField turbulenceSuppressionRate
  (
    alpha*scalar(0.0)
  );
  if (porousAwareTurbulenceTransport)
  {
    const volTensorField gradU(fvc::grad(U));
    const volScalarField sqrtS2(Foam::sqrt(2.0*magSqr(symm(gradU))));
    const volScalarField turbulenceSuppressionMask
    (
      Foam::min
      (
        Foam::max
        (
          Foam::pow
          (
            Foam::max(scalar(1.0) - xh, scalar(0.0)),
            porousTurbulenceSuppressionExponent
          ),
          scalar(0.0)
        ),
        scalar(1.0)
      )
    );
    turbulenceSuppressionRate =
        porousTurbulenceSuppressionCoeff*turbulenceSuppressionMask*sqrtS2
      + porousTurbulenceAlphaSinkCoeff*turbulenceSuppressionMask*alpha;
  }

  const bool positivityPreservingLinearizationEnabled =
    static_cast<bool>(optKOmegaPositivityPreservingLinearizationEnable);
  const volScalarField kTransportDiffusivity
  (
    turbulencePrimalTransportDiffusionScaleActive*(nu + alphaKField*nut)
  );
  // Patankar-style source split: keep explicit source non-negative and treat
  // any negative part implicitly as a sink coefficient.
  const volScalarField kPatankarRef
  (
    static_cast<bool>(optClosureSmoothDenominatorEnable)
    ? Foam::sqrt(sqr(k) + sqr(kMin))
    : Foam::max(k, kMin)
  );
  const volScalarField kProductionRaw(turbulenceProductionMask*PkField);
  const volScalarField kProductionPositive
  (
    Foam::max
    (
      kProductionRaw,
      dimensionedScalar("kProductionZero", kProductionRaw.dimensions(), 0.0)
    )
  );
  const volScalarField kProductionNegative
  (
    Foam::max
    (
      -kProductionRaw,
      dimensionedScalar("kProductionNegZero", kProductionRaw.dimensions(), 0.0)
    )
  );
  const volScalarField kProductionNegativeCoeff(kProductionNegative/kPatankarRef);

  auto solveKWithLinearization =
  [&]()
  {
    if (positivityPreservingLinearizationEnabled)
    {
      fvScalarMatrix kEqn
      (
          fvm::div(phi, k)
        - fvm::laplacian(kTransportDiffusivity, k)
       ==
          kProductionPositive
        - fvm::Sp(kProductionNegativeCoeff, k)
        - fvm::Sp
          (
            betaStarCoeff
           *omegaDestructionState,
            k
          )
        - fvm::Sp(turbulenceDampingCoeff*alpha, k)
        - fvm::Sp(turbulenceSuppressionRate, k)
      );
      kEqn.relax();
      return kEqn.solve();
    }
    else
    {
      fvScalarMatrix kEqn
      (
          fvm::div(phi, k)
        - fvm::laplacian(kTransportDiffusivity, k)
       ==
          kProductionRaw
        - fvm::Sp
          (
            betaStarCoeff
           *omegaDestructionState,
            k
          )
        - fvm::Sp(turbulenceDampingCoeff*alpha, k)
        - fvm::Sp(turbulenceSuppressionRate, k)
      );
      kEqn.relax();
      return kEqn.solve();
    }
  };

  const auto kLinearSolvePerformance = solveKWithLinearization();
  recordLinearSolve
  (
      kLinearSolveStats,
      kLinearSolvePerformance.initialResidual(),
      kLinearSolvePerformance.finalResidual(),
      kLinearSolvePerformance.nIterations()
  );

  const scalar kMinAfter = gMin(k);
  const scalar kMaxAfter = gMax(k);
  const bool kStatePhysical =
       (kMinAfter > 0.0)
    && (kMaxAfter < GREAT)
    && std::isfinite(static_cast<double>(kMinAfter))
    && std::isfinite(static_cast<double>(kMaxAfter));

  if (!kStatePhysical)
  {
    FatalErrorInFunction
      << "Non-physical k state detected after solve: "
      << "min(k)=" << kMinAfter << ", max(k)=" << kMaxAfter
      << ". Positivity regularization is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  #include "CustomKOmegaSSTWallTreatment.H"
}
