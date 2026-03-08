if (customTurbulenceActive)
{
  const bool turbulenceStateFailFast = optDiagnosticMode;
  const bool turbulenceStateRegularizationEnabled =
    static_cast<bool>(optTurbulenceStateRegularizationEnable);
  const label kRecoveryAttempts =
    (turbulenceStateFailFast || !turbulenceStateRegularizationEnabled) ? 1 : 3;
  const dimensionedScalar kRecoveryFloor
  (
    "kRecoveryFloor",
    k.dimensions(),
    Foam::max(kMin.value(), SMALL)
  );
  auto regularizeKState =
  [&]()
  {
    const volScalarField kCentered(k - kRecoveryFloor);
    k =
      kRecoveryFloor
      + 0.5
       *(
          kCentered
        + Foam::sqrt(sqr(kCentered) + sqr(kRecoveryFloor))
        );
    k.correctBoundaryConditions();
  };

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
    if (turbulenceStateFailFast || !turbulenceStateRegularizationEnabled)
    {
      FatalErrorInFunction
        << "Detected non-physical k state before k solve: "
        << "min(k)=" << kMinBefore
        << ". Smooth regularization is disabled; investigate formulation/root cause."
        << exit(FatalError);
    }
    if (Pstream::master())
    {
      Info<< "k positivity recovery: pre-solve state min(k)="
          << kMinBefore << ". Applying smooth positivity regularization.\n";
    }
    regularizeKState();
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

  scalar kMinAfter = gMin(k);
  scalar kMaxAfter = gMax(k);
  bool kStatePhysical = false;
  for (label kAttempt = 0; kAttempt < kRecoveryAttempts; ++kAttempt)
  {
    const dimensionedScalar turbulencePseudoTime
    (
      "turbulencePseudoTimeK",
      dimensionSet(0,0,-1,0,0,0,0),
      primalPseudoTimeCoeffActive*(1.0 + 4.0*kAttempt)
    );

    fvScalarMatrix kEqn
    (
        fvm::div(phi, k)
      - fvm::laplacian(nu + alphaKField*nut, k)
      + fvm::Sp(turbulencePseudoTime, k)
     ==
        turbulenceProductionMask*PkField
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

    const auto kLinearSolvePerformance = kEqn.solve();
    recordLinearSolve
    (
        kLinearSolveStats,
        kLinearSolvePerformance.initialResidual(),
        kLinearSolvePerformance.finalResidual(),
        kLinearSolvePerformance.nIterations()
    );

    kMinAfter = gMin(k);
    kMaxAfter = gMax(k);
    kStatePhysical =
         (kMinAfter > 0.0)
      && (kMaxAfter < GREAT)
      && std::isfinite(static_cast<double>(kMinAfter))
      && std::isfinite(static_cast<double>(kMaxAfter));
    if (kStatePhysical)
    {
      break;
    }
    if (kAttempt + 1 < kRecoveryAttempts)
    {
      if (Pstream::master())
      {
        Info<< "k positivity recovery: solve attempt " << (kAttempt + 1)
            << " produced min(k)=" << kMinAfter
            << ", max(k)=" << kMaxAfter
            << ". Retrying after smooth positivity regularization.\n";
      }
      regularizeKState();
    }
  }

  if (!kStatePhysical)
  {
    if (turbulenceStateFailFast || !turbulenceStateRegularizationEnabled)
    {
      FatalErrorInFunction
        << "Non-physical k state detected after solve: "
        << "min(k)=" << kMinAfter << ", max(k)=" << kMaxAfter
        << ". Smooth regularization is disabled; investigate formulation/root cause."
        << exit(FatalError);
    }
    if (Pstream::master())
    {
      Info<< "k positivity recovery exhausted. "
          << "Applying smooth positivity regularization and continuing in production mode.\n";
    }
    regularizeKState();
  }
  #include "CustomKOmegaSSTWallTreatment.H"
}
