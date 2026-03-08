if (customTurbulenceActive)
{
  const bool turbulenceStateFailFast = optDiagnosticMode;
  const bool turbulenceStateRegularizationEnabled =
    static_cast<bool>(optTurbulenceStateRegularizationEnable);
  const label omegaRecoveryAttempts =
    (turbulenceStateFailFast || !turbulenceStateRegularizationEnabled) ? 1 : 3;
  const dimensionedScalar omegaRecoveryFloor
  (
    "omegaRecoveryFloor",
    omega.dimensions(),
    Foam::max(omegaMin.value(), SMALL)
  );
  auto regularizeOmegaState =
  [&]()
  {
    const volScalarField omegaCentered(omega - omegaRecoveryFloor);
    omega =
      omegaRecoveryFloor
      + 0.5
       *(
          omegaCentered
        + Foam::sqrt(sqr(omegaCentered) + sqr(omegaRecoveryFloor))
        );
    omega.correctBoundaryConditions();
  };

  const scalar kMinBefore = gMin(k);
  const scalar omegaMinBefore = gMin(omega);
  if (!(kMinBefore > 0.0))
  {
    FatalErrorInFunction
      << "Detected non-physical closure state before omega solve: "
      << "min(k)=" << kMinBefore << ", min(omega)=" << omegaMinBefore
      << ". Clipping is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  if (!(omegaMinBefore > 0.0))
  {
    if (turbulenceStateFailFast || !turbulenceStateRegularizationEnabled)
    {
      FatalErrorInFunction
        << "Detected non-physical omega state before omega solve: "
        << "min(omega)=" << omegaMinBefore
        << ". Smooth regularization is disabled; investigate formulation/root cause."
        << exit(FatalError);
    }
    if (Pstream::master())
    {
      Info<< "Omega positivity recovery: pre-solve state min(omega)="
          << omegaMinBefore << ". Applying smooth positivity regularization.\n";
    }
    regularizeOmegaState();
  }
  const scalar nutMinBefore = gMin(nut);
  if (!(nutMinBefore > 0.0))
  {
    FatalErrorInFunction
      << "Detected non-physical eddy viscosity before omega solve: "
      << "min(nut)=" << nutMinBefore
      << ". Clipping is disabled by policy; investigate formulation/root cause."
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

  scalar omegaMinAfter = gMin(omega);
  scalar omegaMaxAfter = gMax(omega);
  bool omegaStatePhysical = false;
  for (label omegaAttempt = 0; omegaAttempt < omegaRecoveryAttempts; ++omegaAttempt)
  {
    const dimensionedScalar turbulencePseudoTime
    (
      "turbulencePseudoTimeOmega",
      dimensionSet(0,0,-1,0,0,0,0),
      primalPseudoTimeCoeffActive*(1.0 + 4.0*omegaAttempt)
    );

    fvScalarMatrix omegaEqn
    (
        fvm::div(phi, omega)
      - fvm::laplacian(nu + alphaOmegaField*nut, omega)
      + fvm::Sp(turbulencePseudoTime, omega)
     ==
        gammaField*turbulenceProductionMask*PkField/nut
      + omegaCrossDiffField
      - fvm::Sp(betaField*omegaDestructionState, omega)
      - fvm::Sp(turbulenceDampingCoeff*alpha, omega)
      - fvm::Sp(turbulenceSuppressionRate, omega)
    );

    omegaEqn.relax();

    const auto omegaLinearSolvePerformance = omegaEqn.solve();
    recordLinearSolve
    (
        omegaLinearSolveStats,
        omegaLinearSolvePerformance.initialResidual(),
        omegaLinearSolvePerformance.finalResidual(),
        omegaLinearSolvePerformance.nIterations()
    );

    omegaMinAfter = gMin(omega);
    omegaMaxAfter = gMax(omega);
    omegaStatePhysical =
         (omegaMinAfter > 0.0)
      && (omegaMaxAfter < GREAT)
      && std::isfinite(static_cast<double>(omegaMinAfter))
      && std::isfinite(static_cast<double>(omegaMaxAfter));
    if (omegaStatePhysical)
    {
      break;
    }
    if (omegaAttempt + 1 < omegaRecoveryAttempts)
    {
      if (Pstream::master())
      {
        Info<< "Omega positivity recovery: solve attempt " << (omegaAttempt + 1)
            << " produced min(omega)=" << omegaMinAfter
            << ", max(omega)=" << omegaMaxAfter
            << ". Retrying after smooth positivity regularization.\n";
      }
      regularizeOmegaState();
    }
  }

  if (!omegaStatePhysical)
  {
    if (turbulenceStateFailFast || !turbulenceStateRegularizationEnabled)
    {
      FatalErrorInFunction
        << "Non-physical omega state detected after solve: "
        << "min(omega)=" << omegaMinAfter << ", max(omega)=" << omegaMaxAfter
        << ". Smooth regularization is disabled; investigate formulation/root cause."
        << exit(FatalError);
    }
    if (Pstream::master())
    {
      Info<< "Omega positivity recovery exhausted. "
          << "Applying smooth positivity regularization and continuing in production mode.\n";
    }
    regularizeOmegaState();
  }
  #include "CustomKOmegaSSTWallTreatment.H"
}
