if (customTurbulenceActive)
{
  const scalar kMinBefore = gMin(k);
  const scalar omegaMinBefore = gMin(omega);
  if (!(kMinBefore > 0.0))
  {
    FatalErrorInFunction
      << "Detected non-physical closure state before omega solve: "
      << "min(k)=" << kMinBefore << ", min(omega)=" << omegaMinBefore
      << ". Hard limiting is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  if (!(omegaMinBefore > 0.0))
  {
    FatalErrorInFunction
      << "Detected non-physical omega state before omega solve: "
      << "min(omega)=" << omegaMinBefore
      << ". Positivity regularization is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  const scalar nutMinBefore = gMin(nut);
  if (!(nutMinBefore > 0.0))
  {
    FatalErrorInFunction
      << "Detected non-physical eddy viscosity before omega solve: "
      << "min(nut)=" << nutMinBefore
      << ". Hard limiting is disabled by policy; investigate formulation/root cause."
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
  const scalar turbulenceSuppressionRampScale =
    turbulenceSuppressionRampScaleActive;
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
      turbulenceSuppressionRampScale
      *(
          porousTurbulenceSuppressionCoeff*turbulenceSuppressionMask*sqrtS2
        + porousTurbulenceAlphaSinkCoeff*turbulenceSuppressionMask*alpha
       );
  }

  const bool positivityPreservingLinearizationEnabled =
    static_cast<bool>(optKOmegaPositivityPreservingLinearizationEnable);
  const volScalarField omegaTransportDiffusivity
  (
    turbulencePrimalTransportDiffusionScaleActive*(nu + alphaOmegaField*nut)
  );

  // Positivity-preserving source treatment:
  // split raw sources into positive explicit part and negative implicit sink.
  const dimensionedScalar nutRegularizationFloor
  (
    "nutRegularizationFloor",
    nut.dimensions(),
    Foam::max
    (
      nutFloor.value(),
      optClosureNutRegularizationScale*Foam::max(nu.value(), SMALL)
    )
  );
  const volScalarField omegaPatankarRef
  (
    static_cast<bool>(optClosureSmoothDenominatorEnable)
    ? Foam::sqrt(sqr(omega) + sqr(omegaMin))
    : Foam::max(omega, omegaMin)
  );
  const volScalarField omegaNutDenominator
  (
    static_cast<bool>(optClosureSmoothDenominatorEnable)
    ? Foam::sqrt(sqr(nut) + sqr(nutRegularizationFloor))
    : Foam::max(nut, nutRegularizationFloor)
  );
  const volScalarField omegaProductionRaw
  (
    gammaField*turbulenceProductionMask*PkField/omegaNutDenominator
  );
  const volScalarField omegaProductionPositive
  (
    Foam::max
    (
      omegaProductionRaw,
      dimensionedScalar("omegaProdZero", omegaProductionRaw.dimensions(), 0.0)
    )
  );
  const volScalarField omegaProductionNegative
  (
    Foam::max
    (
      -omegaProductionRaw,
      dimensionedScalar("omegaProdNegZero", omegaProductionRaw.dimensions(), 0.0)
    )
  );
  const volScalarField omegaProductionNegativeCoeff
  (
    omegaProductionNegative/omegaPatankarRef
  );

  const volScalarField omegaCrossDiffPositive
  (
    Foam::max
    (
      omegaCrossDiffField,
      dimensionedScalar("omegaCrossDiffZero", omegaCrossDiffField.dimensions(), 0.0)
    )
  );
  const volScalarField omegaCrossDiffNegative
  (
    Foam::max
    (
      -omegaCrossDiffField,
      dimensionedScalar("omegaCrossDiffNegZero", omegaCrossDiffField.dimensions(), 0.0)
    )
  );
  const volScalarField omegaCrossDiffNegativeCoeff
  (
    omegaCrossDiffNegative/omegaPatankarRef
  );

  if (positivityPreservingLinearizationEnabled)
  {
    fvScalarMatrix omegaEqn
    (
        fvm::div(phi, omega)
      - fvm::laplacian(omegaTransportDiffusivity, omega)
     ==
        omegaProductionPositive
      + omegaCrossDiffPositive
      - fvm::Sp(omegaProductionNegativeCoeff, omega)
      - fvm::Sp(omegaCrossDiffNegativeCoeff, omega)
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
  }
  else
  {
    fvScalarMatrix omegaEqn
    (
        fvm::div(phi, omega)
      - fvm::laplacian(omegaTransportDiffusivity, omega)
     ==
        omegaProductionRaw
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
  }

  const scalar omegaMinAfter = gMin(omega);
  const scalar omegaMaxAfter = gMax(omega);
  const bool omegaStatePhysical =
       (omegaMinAfter > 0.0)
    && (omegaMaxAfter < GREAT)
    && std::isfinite(static_cast<double>(omegaMinAfter))
    && std::isfinite(static_cast<double>(omegaMaxAfter));

  if (!omegaStatePhysical)
  {
    FatalErrorInFunction
      << "Non-physical omega state detected after solve: "
      << "min(omega)=" << omegaMinAfter << ", max(omega)=" << omegaMaxAfter
      << ". Positivity regularization is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  #include "CustomKOmegaSSTWallTreatment.H"
}
