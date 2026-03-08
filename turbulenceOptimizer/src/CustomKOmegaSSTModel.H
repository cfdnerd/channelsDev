{
  const scalar kMinState = gMin(k);
  const scalar omegaMinState = gMin(omega);
  if (!(kMinState > 0.0) || !(omegaMinState > 0.0))
  {
    FatalErrorInFunction
      << "Non-physical turbulence state entering closure model: "
      << "min(k)=" << kMinState << ", min(omega)=" << omegaMinState
      << ". Positivity regularization is disabled by policy; investigate formulation/root cause."
      << exit(FatalError);
  }
  #include "CustomKOmegaSSTWallTreatment.H"

  if (customTurbulenceActive)
  {
    // Use a single primal/adjoint/closure state map without clip-conditioned drift.
    kClosureState = k;
    omegaClosureState = omega;

    if (static_cast<bool>(optimizationCycleUseLogClosureState))
    {
      const dimensionedScalar kLogFloor
      (
        "kLogFloor",
        kClosureState.dimensions(),
        Foam::max(kMin.value(), SMALL)
      );
      const dimensionedScalar omegaLogFloor
      (
        "omegaLogFloor",
        omegaClosureState.dimensions(),
        Foam::max(omegaMin.value(), SMALL)
      );
      const volScalarField kLogArgument
      (
        kClosureState/kLogFloor
      );
      const volScalarField omegaLogArgument
      (
        omegaClosureState/omegaLogFloor
      );
      logKClosureState = Foam::log(kLogArgument);
      logOmegaClosureState = Foam::log(omegaLogArgument);
    }
    else
    {
      logKClosureState = scalar(0.0);
      logOmegaClosureState = scalar(0.0);
    }

    const volScalarField& kStateRef =
      static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState) ? kClosureState : k;
    const volScalarField& omegaStateRef =
      static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState) ? omegaClosureState : omega;

    const volTensorField gradU(fvc::grad(U));
    const volScalarField S2(2.0*magSqr(symm(gradU)));
    const volScalarField sqrtK(Foam::sqrt(kStateRef));
    const volScalarField yEff(Foam::max(yWall, yMin));
    const volScalarField omegaCrossDiffDen(Foam::max(omegaStateRef, omegaMin));
    const volScalarField CDkOmega
    (
      2.0*alphaOmega2Coeff*(fvc::grad(kStateRef) & fvc::grad(omegaStateRef))/omegaCrossDiffDen
    );
    const dimensionedScalar smallCD("smallCD", CDkOmega.dimensions(), 1e-20);
    const dimensionedScalar smallS2("smallS2", S2.dimensions(), SMALL);
    const volScalarField CDkOmegaPlus(Foam::max(CDkOmega, smallCD));

    const volScalarField arg1Term1(sqrtK/(betaStarCoeff*omegaStateRef*yEff));
    const volScalarField arg1Term2(500.0*nu/(sqr(yEff)*omegaStateRef));
    const volScalarField arg1Term3(4.0*alphaOmega2Coeff*kStateRef/(CDkOmegaPlus*sqr(yEff)));
    const volScalarField arg1
    (
      Foam::min
      (
        Foam::min(Foam::max(arg1Term1, arg1Term2), arg1Term3),
        scalar(10.0)
      )
    );
    const volScalarField arg2
    (
      Foam::min
      (
        Foam::max(2.0*sqrtK/(betaStarCoeff*omegaStateRef*yEff), 500.0*nu/(sqr(yEff)*omegaStateRef)),
        scalar(100.0)
      )
    );

    F1SST = Foam::tanh(Foam::pow(arg1, 4.0));
    F2SST = Foam::tanh(sqr(arg2));

    alphaKField = F1SST*(alphaK1Coeff - alphaK2Coeff) + alphaK2Coeff;
    alphaOmegaField = F1SST*(alphaOmega1Coeff - alphaOmega2Coeff) + alphaOmega2Coeff;
    betaField = F1SST*(beta1Coeff - beta2Coeff) + beta2Coeff;
    gammaField = F1SST*(gamma1Coeff - gamma2Coeff) + gamma2Coeff;

    const volScalarField sqrtS2(Foam::sqrt(Foam::max(S2, smallS2)));
    const volScalarField omegaLimiterBranchK(a1Coeff*omegaStateRef);
    const volScalarField omegaLimiterBranchS(b1Coeff*F2SST*sqrtS2);
    const dimensionedScalar omegaLimiterEps
    (
      "omegaLimiterEps",
      omegaLimiterBranchK.dimensions(),
      Foam::max(1e-3*Foam::max(omegaMin.value(), SMALL), SMALL)
    );
    const volScalarField omegaLimiterBlendDen
    (
      Foam::sqrt
      (
        sqr(omegaLimiterBranchK - omegaLimiterBranchS)
        + sqr(omegaLimiterEps)
      )
    );
    const volScalarField omegaLimiter
    (
      0.5*(omegaLimiterBranchK + omegaLimiterBranchS + omegaLimiterBlendDen)
    );
    const volScalarField omegaLimiterWeightK
    (
      0.5*(1.0 + (omegaLimiterBranchK - omegaLimiterBranchS)/omegaLimiterBlendDen)
    );
    const dimensionedScalar smallOmegaLimiter
    (
      "smallOmegaLimiter",
      omegaLimiter.dimensions(),
      SMALL
    );
    const volScalarField omegaLimiterPositive
    (
      Foam::max(omegaLimiter, smallOmegaLimiter)
    );
    const volScalarField omegaNutDenominator(omegaLimiterPositive);
    const volScalarField rawNut(a1Coeff*kStateRef/omegaNutDenominator);
    nut = rawNut;
    const scalar nutMinValue = gMin(nut);
    const scalar nutMaxValue = gMax(nut);
    if (!(nutMinValue > 0.0) || !(nutMaxValue < GREAT))
    {
      FatalErrorInFunction
        << "Non-physical nut predicted by closure model: "
        << "min(nut)=" << nutMinValue << ", max(nut)=" << nutMaxValue
        << ". Clipping is disabled by policy; investigate formulation/root cause."
        << exit(FatalError);
    }

    auto computeNutFromState =
    [&](const volScalarField& kState, const volScalarField& omegaState) -> volScalarField
    {
      const volScalarField sqrtKState(Foam::sqrt(kState));
      const volScalarField arg2State
      (
        Foam::min
        (
          Foam::max
          (
            2.0*sqrtKState/(betaStarCoeff*omegaState*yEff),
            500.0*nu/(sqr(yEff)*omegaState)
          ),
          scalar(100.0)
        )
      );
      const volScalarField F2State(Foam::tanh(sqr(arg2State)));
      const volScalarField omegaLimiterBranchKState(a1Coeff*omegaState);
      const volScalarField omegaLimiterBranchSState(b1Coeff*F2State*sqrtS2);
      const volScalarField omegaLimiterBlendDenState
      (
        Foam::sqrt
        (
          sqr(omegaLimiterBranchKState - omegaLimiterBranchSState)
          + sqr(omegaLimiterEps)
        )
      );
      const volScalarField omegaLimiterState
      (
        0.5
        *(
          omegaLimiterBranchKState
          + omegaLimiterBranchSState
          + omegaLimiterBlendDenState
        )
      );
      const volScalarField omegaLimiterStatePositive
      (
        Foam::max(omegaLimiterState, smallOmegaLimiter)
      );
      const volScalarField omegaNutDenominatorState(omegaLimiterStatePositive);
      const volScalarField rawNutState(a1Coeff*kState/omegaNutDenominatorState);
      return volScalarField
      (
        rawNutState
      );
    };

    if (static_cast<bool>(optClosureDerivativeFDEnable))
    {
      // Use multiplicative (log-space) perturbations to preserve positivity.
      const scalar closureFdLogStep =
        Foam::min(Foam::max(optClosureDerivativeFDRelStep, scalar(1e-7)), scalar(0.25));
      const scalar closureFdExpPlus = Foam::exp(closureFdLogStep);
      const scalar closureFdExpMinus = 1.0/closureFdExpPlus;
      const volScalarField kPlus(kStateRef*closureFdExpPlus);
      const volScalarField kMinus(kStateRef*closureFdExpMinus);
      const volScalarField nutKPlus(computeNutFromState(kPlus, omegaStateRef));
      const volScalarField nutKMinus(computeNutFromState(kMinus, omegaStateRef));
      const dimensionedScalar kSpanFloor
      (
        "kSpanFloor",
        k.dimensions(),
        Foam::max(1e-10*Foam::max(kMin.value(), SMALL), SMALL)
      );
      dNutDkBlock =
        (nutKPlus - nutKMinus)/Foam::max(kPlus - kMinus, kSpanFloor);

      const volScalarField omegaPlus(omegaStateRef*closureFdExpPlus);
      const volScalarField omegaMinus(omegaStateRef*closureFdExpMinus);
      const volScalarField nutOmegaPlus(computeNutFromState(kStateRef, omegaPlus));
      const volScalarField nutOmegaMinus(computeNutFromState(kStateRef, omegaMinus));
      const dimensionedScalar omegaSpanFloor
      (
        "omegaSpanFloor",
        omega.dimensions(),
        Foam::max(1e-10*Foam::max(omegaMin.value(), SMALL), SMALL)
      );
      dNutDOmegaBlock =
        (nutOmegaPlus - nutOmegaMinus)/Foam::max(omegaPlus - omegaMinus, omegaSpanFloor);
    }
    else
    {
      dNutDkBlock = a1Coeff/omegaNutDenominator;
      dNutDOmegaBlock =
        -sqr(a1Coeff)*kStateRef/sqr(omegaNutDenominator)
        *omegaLimiterWeightK;
    }

    dAlphaEffDkBlock = dNutDkBlock/Prt;
    dAlphaEffDOmegaBlock = dNutDOmegaBlock/Prt;

    const volScalarField rawPk(nut*S2);
    PkField = rawPk;

    omegaCrossDiffField =
      2.0
     *(1.0 - F1SST)
     *alphaOmega2Coeff
     *(fvc::grad(kStateRef) & fvc::grad(omegaStateRef))
     /omegaCrossDiffDen;

    alphat = nut/Prt;
    nuEff = nu + nut;
    alphaEff = DT + alphat;
  }
  else
  {
    F1SST = scalar(1.0);
    F2SST = scalar(1.0);
    alphaKField = alphaK2Coeff;
    alphaOmegaField = alphaOmega2Coeff;
    betaField = beta2Coeff;
    gammaField = gamma2Coeff;
    PkField = dimensionedScalar("zeroPkField", PkField.dimensions(), 0.0);
    omegaCrossDiffField =
      dimensionedScalar
      (
        "zeroOmegaCrossDiffField",
        omegaCrossDiffField.dimensions(),
        0.0
      );
    nut = dimensionedScalar("zeroNutField", nut.dimensions(), 0.0);
    alphat = dimensionedScalar("zeroAlphatField", alphat.dimensions(), 0.0);
    dNutDkBlock =
      dimensionedScalar("zeroDNutDk", dNutDkBlock.dimensions(), 0.0);
    dNutDOmegaBlock =
      dimensionedScalar("zeroDNutDOmega", dNutDOmegaBlock.dimensions(), 0.0);
    dAlphaEffDkBlock =
      dimensionedScalar("zeroDAlphaEffDk", dAlphaEffDkBlock.dimensions(), 0.0);
    dAlphaEffDOmegaBlock =
      dimensionedScalar("zeroDAlphaEffDOmega", dAlphaEffDOmegaBlock.dimensions(), 0.0);
    kClosureState = k;
    omegaClosureState = omega;
    logKClosureState = scalar(0.0);
    logOmegaClosureState = scalar(0.0);
    nuEff = nu;
    alphaEff = DT;
  }
}
