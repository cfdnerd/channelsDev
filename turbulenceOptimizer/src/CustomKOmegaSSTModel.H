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

    const bool closureLogRegularizationActive =
      static_cast<bool>(optClosureLogRegularizationEnable);
    const bool closureSmoothDenominatorActive =
      static_cast<bool>(optClosureSmoothDenominatorEnable);
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
      if (closureLogRegularizationActive)
      {
        // Log-domain regularization keeps states strictly positive without hard clipping.
        kClosureState = kLogFloor*(Foam::exp(logKClosureState) + scalar(1.0));
        omegaClosureState = omegaLogFloor*(Foam::exp(logOmegaClosureState) + scalar(1.0));
      }
    }
    else
    {
      logKClosureState = scalar(0.0);
      logOmegaClosureState = scalar(0.0);
      if (closureLogRegularizationActive)
      {
        kClosureState += kMin;
        omegaClosureState += omegaMin;
      }
    }

    const volScalarField& kStateRef =
      static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState) ? kClosureState : k;
    const volScalarField& omegaStateRef =
      static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState) ? omegaClosureState : omega;

    const volTensorField gradU(fvc::grad(U));
    const volScalarField S2(2.0*magSqr(symm(gradU)));
    const volScalarField sqrtK(Foam::sqrt(kStateRef));
    const volScalarField yEff(Foam::max(yWall, yMin));
    const volScalarField omegaCrossDiffDen
    (
      closureSmoothDenominatorActive
      ? Foam::sqrt(sqr(omegaStateRef) + sqr(omegaMin))
      : Foam::max(omegaStateRef, omegaMin)
    );
    auditDivisionDenominator
    (
      omegaCrossDiffDen,
      "CustomKOmegaSSTModel.H::omegaCrossDiffDen",
      SMALL
    );
    const volScalarField CDkOmega
    (
      2.0*alphaOmega2Coeff*(fvc::grad(kStateRef) & fvc::grad(omegaStateRef))/omegaCrossDiffDen
    );
    const dimensionedScalar smallCD("smallCD", CDkOmega.dimensions(), 1e-20);
    const dimensionedScalar smallS2("smallS2", S2.dimensions(), SMALL);
    const volScalarField CDkOmegaPlus
    (
      closureSmoothDenominatorActive
      ? 0.5
       *(
          CDkOmega
        + smallCD
        + Foam::sqrt(sqr(CDkOmega - smallCD) + sqr(smallCD))
        )
      : Foam::max(CDkOmega, smallCD)
    );
    auditDivisionDenominator
    (
      CDkOmegaPlus,
      "CustomKOmegaSSTModel.H::CDkOmegaPlus",
      SMALL
    );

    const volScalarField arg1Term1Den(betaStarCoeff*omegaStateRef*yEff);
    const volScalarField arg1Term2Den(sqr(yEff)*omegaStateRef);
    const volScalarField arg1Term3Den(CDkOmegaPlus*sqr(yEff));
    auditDivisionDenominator
    (
      arg1Term1Den,
      "CustomKOmegaSSTModel.H::arg1Term1 denominator",
      SMALL
    );
    auditDivisionDenominator
    (
      arg1Term2Den,
      "CustomKOmegaSSTModel.H::arg1Term2 denominator",
      SMALL
    );
    auditDivisionDenominator
    (
      arg1Term3Den,
      "CustomKOmegaSSTModel.H::arg1Term3 denominator",
      SMALL
    );
    const volScalarField arg1Term1(sqrtK/arg1Term1Den);
    const volScalarField arg1Term2(500.0*nu/arg1Term2Den);
    const volScalarField arg1Term3(4.0*alphaOmega2Coeff*kStateRef/arg1Term3Den);
    const volScalarField arg1
    (
      Foam::min
      (
        Foam::min(Foam::max(arg1Term1, arg1Term2), arg1Term3),
        scalar(10.0)
      )
    );
    const volScalarField arg2Term1(2.0*sqrtK/arg1Term1Den);
    const volScalarField arg2Term2(500.0*nu/arg1Term2Den);
    const volScalarField arg2
    (
      Foam::min
      (
        Foam::max(arg2Term1, arg2Term2),
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
    auditDivisionDenominator
    (
      omegaLimiterBlendDen,
      "CustomKOmegaSSTModel.H::omegaLimiterBlendDen",
      SMALL
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
      closureSmoothDenominatorActive
      ? Foam::sqrt(sqr(omegaLimiter) + sqr(smallOmegaLimiter))
      : Foam::max(omegaLimiter, smallOmegaLimiter)
    );
    const volScalarField omegaNutDenominator(omegaLimiterPositive);
    auditDivisionDenominator
    (
      omegaNutDenominator,
      "CustomKOmegaSSTModel.H::omegaNutDenominator",
      SMALL
    );
    auto computeNutFromState =
    [&](const volScalarField& kState, const volScalarField& omegaState) -> volScalarField
    {
      const volScalarField omegaStatePositive
      (
        closureSmoothDenominatorActive
        ? Foam::sqrt(sqr(omegaState) + sqr(omegaMin))
        : Foam::max(omegaState, omegaMin)
      );
      const volScalarField sqrtKState(Foam::sqrt(kState));
      const volScalarField arg2StateTerm1Den(betaStarCoeff*omegaStatePositive*yEff);
      const volScalarField arg2StateTerm2Den(sqr(yEff)*omegaStatePositive);
      auditDivisionDenominator
      (
        arg2StateTerm1Den,
        "CustomKOmegaSSTModel.H::computeNutFromState arg2 term1 denominator",
        SMALL
      );
      auditDivisionDenominator
      (
        arg2StateTerm2Den,
        "CustomKOmegaSSTModel.H::computeNutFromState arg2 term2 denominator",
        SMALL
      );
      const volScalarField arg2StateTerm1(2.0*sqrtKState/arg2StateTerm1Den);
      const volScalarField arg2StateTerm2(500.0*nu/arg2StateTerm2Den);
      const volScalarField arg2State
      (
        Foam::min
        (
          Foam::max
          (
            arg2StateTerm1,
            arg2StateTerm2
          ),
          scalar(100.0)
        )
      );
      const volScalarField F2State(Foam::tanh(sqr(arg2State)));
      const volScalarField omegaLimiterBranchKState(a1Coeff*omegaStatePositive);
      const volScalarField omegaLimiterBranchSState(b1Coeff*F2State*sqrtS2);
      const volScalarField omegaLimiterBlendDenState
      (
        Foam::sqrt
        (
          sqr(omegaLimiterBranchKState - omegaLimiterBranchSState)
          + sqr(omegaLimiterEps)
        )
      );
      auditDivisionDenominator
      (
        omegaLimiterBlendDenState,
        "CustomKOmegaSSTModel.H::computeNutFromState omegaLimiterBlendDenState",
        SMALL
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
        closureSmoothDenominatorActive
        ? Foam::sqrt(sqr(omegaLimiterState) + sqr(smallOmegaLimiter))
        : Foam::max(omegaLimiterState, smallOmegaLimiter)
      );
      const volScalarField omegaNutDenominatorState(omegaLimiterStatePositive);
      auditDivisionDenominator
      (
        omegaNutDenominatorState,
        "CustomKOmegaSSTModel.H::computeNutFromState omegaNutDenominatorState",
        SMALL
      );
      const volScalarField rawNutState(a1Coeff*kState/omegaNutDenominatorState);
      return volScalarField
      (
        rawNutState
      );
    };
    const volScalarField rawNut(computeNutFromState(kStateRef, omegaStateRef));
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
      const volScalarField kFdDenominator(Foam::max(kPlus - kMinus, kSpanFloor));
      auditDivisionDenominator
      (
        kFdDenominator,
        "CustomKOmegaSSTModel.H::dNutDk FD denominator",
        SMALL
      );
      dNutDkBlock =
        (nutKPlus - nutKMinus)/kFdDenominator;

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
      const volScalarField omegaFdDenominator
      (
        Foam::max(omegaPlus - omegaMinus, omegaSpanFloor)
      );
      auditDivisionDenominator
      (
        omegaFdDenominator,
        "CustomKOmegaSSTModel.H::dNutDOmega FD denominator",
        SMALL
      );
      dNutDOmegaBlock =
        (nutOmegaPlus - nutOmegaMinus)/omegaFdDenominator;
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
