//***********************************//update mu and lambda
opt++;
const label continuationScheduleIterationForRamps = Foam::max(opt, label(1));
auto linearRampValue =
  [&](const scalar startValue, const scalar endValue, const label startIter, const label endIter) -> scalar
  {
    if (continuationScheduleIterationForRamps <= startIter)
    {
      return startValue;
    }
    if (continuationScheduleIterationForRamps >= endIter)
    {
      return endValue;
    }
    const scalar rampSpan =
      Foam::max
      (
        scalar(endIter - startIter),
        scalar(1.0)
      );
    const scalar rampProgress =
      scalar(continuationScheduleIterationForRamps - startIter)/rampSpan;
    return startValue + (endValue - startValue)*rampProgress;
  };
projectionBetaScaleActive =
  linearRampValue
  (
    optProjectionBetaScaleStart,
    optProjectionBetaScaleEnd,
    optProjectionBetaRampStartIteration,
    optProjectionBetaRampEndIteration
  );
projectionDelStepActive =
  linearRampValue
  (
    optProjectionDelStepStart,
    optProjectionDelStepEnd,
    optProjectionDelStepRampStartIteration,
    optProjectionDelStepRampEndIteration
  );
filterRadiusScaleActive =
  static_cast<bool>(optFilterStrengthRampEnable)
  ? linearRampValue
    (
      scalar(1.0),
      optFilterRadiusScaleMin,
      optFilterStrengthRampStartIteration,
      optFilterStrengthRampEndIteration
    )
  : scalar(1.0);
materialInterpolationExponentActive =
  static_cast<bool>(optMaterialInterpolationExponentEnable)
  ? linearRampValue
    (
      optMaterialInterpolationExponentStart,
      optMaterialInterpolationExponentEnd,
      optMaterialInterpolationExponentRampStartIteration,
      optMaterialInterpolationExponentRampEndIteration
    )
  : scalar(1.0);
grayPenaltyWeightActive =
  static_cast<bool>(optGrayPenaltyEnable)
  ? linearRampValue
    (
      optGrayPenaltyWeightStart,
      optGrayPenaltyWeightEnd,
      optGrayPenaltyRampStartIteration,
      optGrayPenaltyRampEndIteration
    )
  : scalar(0.0);
b.value() =
  filterBCoeffBaseValue/sqr(Foam::max(filterRadiusScaleActive, scalar(1e-6)));
const scalar governorSlowFactorClamped =
  Foam::max
  (
    Foam::min(optStabilityContinuationGovernorSlowFactor, scalar(1.0)),
    scalar(0.05)
  );
const scalar governorMoveLimitScaleClamped =
  Foam::max
  (
    Foam::min(optStabilityContinuationGovernorMoveLimitScale, scalar(1.0)),
    scalar(0.05)
  );
const bool governorPenaltyStress =
  hasPreviousGovernorMetrics
  &&
  (
       previousConstraintPenaltyForGovernor >= optStabilityContinuationGovernorConstraintPenaltyThreshold
    || previousLinearPenaltyForGovernor >= optStabilityContinuationGovernorLinearPenaltyThreshold
  );
const bool governorTrustUndefinedStress =
  hasPreviousGovernorMetrics
  && static_cast<bool>(optStabilityContinuationGovernorTrustUndefinedBlocks)
  && !previousTrustRatioDefinedForGovernor;
continuationGovernorStressActive =
    static_cast<bool>(optStabilityContinuationGovernorEnable)
 && (governorPenaltyStress || governorTrustUndefinedStress);
const scalar continuationRampScale =
  continuationGovernorStressActive
  ? governorSlowFactorClamped
  : scalar(1.0);
const bool mmaStagnationRecoveryEligible =
    static_cast<bool>(optMMAStagnationRecoveryEnable)
 && optimizerTinyStepStreak >= optMMAStagnationWindow;
optimizerMMAStagnationRecoveryActive = mmaStagnationRecoveryEligible;
if (mmaStagnationRecoveryEligible && static_cast<bool>(optMMAStagnationAsymptoteResetEnable))
{
  optimizerResetMMAAsymptotesNextIter = true;
}
if (turbulenceProModeActive)
{
  continuationIteration = Foam::max(opt, label(1));
  continuationScheduleIterationActive = continuationIteration;
  optimizerMoveLimitActive =
    continuationGovernorStressActive
    ? Foam::max(movlim*governorMoveLimitScaleClamped, scalar(1e-4))
    : movlim;
  if (mmaStagnationRecoveryEligible)
  {
    const scalar stagnationMoveLimitCap =
      Foam::max(Foam::min(optMMAStagnationMoveLimitCap, scalar(1.0)), scalar(1e-4));
    optimizerMoveLimitActive =
      Foam::min
      (
        Foam::max(optimizerMoveLimitActive*optMMAStagnationMoveLimitBoost, scalar(1e-4)),
        stagnationMoveLimitCap
      );
  }
  if (static_cast<bool>(optTurbulenceSuppressionRampEnable))
  {
    if (continuationScheduleIterationActive <= optTurbulenceSuppressionRampStartIteration)
    {
      turbulenceSuppressionRampScaleActive = optTurbulenceSuppressionRampMinScale;
    }
    else if (continuationScheduleIterationActive >= optTurbulenceSuppressionRampEndIteration)
    {
      turbulenceSuppressionRampScaleActive = scalar(1.0);
    }
    else
    {
      const scalar suppressionRampSpan =
        Foam::max
        (
          scalar(optTurbulenceSuppressionRampEndIteration - optTurbulenceSuppressionRampStartIteration),
          scalar(1.0)
        );
      const scalar suppressionRampProgress =
        (continuationScheduleIterationActive - optTurbulenceSuppressionRampStartIteration)
       /suppressionRampSpan;
      turbulenceSuppressionRampScaleActive =
        optTurbulenceSuppressionRampMinScale
      + (scalar(1.0) - optTurbulenceSuppressionRampMinScale)*suppressionRampProgress;
    }
  }
  else
  {
    turbulenceSuppressionRampScaleActive = scalar(1.0);
  }
  if(solid_area)
  {
     forAll(cells_solid, i)
     {
        xh[cells_solid[i]]=0;
     }
  }
  if(fluid_area)
  {
     forAll(cells_fluid, i)
     {
        xh[cells_fluid[i]]=1.0;
     }
  }
  xh.correctBoundaryConditions();

  if(opt<=63)
  {
    alphaMax=alphamax*(opt*continuationRampScale/7.0+1.0);
  }
  if(opt>100)
  {
    alphaMax=alphaMax*(1.0 + 0.05*continuationRampScale);
  }
  alphaMax.value()=Foam::min(alphaMax.value(),alphaMAX.value());
  if(opt>=80)
  {
    qu=0.005+(opt-80)*1e-4*continuationRampScale;
  }
  qu=Foam::min(qu,0.01);
  const volScalarField xhMaterial
  (
    Foam::pow
    (
      Foam::max(xh, scalar(0.0)),
      materialInterpolationExponentActive
    )
  );
  const volScalarField xhBlendDenominator(qu + xhMaterial);
  auditDivisionDenominator
  (
    xhBlendDenominator,
    "update.H::alpha/DT denominator (qu+xhMaterial)",
    SMALL
  );
  alpha=alphaMax*qu*(1-xhMaterial)/xhBlendDenominator;
  DT=(ks+(kf-ks)*xhMaterial*(1+qu)/xhBlendDenominator)/rhoCp;
}
else
{
  continuationIteration = Foam::max(opt, label(1));
  continuationScheduleIterationActive = continuationIteration;
  if (static_cast<bool>(optTurbulenceSuppressionRampEnable))
  {
    if (continuationScheduleIterationActive <= optTurbulenceSuppressionRampStartIteration)
    {
      turbulenceSuppressionRampScaleActive = optTurbulenceSuppressionRampMinScale;
    }
    else if (continuationScheduleIterationActive >= optTurbulenceSuppressionRampEndIteration)
    {
      turbulenceSuppressionRampScaleActive = scalar(1.0);
    }
    else
    {
      const scalar suppressionRampSpan =
        Foam::max
        (
          scalar(optTurbulenceSuppressionRampEndIteration - optTurbulenceSuppressionRampStartIteration),
          scalar(1.0)
        );
      const scalar suppressionRampProgress =
        (continuationScheduleIterationActive - optTurbulenceSuppressionRampStartIteration)
       /suppressionRampSpan;
      turbulenceSuppressionRampScaleActive =
        optTurbulenceSuppressionRampMinScale
      + (scalar(1.0) - optTurbulenceSuppressionRampMinScale)*suppressionRampProgress;
    }
  }
  else
  {
    turbulenceSuppressionRampScaleActive = scalar(1.0);
  }
  const scalar adaptiveMoveLimitMinClamped =
    Foam::max(Foam::min(optAdaptiveMoveLimitMin, scalar(1.0)), scalar(1e-4));
  const scalar adaptiveMoveLimitMaxClamped =
    Foam::max
    (
      Foam::min(optAdaptiveMoveLimitMax, scalar(1.0)),
      adaptiveMoveLimitMinClamped
    );
  const scalar adaptiveTrustRatioShrinkClamped =
    Foam::max(optAdaptiveMoveLimitTrustRatioShrink, scalar(0.0));
  const scalar adaptiveTrustRatioGrowClamped =
    Foam::max(optAdaptiveMoveLimitTrustRatioGrow, adaptiveTrustRatioShrinkClamped);
  const scalar adaptiveShrinkFactorClamped =
    Foam::max(Foam::min(optAdaptiveMoveLimitShrinkFactor, scalar(1.0)), scalar(0.05));
  const scalar adaptiveGrowFactorClamped =
    Foam::max(optAdaptiveMoveLimitGrowFactor, scalar(1.0));
  const scalar adaptiveNegativeTrustDampingClamped =
    Foam::max(Foam::min(optAdaptiveMoveLimitNegativeTrustDamping, scalar(1.0)), scalar(0.05));
  if (static_cast<bool>(optAdaptiveMoveLimitEnable))
  {
    scalar adaptiveMoveLimitCandidate =
      Foam::max(optimizerMoveLimitActive, adaptiveMoveLimitMinClamped);
    if (!hasPreviousTrustRatio)
    {
      adaptiveMoveLimitCandidate = movlim;
    }
    if (hasPreviousTrustRatio)
    {
      if (previousTrustRatio < adaptiveTrustRatioShrinkClamped)
      {
        adaptiveMoveLimitCandidate *= adaptiveShrinkFactorClamped;
      }
      else if (previousTrustRatio > adaptiveTrustRatioGrowClamped)
      {
        adaptiveMoveLimitCandidate *= adaptiveGrowFactorClamped;
      }
      if (previousTrustRatio < scalar(0.0))
      {
        adaptiveMoveLimitCandidate *= adaptiveNegativeTrustDampingClamped;
      }
    }
    optimizerMoveLimitActive =
      Foam::max
      (
        Foam::min(adaptiveMoveLimitCandidate, adaptiveMoveLimitMaxClamped),
        adaptiveMoveLimitMinClamped
      );
  }
  else
  {
    optimizerMoveLimitActive = movlim;
  }
  if (continuationGovernorStressActive)
  {
    optimizerMoveLimitActive =
      Foam::max(optimizerMoveLimitActive*governorMoveLimitScaleClamped, scalar(1e-4));
  }
  if (mmaStagnationRecoveryEligible)
  {
    const scalar stagnationMoveLimitCap =
      Foam::max(Foam::min(optMMAStagnationMoveLimitCap, scalar(1.0)), scalar(1e-4));
    optimizerMoveLimitActive =
      Foam::min
      (
        Foam::max(optimizerMoveLimitActive*optMMAStagnationMoveLimitBoost, scalar(1e-4)),
        stagnationMoveLimitCap
      );
  }
  if(solid_area)
  {
     forAll(cells_solid, i)
     {
        xh[cells_solid[i]]=0;
     }
  }
  if(fluid_area)
  {
     forAll(cells_fluid, i)
     {
        xh[cells_fluid[i]]=1.0;
     }
  }
  xh.correctBoundaryConditions();

  const scalar continuationAlphaBase = alphamax.value();
  const scalar continuationAlphaAt63 = continuationAlphaBase*(63.0/7.0 + 1.0);
  scalar scheduledAlphaMaxValue = continuationAlphaAt63;
  if (continuationScheduleIterationActive <= 63)
  {
    scheduledAlphaMaxValue =
      continuationAlphaBase*(continuationScheduleIterationActive*continuationRampScale/7.0 + 1.0);
  }
  else if
  (
      continuationScheduleIterationActive > 100
   && !static_cast<bool>(optDisableLateStageContinuationAcceleration)
  )
  {
    const scalar continuationGrowthSteps =
      scalar(continuationScheduleIterationActive - 100);
    const scalar continuationGrowthBase =
      1.0 + 0.05*continuationRampScale;
    const scalar continuationGrowthFactor =
      Foam::pow(continuationGrowthBase, continuationGrowthSteps);
    scheduledAlphaMaxValue = continuationAlphaAt63*continuationGrowthFactor;
  }
  alphaMax.value() = Foam::min(scheduledAlphaMaxValue, alphaMAX.value());

  const scalar scheduledQu =
    continuationScheduleIterationActive >= 80
    ? 0.005 + (continuationScheduleIterationActive-80)*1e-4*continuationRampScale
    : 0.005;
  qu = scheduledQu;
  alphaMax.value()=Foam::min(alphaMax.value(),alphaMAX.value());
  qu=Foam::min(qu,0.01);
  const volScalarField xhMaterial
  (
    Foam::pow
    (
      Foam::max(xh, scalar(0.0)),
      materialInterpolationExponentActive
    )
  );
  const volScalarField xhBlendDenominator(qu + xhMaterial);
  auditDivisionDenominator
  (
    xhBlendDenominator,
    "update.H::alpha/DT denominator (qu+xhMaterial)",
    SMALL
  );
  if (static_cast<bool>(optDivisionAuditEnable))
  {
    const scalar rhoCpValue = rhoCp.value();
    if (!std::isfinite(static_cast<double>(rhoCpValue)) || Foam::mag(rhoCpValue) <= SMALL)
    {
      FatalErrorInFunction
        << "Division denominator audit failed."
        << " Context=update.H::DT denominator (rhoCp)"
        << " rhoCp=" << rhoCpValue
        << exit(FatalError);
    }
  }
  alpha=alphaMax*qu*(1-xhMaterial)/xhBlendDenominator;
  DT=(ks+(kf-ks)*xhMaterial*(1+qu)/xhBlendDenominator)/rhoCp;
}
