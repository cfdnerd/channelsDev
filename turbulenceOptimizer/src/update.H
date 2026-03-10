//***********************************//update mu and lambda
opt++;
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
if (turbulenceProModeActive)
{
  continuationIteration = Foam::max(opt, label(1));
  continuationScheduleIterationActive = continuationIteration;
  optimizerMoveLimitActive =
    continuationGovernorStressActive
    ? Foam::max(movlim*governorMoveLimitScaleClamped, scalar(1e-4))
    : movlim;
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
  alpha=alphaMax*qu*(1-xh)/(qu+xh);
  DT=(ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoCp;
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
  const volScalarField xhBlendDenominator(qu + xh);
  auditDivisionDenominator
  (
    xhBlendDenominator,
    "update.H::alpha/DT denominator (qu+xh)",
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
  alpha=alphaMax*qu*(1-xh)/xhBlendDenominator;
  DT=(ks+(kf-ks)*xh*(1+qu)/xhBlendDenominator)/rhoCp;
}
