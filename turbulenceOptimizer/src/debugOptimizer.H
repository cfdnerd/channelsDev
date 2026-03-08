scalar xLocalMin = 1e300, xLocalMax = -1e300;
scalar xpLocalMin = 1e300, xpLocalMax = -1e300;
scalar xhLocalMin = 1e300, xhLocalMax = -1e300;
scalar alphaLocalMin = 1e300, alphaLocalMax = -1e300;
scalar DTLocalMin = 1e300, DTLocalMax = -1e300;
scalar TLocalMin = 1e300, TLocalMax = -1e300;
scalar pLocalMin = 1e300, pLocalMax = -1e300;
scalar paLocalMin = 1e300, paLocalMax = -1e300;
scalar pbLocalMin = 1e300, pbLocalMax = -1e300;
scalar xVolumeSumLocal = 0;
scalar xhVolumeSumLocal = 0;
scalar grayVolumeLocal = 0;
scalar solidVolumeLocal = 0;
scalar fluidVolumeLocal = 0;
scalar maxAbsXChangeLocal = 0;
scalar l2XChangeLocal = 0;
scalar maxAbsXhChangeLocal = 0;
scalar l2XhChangeLocal = 0;
label xLowerBoundCountLocal = 0;
label xUpperBoundCountLocal = 0;
label xMoveLimitHitCountLocal = 0;
label xNonFiniteCountLocal = 0;
scalar drhoLocalMin = 1e300, drhoLocalMax = -1e300;
scalar fsensLocalMin = 1e300, fsensLocalMax = -1e300, fsensL2Local = 0;
scalar gsensPowerLocalMin = 1e300, gsensPowerLocalMax = -1e300, gsensPowerL2Local = 0;
scalar gsensVolLocalMin = 1e300, gsensVolLocalMax = -1e300, gsensVolL2Local = 0;
label fsensNonFiniteCountLocal = 0;
label gsensPowerNonFiniteCountLocal = 0;
label gsensVolNonFiniteCountLocal = 0;
scalar UMaxMagLocal = 0, UaMaxMagLocal = 0, UbMaxMagLocal = 0;
scalar UMeanMagLocal = 0, UaMeanMagLocal = 0, UbMeanMagLocal = 0;
scalar objectivePnormSumLocal = 0;
label kFloorCountLocal = 0;
label omegaFloorCountLocal = 0;
label nutZeroCountLocal = 0;
label nutCapCountLocal = 0;
label solidCellCountLocal = 0;
label grayCellCountLocal = 0;
label fluidCellCountLocal = 0;
label kFloorSolidCountLocal = 0;
label kFloorGrayCountLocal = 0;
label kFloorFluidCountLocal = 0;
label omegaFloorSolidCountLocal = 0;
label omegaFloorGrayCountLocal = 0;
label omegaFloorFluidCountLocal = 0;

for(i=0;i<n;i++)
{
  const scalar cellVolume = mesh.V()[i];
  const scalar xValue = x[i];
  const scalar xpValue = xp[i];
  const scalar xhValue = xh[i];
  const scalar alphaValue = alpha[i];
  const scalar DTValue = DT[i];
  const scalar TValue = T[i];
  const scalar pValue = p[i];
  const scalar paValue = pa[i];
  const scalar pbValue = pb[i];
  const scalar xChange = xValue - xOld[i];
  const scalar xhChange = xhValue - xhOld[i];
  const scalar absXChange = Foam::mag(xChange);
  const scalar absXhChange = Foam::mag(xhChange);
  const scalar UMag = Foam::mag(U[i]);
  const scalar UaMag = Foam::mag(Ua[i]);
  const scalar UbMag = Foam::mag(Ub[i]);
  const scalar kValue =
    static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
    ? kClosureState[i]
    : k[i];
  const scalar omegaValue =
    static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState)
    ? omegaClosureState[i]
    : omega[i];
  const scalar nutValue = nut[i];
  const scalar drhoValue = drho[i];
  const scalar fsensValue = fsensMeanT[i];
  const scalar gsensPowerValue = gsensPowerDiss[i];
  const scalar gsensVolValue = gsensVol[i];
  const bool isSolid = xhValue <= 0.05;
  const bool isFluid = xhValue >= 0.95;
  const bool isGray = !isSolid && !isFluid;

  xLocalMin = Foam::min(xLocalMin, xValue);
  xLocalMax = Foam::max(xLocalMax, xValue);
  xpLocalMin = Foam::min(xpLocalMin, xpValue);
  xpLocalMax = Foam::max(xpLocalMax, xpValue);
  xhLocalMin = Foam::min(xhLocalMin, xhValue);
  xhLocalMax = Foam::max(xhLocalMax, xhValue);
  alphaLocalMin = Foam::min(alphaLocalMin, alphaValue);
  alphaLocalMax = Foam::max(alphaLocalMax, alphaValue);
  DTLocalMin = Foam::min(DTLocalMin, DTValue);
  DTLocalMax = Foam::max(DTLocalMax, DTValue);
  TLocalMin = Foam::min(TLocalMin, TValue);
  TLocalMax = Foam::max(TLocalMax, TValue);
  pLocalMin = Foam::min(pLocalMin, pValue);
  pLocalMax = Foam::max(pLocalMax, pValue);
  paLocalMin = Foam::min(paLocalMin, paValue);
  paLocalMax = Foam::max(paLocalMax, paValue);
  pbLocalMin = Foam::min(pbLocalMin, pbValue);
  pbLocalMax = Foam::max(pbLocalMax, pbValue);
  UMaxMagLocal = Foam::max(UMaxMagLocal, UMag);
  UaMaxMagLocal = Foam::max(UaMaxMagLocal, UaMag);
  UbMaxMagLocal = Foam::max(UbMaxMagLocal, UbMag);

  xVolumeSumLocal += xValue*cellVolume;
  xhVolumeSumLocal += xhValue*cellVolume;
  UMeanMagLocal += UMag*cellVolume;
  UaMeanMagLocal += UaMag*cellVolume;
  UbMeanMagLocal += UbMag*cellVolume;

  if (isSolid)
  {
    solidCellCountLocal++;
  }
  else if (isFluid)
  {
    fluidCellCountLocal++;
  }
  else
  {
    grayCellCountLocal++;
  }

  if (customTurbulenceActive)
  {
    if (kValue <= 1.000001*kMin.value())
    {
      kFloorCountLocal++;
      if (isSolid)
      {
        kFloorSolidCountLocal++;
      }
      else if (isFluid)
      {
        kFloorFluidCountLocal++;
      }
      else
      {
        kFloorGrayCountLocal++;
      }
    }
    if (omegaValue <= 1.000001*omegaMin.value())
    {
      omegaFloorCountLocal++;
      if (isSolid)
      {
        omegaFloorSolidCountLocal++;
      }
      else if (isFluid)
      {
        omegaFloorFluidCountLocal++;
      }
      else
      {
        omegaFloorGrayCountLocal++;
      }
    }
    if (nutValue <= nutFloor.value())
    {
      nutZeroCountLocal++;
    }
    if (nutValue >= 0.999999*nutMax.value())
    {
      nutCapCountLocal++;
    }
  }

  if (isGray)
  {
    grayVolumeLocal += cellVolume;
  }
  if (isSolid)
  {
    solidVolumeLocal += cellVolume;
  }
  if (isFluid)
  {
    fluidVolumeLocal += cellVolume;
  }

  maxAbsXChangeLocal = Foam::max(maxAbsXChangeLocal, absXChange);
  l2XChangeLocal += xChange*xChange;
  maxAbsXhChangeLocal = Foam::max(maxAbsXhChangeLocal, absXhChange);
  l2XhChangeLocal += xhChange*xhChange;

  if (xValue <= 1e-8)
  {
    xLowerBoundCountLocal++;
  }
  if (xValue >= 1.0 - 1e-8)
  {
    xUpperBoundCountLocal++;
  }
  if (absXChange >= 0.999*optimizerMoveLimitActive)
  {
    xMoveLimitHitCountLocal++;
  }
  if (!std::isfinite(static_cast<double>(xValue)))
  {
    xNonFiniteCountLocal++;
  }

  drhoLocalMin = Foam::min(drhoLocalMin, drhoValue);
  drhoLocalMax = Foam::max(drhoLocalMax, drhoValue);

  fsensLocalMin = Foam::min(fsensLocalMin, fsensValue);
  fsensLocalMax = Foam::max(fsensLocalMax, fsensValue);
  fsensL2Local += fsensValue*fsensValue;
  if (!std::isfinite(static_cast<double>(fsensValue)))
  {
    fsensNonFiniteCountLocal++;
  }

  gsensPowerLocalMin = Foam::min(gsensPowerLocalMin, gsensPowerValue);
  gsensPowerLocalMax = Foam::max(gsensPowerLocalMax, gsensPowerValue);
  gsensPowerL2Local += gsensPowerValue*gsensPowerValue;
  if (!std::isfinite(static_cast<double>(gsensPowerValue)))
  {
    gsensPowerNonFiniteCountLocal++;
  }

  gsensVolLocalMin = Foam::min(gsensVolLocalMin, gsensVolValue);
  gsensVolLocalMax = Foam::max(gsensVolLocalMax, gsensVolValue);
  gsensVolL2Local += gsensVolValue*gsensVolValue;
  if (!std::isfinite(static_cast<double>(gsensVolValue)))
  {
    gsensVolNonFiniteCountLocal++;
  }

  if (objFunction == 2 && maxT > SMALL)
  {
    objectivePnormSumLocal += std::pow(TValue/maxT, Pnorm);
  }
}

reduce(xLocalMin, minOp<scalar>());
reduce(xLocalMax, maxOp<scalar>());
reduce(xpLocalMin, minOp<scalar>());
reduce(xpLocalMax, maxOp<scalar>());
reduce(xhLocalMin, minOp<scalar>());
reduce(xhLocalMax, maxOp<scalar>());
reduce(alphaLocalMin, minOp<scalar>());
reduce(alphaLocalMax, maxOp<scalar>());
reduce(DTLocalMin, minOp<scalar>());
reduce(DTLocalMax, maxOp<scalar>());
reduce(TLocalMin, minOp<scalar>());
reduce(TLocalMax, maxOp<scalar>());
reduce(pLocalMin, minOp<scalar>());
reduce(pLocalMax, maxOp<scalar>());
reduce(paLocalMin, minOp<scalar>());
reduce(paLocalMax, maxOp<scalar>());
reduce(pbLocalMin, minOp<scalar>());
reduce(pbLocalMax, maxOp<scalar>());
reduce(UMaxMagLocal, maxOp<scalar>());
reduce(UaMaxMagLocal, maxOp<scalar>());
reduce(UbMaxMagLocal, maxOp<scalar>());
reduce(xVolumeSumLocal, sumOp<scalar>());
reduce(xhVolumeSumLocal, sumOp<scalar>());
reduce(grayVolumeLocal, sumOp<scalar>());
reduce(solidVolumeLocal, sumOp<scalar>());
reduce(fluidVolumeLocal, sumOp<scalar>());
reduce(UMeanMagLocal, sumOp<scalar>());
reduce(UaMeanMagLocal, sumOp<scalar>());
reduce(UbMeanMagLocal, sumOp<scalar>());
reduce(kFloorCountLocal, sumOp<label>());
reduce(omegaFloorCountLocal, sumOp<label>());
reduce(nutZeroCountLocal, sumOp<label>());
reduce(nutCapCountLocal, sumOp<label>());
reduce(solidCellCountLocal, sumOp<label>());
reduce(grayCellCountLocal, sumOp<label>());
reduce(fluidCellCountLocal, sumOp<label>());
reduce(kFloorSolidCountLocal, sumOp<label>());
reduce(kFloorGrayCountLocal, sumOp<label>());
reduce(kFloorFluidCountLocal, sumOp<label>());
reduce(omegaFloorSolidCountLocal, sumOp<label>());
reduce(omegaFloorGrayCountLocal, sumOp<label>());
reduce(omegaFloorFluidCountLocal, sumOp<label>());
reduce(maxAbsXChangeLocal, maxOp<scalar>());
reduce(l2XChangeLocal, sumOp<scalar>());
reduce(maxAbsXhChangeLocal, maxOp<scalar>());
reduce(l2XhChangeLocal, sumOp<scalar>());
reduce(xLowerBoundCountLocal, sumOp<label>());
reduce(xUpperBoundCountLocal, sumOp<label>());
reduce(xMoveLimitHitCountLocal, sumOp<label>());
reduce(xNonFiniteCountLocal, sumOp<label>());
reduce(drhoLocalMin, minOp<scalar>());
reduce(drhoLocalMax, maxOp<scalar>());
reduce(fsensLocalMin, minOp<scalar>());
reduce(fsensLocalMax, maxOp<scalar>());
reduce(fsensL2Local, sumOp<scalar>());
reduce(fsensNonFiniteCountLocal, sumOp<label>());
reduce(gsensPowerLocalMin, minOp<scalar>());
reduce(gsensPowerLocalMax, maxOp<scalar>());
reduce(gsensPowerL2Local, sumOp<scalar>());
reduce(gsensPowerNonFiniteCountLocal, sumOp<label>());
reduce(gsensVolLocalMin, minOp<scalar>());
reduce(gsensVolLocalMax, maxOp<scalar>());
reduce(gsensVolL2Local, sumOp<scalar>());
reduce(gsensVolNonFiniteCountLocal, sumOp<label>());
reduce(objectivePnormSumLocal, sumOp<scalar>());

const scalar xMean = xVolumeSumLocal/area;
const scalar xhMean = xhVolumeSumLocal/area;
const scalar grayVolumeFraction = grayVolumeLocal/area;
const scalar solidVolumeFraction = solidVolumeLocal/area;
const scalar fluidVolumeFraction = fluidVolumeLocal/area;
const scalar l2XChange = std::sqrt(l2XChangeLocal);
const scalar l2XhChange = std::sqrt(l2XhChangeLocal);
const scalar normalizedL2XChange = l2XChange/std::sqrt(Foam::max(N, scalar(1)));
const scalar normalizedL2XhChange = l2XhChange/std::sqrt(Foam::max(N, scalar(1)));
const scalar UMeanMag = UMeanMagLocal/area;
const scalar UaMeanMag = UaMeanMagLocal/area;
const scalar UbMeanMag = UbMeanMagLocal/area;
const scalar turbulenceCellCountDenom = Foam::max(N, scalar(1));
const scalar kFloorFraction = static_cast<scalar>(kFloorCountLocal)/turbulenceCellCountDenom;
const scalar omegaFloorFraction = static_cast<scalar>(omegaFloorCountLocal)/turbulenceCellCountDenom;
const scalar nutZeroFraction = static_cast<scalar>(nutZeroCountLocal)/turbulenceCellCountDenom;
const scalar nutCapFraction = static_cast<scalar>(nutCapCountLocal)/turbulenceCellCountDenom;
const scalar solidCellCountDenom = Foam::max(static_cast<scalar>(solidCellCountLocal), scalar(1));
const scalar grayCellCountDenom = Foam::max(static_cast<scalar>(grayCellCountLocal), scalar(1));
const scalar fluidCellCountDenom = Foam::max(static_cast<scalar>(fluidCellCountLocal), scalar(1));
const scalar kFloorSolidFraction = static_cast<scalar>(kFloorSolidCountLocal)/solidCellCountDenom;
const scalar kFloorGrayFraction = static_cast<scalar>(kFloorGrayCountLocal)/grayCellCountDenom;
const scalar kFloorFluidFraction = static_cast<scalar>(kFloorFluidCountLocal)/fluidCellCountDenom;
const scalar omegaFloorSolidFraction = static_cast<scalar>(omegaFloorSolidCountLocal)/solidCellCountDenom;
const scalar omegaFloorGrayFraction = static_cast<scalar>(omegaFloorGrayCountLocal)/grayCellCountDenom;
const scalar omegaFloorFluidFraction = static_cast<scalar>(omegaFloorFluidCountLocal)/fluidCellCountDenom;
const scalar fsensL2 = std::sqrt(fsensL2Local);
const scalar gsensPowerL2 = std::sqrt(gsensPowerL2Local);
const scalar gsensVolL2 = std::sqrt(gsensVolL2Local);

scalar objectiveValue = MeanT;
word objectiveLabel("MeanT");
if (objFunction == 2)
{
  objectiveLabel = "PnormT";
  if (maxT > SMALL)
  {
    objectiveValue = maxT*std::pow(objectivePnormSumLocal, 1.0/Pnorm);
  }
  else
  {
    objectiveValue = 0;
  }
}
if (objFunction == 3 && static_cast<bool>(optThermalBoundedObjectiveEnable))
{
  objectiveLabel = "SmoothExcessT";
  objectiveValue = thermalExcessObjective;
}
const scalar objectiveDelta = hasPreviousObjectiveValue ? objectiveValue - previousObjectiveValue : 0;
const scalar powerConstraintWeightSafe = Foam::max(optPowerConstraintWeight, SMALL);
const scalar powerConstraintMargin = gx[mmaPowerConstraintIndex]/powerConstraintWeightSafe;
const scalar volumeConstraintWeightSafe = Foam::max(optVolumeConstraintWeight, SMALL);
const scalar volumeConstraintMargin = gx[mmaVolumeConstraintIndex]/volumeConstraintWeightSafe;
const scalar thermalTmaxConstraintMargin =
  static_cast<bool>(optThermalEnvelopeConstraintEnable) ? thermalTmaxKSConstraint : scalar(0.0);
const scalar thermalTminConstraintMargin =
  static_cast<bool>(optThermalEnvelopeConstraintEnable) ? thermalTminKSConstraint : scalar(0.0);
const scalar powerConstraintMarginDelta =
  hasPreviousPowerConstraintMargin ? powerConstraintMargin - previousPowerConstraintMargin : 0;
const scalar powerConstraintFinalMargin =
  PowerDiss/Foam::max(PowerDissMax, SMALL) - 1.0;
const scalar powerConstraintFinalMarginDelta =
  hasPreviousPowerConstraintFinalMargin
  ? powerConstraintFinalMargin - previousPowerConstraintFinalMargin
  : 0;
const scalar volumeConstraintMarginDelta =
  hasPreviousVolumeConstraintMargin ? volumeConstraintMargin - previousVolumeConstraintMargin : 0;
const scalar thermalTmaxConstraintMarginDelta =
  hasPreviousThermalTmaxConstraintMargin
  ? thermalTmaxConstraintMargin - previousThermalTmaxConstraintMargin
  : 0;
const scalar thermalTminConstraintMarginDelta =
  hasPreviousThermalTminConstraintMargin
  ? thermalTminConstraintMargin - previousThermalTminConstraintMargin
  : 0;

label solverSeverity(0);
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(ULinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(pLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(kLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(omegaLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(kaLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(omegaaLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(TLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(UbLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(pbLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(UaLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(paLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(xpFilterLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(fsensFilterLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(gsensPowerFilterLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(gsensVolFilterLinearSolveStats));
reduce(solverSeverity, maxOp<label>());

auto linearSolveMaxFinalResidual =
[&](const LinearSolveStats& stats) -> scalar
{
  return stats.callCount > 0 ? stats.finalResidualMax : scalar(0.0);
};
scalar maxLinearFinalResidual = 0.0;
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(ULinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(pLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(kLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(omegaLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(kaLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(omegaaLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(TLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(UbLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(pbLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(UaLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(paLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(xpFilterLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(fsensFilterLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(gsensPowerFilterLinearSolveStats));
maxLinearFinalResidual = Foam::max(maxLinearFinalResidual, linearSolveMaxFinalResidual(gsensVolFilterLinearSolveStats));
reduce(maxLinearFinalResidual, maxOp<scalar>());
const scalar linearResidualGrowthRatio =
  hasPreviousMaxLinearFinalResidual
  ? maxLinearFinalResidual/Foam::max(previousMaxLinearFinalResidual, SMALL)
  : scalar(1.0);

const scalar kGlobalMin =
  customTurbulenceActive
  ? gMin(static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState) ? kClosureState : k)
  : 0;
const scalar kGlobalMax =
  customTurbulenceActive
  ? gMax(static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState) ? kClosureState : k)
  : 0;
const scalar omegaGlobalMin =
  customTurbulenceActive
  ? gMin(static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState) ? omegaClosureState : omega)
  : 0;
const scalar omegaGlobalMax =
  customTurbulenceActive
  ? gMax(static_cast<bool>(optimizationCycleClosureConsistencyUseBoundedState) ? omegaClosureState : omega)
  : 0;
const scalar nutGlobalMax = customTurbulenceActive ? gMax(nut) : 0;

const scalar objectiveRelativeDelta = Foam::mag(objectiveDelta)/Foam::max(Foam::mag(objectiveValue), SMALL);
const bool objectiveWorsening =
  opt > 0
  && objectiveDelta > optObjectiveWorsenTol;
const scalar maxConstraintViolationCore =
  Foam::max(Foam::max(powerConstraintMargin, volumeConstraintMargin), scalar(0));
const scalar thermalConstraintViolation =
  static_cast<bool>(optThermalEnvelopeConstraintEnable)
  ? Foam::max(Foam::max(thermalTmaxConstraintMargin, thermalTminConstraintMargin), scalar(0))
  : scalar(0);
const scalar maxConstraintViolation =
  Foam::max(maxConstraintViolationCore, thermalConstraintViolation);
const scalar maxConstraintViolationDelta =
  hasPreviousMaxConstraintViolation ? maxConstraintViolation - previousMaxConstraintViolation : 0;
const scalar guardConstraintViolation =
  Foam::max
  (
    Foam::max(powerConstraintFinalMargin, volumeConstraintMargin),
    thermalConstraintViolation
  );
const scalar guardConstraintViolationDelta =
  hasPreviousGuardConstraintViolation
  ? guardConstraintViolation - previousGuardConstraintViolation
  : 0;
const bool powerConstraintMarginWorsening =
  hasPreviousPowerConstraintMargin
  && powerConstraintMargin > optConstraintTol
  && powerConstraintMarginDelta > optConstraintWorsenTol;
const bool powerConstraintFinalWorsening =
  hasPreviousPowerConstraintFinalMargin
  && powerConstraintFinalMargin > optConstraintTol
  && powerConstraintFinalMarginDelta > optConstraintWorsenTol;
const bool powerConstraintWorsening =
  powerConstraintMarginWorsening;
const bool maxConstraintViolationWorsening =
  hasPreviousMaxConstraintViolation
  && maxConstraintViolation > optConstraintTol
  && maxConstraintViolationDelta > optConstraintWorsenTol;
const bool guardConstraintViolationWorsening =
  hasPreviousGuardConstraintViolation
  && guardConstraintViolation > optConstraintTol
  && guardConstraintViolationDelta > optConstraintWorsenTol;
const bool constraintWorsenGuardEnabled = optConstraintWorsenWindow > 0;
const label scheduleIterationCurrent =
  Foam::max(opt, label(1));
const label scheduleIterationPrevious =
  Foam::max(opt - 1, label(0));
const scalar powerConstraintScheduleDenomCurrent =
  Foam::max(PowerDissRelax - scalar(0.2)*scheduleIterationCurrent, PowerDissMax);
const scalar powerConstraintScheduleDenomPrevious =
  Foam::max(PowerDissRelax - scalar(0.2)*scheduleIterationPrevious, PowerDissMax);
const bool powerConstraintScheduleSettled =
  powerConstraintScheduleDenomCurrent <= PowerDissMax + SMALL;
const bool powerConstraintScheduleTightening =
  powerConstraintScheduleDenomCurrent + SMALL < powerConstraintScheduleDenomPrevious;
const bool constraintWorsenScheduleReady =
  powerConstraintScheduleSettled && !powerConstraintScheduleTightening;
const bool constraintWorsenScheduleReadyTransition =
  constraintWorsenScheduleReady
  && hasPreviousConstraintWorsenScheduleReady
  && !previousConstraintWorsenScheduleReady;
const scalar powerConstraintScheduleStep = 0.2;
const scalar powerConstraintScheduleGap =
  Foam::max(powerConstraintScheduleDenomCurrent - PowerDissMax, scalar(0));
const label powerConstraintScheduleIterationsToSettle =
  powerConstraintScheduleGap > SMALL
  ? static_cast<label>(powerConstraintScheduleGap/powerConstraintScheduleStep + 0.999999)
  : 0;
const scalar optimizationIterationsRemainingRaw =
  Foam::max(runTime.endTime().value() - static_cast<scalar>(opt), scalar(0));
const label optimizationIterationsRemaining =
  static_cast<label>(optimizationIterationsRemainingRaw + 0.000001);
const label scheduleHorizonLagIterations =
  Foam::max(powerConstraintScheduleIterationsToSettle - optimizationIterationsRemaining, label(0));
const label constraintScheduleHorizonLagTolerance =
  Foam::max(optConstraintScheduleHorizonLagTolerance, label(0));
const bool powerConstraintScheduleWillSettleInRun =
  powerConstraintScheduleIterationsToSettle <= optimizationIterationsRemaining;
const bool continuationFreezePhaseActive =
  powerConstraintScheduleSettled || !powerConstraintScheduleWillSettleInRun;
const scalar convergenceConstraintViolation =
  powerConstraintScheduleSettled ? maxConstraintViolation : guardConstraintViolation;
const scalar scheduleHorizonConstraintViolation =
  powerConstraintScheduleSettled ? maxConstraintViolation : guardConstraintViolation;
const bool constraintWorsenCountingActive =
  constraintWorsenGuardEnabled
  && scheduleIterationCurrent >= optConstraintWorsenMinIterations
  && constraintWorsenScheduleReady;
if (constraintWorsenCountingActive && !hasGuardConstraintViolationAtCountStart)
{
  guardConstraintViolationAtCountStart = guardConstraintViolation;
  guardConstraintViolationMinSinceCountStart = guardConstraintViolation;
  hasGuardConstraintViolationAtCountStart = true;
  guardConstraintViolationCountStartIter = opt;
}
else if (constraintWorsenCountingActive)
{
  guardConstraintViolationMinSinceCountStart =
    Foam::min(guardConstraintViolationMinSinceCountStart, guardConstraintViolation);
}
else if (!constraintWorsenCountingActive)
{
  guardConstraintViolationAtCountStart = 0;
  guardConstraintViolationMinSinceCountStart = 0;
  hasGuardConstraintViolationAtCountStart = false;
  guardConstraintViolationCountStartIter = -1;
}
const scalar guardConstraintViolationNetRise =
  hasGuardConstraintViolationAtCountStart
  ? guardConstraintViolation - guardConstraintViolationAtCountStart
  : 0;
const scalar guardConstraintViolationReboundRise =
  hasGuardConstraintViolationAtCountStart
  ? guardConstraintViolation - guardConstraintViolationMinSinceCountStart
  : 0;
const label guardConstraintViolationCountAge =
  hasGuardConstraintViolationAtCountStart && guardConstraintViolationCountStartIter >= 0
  ? Foam::max(opt - guardConstraintViolationCountStartIter, 0)
  : 0;
const label maxGuardConstraintViolationAgeInRun =
  guardConstraintViolationCountAge + optimizationIterationsRemaining;
const scalar effectiveConstraintWorsenNetRiseTol =
  Foam::max
  (
    optConstraintWorsenNetRiseTol,
    optConstraintWorsenTol*Foam::max(optConstraintWorsenWindow, label(1))
  );
const label effectiveConstraintWorsenNetRiseMinAge =
  Foam::max
  (
    label(1),
    Foam::min
    (
      Foam::min
      (
        optConstraintWorsenNetRiseMinAge,
        Foam::max(optConstraintWorsenWindow, label(1))
      ),
      Foam::max(maxGuardConstraintViolationAgeInRun, label(1))
    )
  );
const label effectiveConstraintReboundMinAge =
  Foam::max
  (
    label(1),
    Foam::min
    (
      Foam::max
      (
        Foam::max(optConstraintReboundMinAge, label(1)),
        effectiveConstraintWorsenNetRiseMinAge
      ),
      Foam::max(maxGuardConstraintViolationAgeInRun, label(1))
    )
  );
const label effectiveConstraintReboundWindow =
  Foam::max(optConstraintReboundWindow, label(1));
const scalar effectiveConstraintReboundRiseTol =
  Foam::max
  (
    optConstraintReboundRiseTol,
    optConstraintWorsenTol*Foam::max(effectiveConstraintReboundWindow, label(1))
  );
const label candidatePowerConstraintWorsenStreak =
  (constraintWorsenCountingActive && powerConstraintFinalWorsening)
  ? powerConstraintWorsenStreak + 1
  : 0;
const label candidateMaxViolationWorsenStreak =
  (constraintWorsenCountingActive && guardConstraintViolationWorsening)
  ? maxViolationWorsenStreak + 1
  : 0;
const label candidateObjectiveWorsenStreak =
  (
      constraintWorsenCountingActive
   && maxConstraintViolationWorsening
   && objectiveWorsening
  )
  ? objectiveWorsenStreak + 1
  : 0;
const bool constraintWorsenGuardTriggered =
  constraintWorsenCountingActive
  && candidateMaxViolationWorsenStreak >= optConstraintWorsenWindow;
const bool constraintWorsenNetRiseTriggered =
  constraintWorsenCountingActive
  && hasGuardConstraintViolationAtCountStart
  && guardConstraintViolationCountAge >= effectiveConstraintWorsenNetRiseMinAge
  && guardConstraintViolation > optConstraintTol
  && guardConstraintViolationDelta > optConstraintWorsenTol
  && guardConstraintViolationNetRise >= effectiveConstraintWorsenNetRiseTol;
const bool constraintWorsenReboundTriggered =
  constraintWorsenCountingActive
  && optConstraintReboundWindow > 0
  && hasGuardConstraintViolationAtCountStart
  && guardConstraintViolationCountAge >= effectiveConstraintReboundMinAge
  && candidateMaxViolationWorsenStreak >= effectiveConstraintReboundWindow
  && guardConstraintViolation > optConstraintTol
  && guardConstraintViolationDelta > optConstraintWorsenTol
  && guardConstraintViolationReboundRise >= effectiveConstraintReboundRiseTol;
const scalar constraintReboundNetRiseCorroborationTol =
  0.5*effectiveConstraintWorsenNetRiseTol;
const bool constraintWorsenReboundCorroborated =
  objectiveWorsening
  || candidateMaxViolationWorsenStreak > effectiveConstraintReboundWindow
  || guardConstraintViolationNetRise >= constraintReboundNetRiseCorroborationTol;
const bool constraintWorsenReboundStopTriggered =
  constraintWorsenReboundTriggered && constraintWorsenReboundCorroborated;
const bool constraintScheduleHorizonGuardTriggered =
  static_cast<bool>(optConstraintScheduleHorizonGuard)
  && constraintWorsenGuardEnabled
  && opt >= optConstraintWorsenMinIterations
  && !constraintWorsenScheduleReady
  && scheduleHorizonLagIterations > constraintScheduleHorizonLagTolerance
  && scheduleHorizonConstraintViolation > optConstraintTol;
const bool objectiveWorsenGuardTriggered =
  optObjectiveWorsenWindow > 0
  && candidateObjectiveWorsenStreak >= optObjectiveWorsenWindow;
const scalar turbulenceClipObsFraction = 0.25;
const scalar turbulenceClipWarnFraction = 0.70;
const scalar maxClipFraction = Foam::max
(
  Foam::max(kFloorFraction, omegaFloorFraction),
  Foam::max(nutZeroFraction, nutCapFraction)
);
const scalar maxFluidClipFraction = Foam::max(kFloorFluidFraction, omegaFloorFluidFraction);
const scalar maxGrayClipFraction = Foam::max(kFloorGrayFraction, omegaFloorGrayFraction);
const scalar maxTransitionClipFraction = Foam::max(maxGrayClipFraction, maxFluidClipFraction);
const scalar clipGuardFluidFractionThreshold = scalar(1.0);
const label clipGuardWindowDisabled(0);
const scalar clipGuardNoProgressTolClamped = scalar(0.0);
const scalar clipGuardNoProgressTolReference(5e-4);
const bool clipGuardEnabled = false;
const bool clipGuardActive = false;
const label candidateClipGuardStreak = 0;
const bool clipGuardNoProgressActive = false;
const label candidateClipGuardNoProgressStreak = 0;
const bool clipGuardTriggered = false;
const bool clipGuardRejectMode = false;
const bool clipGuardForcedReject = false;
const scalar clipHardRejectFluidFractionClamped = scalar(1.0);
const scalar clipHardRejectGrayFractionClamped = scalar(1.0);
const scalar clipHardRejectFluidJumpClamped = GREAT;
const bool clipHardRejectByFraction = false;
const bool clipHardRejectByJump = false;
const bool clipHardForcedReject = false;
const bool clipForcedReject = false;
const scalar clipMeritExcess = scalar(0.0);
const scalar clipPenaltyPreScheduleScale = scalar(1.0);
const scalar clipPenaltyRejectScale = scalar(1.0);
const scalar meritClipPenaltyScheduleScale = scalar(1.0);
const scalar meritClipPenaltyRejectScale = scalar(1.0);
const scalar meritClipPenaltyActive = scalar(0.0);
const scalar powerViolationPositive =
  Foam::max(powerConstraintFinalMargin, scalar(0.0));
const scalar volumeViolationPositive =
  Foam::max(volumeConstraintMargin, scalar(0.0));
const scalar thermalTmaxViolationPositive =
  static_cast<bool>(optThermalEnvelopeConstraintEnable)
  ? Foam::max(thermalTmaxConstraintMargin, scalar(0.0))
  : scalar(0.0);
const scalar thermalTminViolationPositive =
  static_cast<bool>(optThermalEnvelopeConstraintEnable)
  ? Foam::max(thermalTminConstraintMargin, scalar(0.0))
  : scalar(0.0);
const scalar thermalEnvelopeViolation =
  Foam::max(thermalTmaxViolationPositive, thermalTminViolationPositive);
const scalar meritConstraintTermLinear =
  scalar(0.0)*guardConstraintViolation;
const scalar meritConstraintTerm = meritConstraintTermLinear;
const scalar meritThermalEnvelopeTerm =
  static_cast<bool>(optThermalEnvelopeConstraintEnable)
  ? optThermalEnvelopePenalty*thermalEnvelopeViolation
  : scalar(0.0);
const scalar meritValue =
  objectiveValue
  + meritConstraintTerm
  + meritThermalEnvelopeTerm
  + meritClipPenaltyActive*clipMeritExcess;
const scalar meritAcceptanceTol =
  hasPreviousAcceptedMerit
  ? (scalar(0.0) + scalar(0.0)*Foam::mag(previousAcceptedMerit))
  : 0;
const scalar rejectMoveScaleClamped =
  Foam::max(Foam::min(scalar(1.0), scalar(0.99)), scalar(0.1));
const scalar recoverMoveScaleClamped =
  Foam::max(scalar(1.0), scalar(1.0));
(void)recoverMoveScaleClamped;
const scalar meritTrustMoveMinClamped =
  Foam::max(Foam::min(scalar(1.0), scalar(1.0)), scalar(0.01));
const scalar continuationLowTrustThresholdReport =
  Foam::max(Foam::min(optContinuationLowTrustThreshold, scalar(1.0)), meritTrustMoveMinClamped);
const bool continuationFrozenAny =
  continuationFrozenByReject || continuationFrozenByLowTrust || continuationFrozenByClipHold;

scalar predictedObjectiveDeltaLocal = 0;
scalar predictedPowerDeltaLocal = 0;
scalar predictedVolumeDeltaLocal = 0;
scalar predictedThermalTmaxDeltaLocal = 0;
scalar predictedThermalTminDeltaLocal = 0;
for (i=0; i<n; i++)
{
  const scalar dx = x[i] - xOld[i];
  predictedObjectiveDeltaLocal += static_cast<scalar>(dfdx[i])*dx;
  predictedPowerDeltaLocal += static_cast<scalar>(dgdx[mmaPowerConstraintIndex][i])*dx;
  predictedVolumeDeltaLocal += static_cast<scalar>(dgdx[mmaVolumeConstraintIndex][i])*dx;
  if (mmaUseThermalEnvelopeConstraints && m > mmaThermalTminConstraintIndex)
  {
    predictedThermalTmaxDeltaLocal +=
      static_cast<scalar>(dgdx[mmaThermalTmaxConstraintIndex][i])*dx;
    predictedThermalTminDeltaLocal +=
      static_cast<scalar>(dgdx[mmaThermalTminConstraintIndex][i])*dx;
  }
}
reduce(predictedObjectiveDeltaLocal, sumOp<scalar>());
reduce(predictedPowerDeltaLocal, sumOp<scalar>());
reduce(predictedVolumeDeltaLocal, sumOp<scalar>());
reduce(predictedThermalTmaxDeltaLocal, sumOp<scalar>());
reduce(predictedThermalTminDeltaLocal, sumOp<scalar>());
const scalar predictedObjectiveDelta = predictedObjectiveDeltaLocal;
const scalar predictedPowerDelta = predictedPowerDeltaLocal/powerConstraintWeightSafe;
const scalar predictedVolumeDelta = predictedVolumeDeltaLocal/volumeConstraintWeightSafe;
const scalar predictedThermalTmaxDelta = predictedThermalTmaxDeltaLocal;
const scalar predictedThermalTminDelta = predictedThermalTminDeltaLocal;
const scalar predictedThermalTmaxViolationPositive =
  Foam::max(thermalTmaxViolationPositive + predictedThermalTmaxDelta, scalar(0.0));
const scalar predictedThermalTminViolationPositive =
  Foam::max(thermalTminViolationPositive + predictedThermalTminDelta, scalar(0.0));
const scalar predictedGuardConstraintDelta =
  Foam::max
  (
    Foam::max(predictedPowerDelta, predictedVolumeDelta),
    Foam::max(predictedThermalTmaxDelta, predictedThermalTminDelta)
  );
const scalar predictedConstraintMeritDeltaLinear =
  scalar(0.0)*predictedGuardConstraintDelta;
const scalar predictedConstraintMeritDelta = predictedConstraintMeritDeltaLinear;
const scalar predictedThermalEnvelopeMeritDelta =
  static_cast<bool>(optThermalEnvelopeConstraintEnable)
  ? optThermalEnvelopePenalty
    *(
       Foam::max(predictedThermalTmaxViolationPositive, predictedThermalTminViolationPositive)
      - thermalEnvelopeViolation
     )
  : scalar(0.0);
const scalar predictedMeritDelta =
  predictedObjectiveDelta + predictedConstraintMeritDelta + predictedThermalEnvelopeMeritDelta;
const scalar predictedMeritReduction =
  -predictedMeritDelta;
const scalar actualMeritReduction =
  hasPreviousAcceptedMerit ? previousAcceptedMerit - meritValue : 0.0;
const bool trustRatioDefined =
  hasPreviousAcceptedMerit && predictedMeritReduction > SMALL;
const scalar trustRatio =
  trustRatioDefined
  ? actualMeritReduction/Foam::max(predictedMeritReduction, SMALL)
  : 1.0;
const bool trustRatioAccepted =
  !static_cast<bool>(false)
  || !hasPreviousAcceptedMerit
  || !trustRatioDefined
  || trustRatio >= scalar(0.0);
const bool trustRatioAcceptedWithSensitivityGuard =
  trustRatioAccepted
  && (
       !static_cast<bool>(optSensitivityQualityGateEnable)
    || sensitivityQualityScore > optSensitivityQualityWarn
    || !trustRatioDefined
    || trustRatio >= Foam::max(scalar(0.0), scalar(0.2))
     );
(void)trustRatioAcceptedWithSensitivityGuard;
const bool filterAcceptanceEnabled =
  static_cast<bool>(false)
  && optimizationCycleTurbulentOptimizationMode
  && label(0) > 0;
bool filterDominated = false;
if (filterAcceptanceEnabled)
{
  for (std::size_t fi = 0; fi < filterConstraintArchive.size(); ++fi)
  {
    if
    (
         guardConstraintViolation >= filterConstraintArchive[fi] - scalar(0.0)
      && objectiveValue >= filterObjectiveArchive[fi] - scalar(0.0)
    )
    {
      filterDominated = true;
      break;
    }
  }
}
const bool filterAccepted = !filterAcceptanceEnabled || !filterDominated;
const label filterDeadlockRejectThreshold =
  Foam::max(label(0), label(4));
const bool filterDeadlockBypass =
(
     filterAcceptanceEnabled
  && !optDiagnosticMode
  && !filterAccepted
  && hasPreviousAcceptedMerit
  && static_cast<bool>(false)
  && meritRejectStreak >= filterDeadlockRejectThreshold
  && constraintWorsenScheduleReady
  && !clipForcedReject
  && !objectiveWorsening
  && guardConstraintViolationDelta <= optConstraintWorsenTol
  && meritValue <= previousAcceptedMerit + scalar(2.0)*meritAcceptanceTol
)
||
(
     filterAcceptanceEnabled
  && !optDiagnosticMode
  && !filterAccepted
  && hasPreviousAcceptedMerit
  && static_cast<bool>(false)
  && meritRejectStreak >= 3*filterDeadlockRejectThreshold
  && meritTrustMoveScale <= scalar(1.05)*meritTrustMoveMinClamped
  && !clipForcedReject
  && guardConstraintViolationDelta
     < -Foam::max(scalar(0.05)*Foam::max(optConstraintTol, SMALL), scalar(1e-5))
  && (
       !static_cast<bool>(optThermalEnvelopeConstraintEnable)
    || thermalTmaxConstraintMarginDelta <= scalar(0.0)
     )
);
const bool filterAcceptedForMerit =
  filterAccepted || filterDeadlockBypass;
const bool meritAcceptedPreGate =
  !static_cast<bool>(false)
  || !hasPreviousAcceptedMerit
  || meritValue <= previousAcceptedMerit + meritAcceptanceTol;
(void)meritAcceptedPreGate;
const scalar objectiveAbsValue = Foam::mag(objectiveValue);
const scalar thermalFieldAbsMaxObserved =
  Foam::max(Foam::mag(TLocalMin), Foam::mag(TLocalMax));
const bool objectiveEnvelopeViolated =
     !std::isfinite(static_cast<double>(objectiveValue))
  || objectiveAbsValue > Foam::max(optObjectiveAbsMax, scalar(0.0));
const bool thermalFieldEnvelopeViolated =
     !std::isfinite(static_cast<double>(TLocalMin))
  || !std::isfinite(static_cast<double>(TLocalMax))
  || thermalFieldAbsMaxObserved > Foam::max(optThermalFieldAbsMax, scalar(0.0));
const bool fieldEnvelopeForcedReject =
  optimizationCycleTurbulentOptimizationMode
  && (objectiveEnvelopeViolated || thermalFieldEnvelopeViolated);
const bool constraintNonDeteriorationAccepted =
     !hasPreviousGuardConstraintViolation
  || guardConstraintViolation <= previousGuardConstraintViolation + optConstraintWorsenTol;
const bool clipJumpHealthGatePassed =
     !hasPreviousClipDiagnostics
  || !clipHardRejectByJump;
const scalar residualGrowthAcceptLimit =
  optDiagnosticMode
  ? scalar(1.2)
  : scalar(2.0);
const bool residualGrowthGatePassed =
     !hasPreviousMaxLinearFinalResidual
  || maxLinearFinalResidual <= linearSolveWarnResidual
  || linearResidualGrowthRatio <= residualGrowthAcceptLimit;
const bool healthGatesPassed =
     !fieldEnvelopeForcedReject
  && clipJumpHealthGatePassed
  && residualGrowthGatePassed;
const bool sensitivityQualityGateRejectEnabled =
     static_cast<bool>(optSensitivityQualityGateEnable)
  && (
       optDiagnosticMode
    || !static_cast<bool>(optProductionIgnoreSensitivityReject)
     );
const bool sensitivityQualityGatePassed =
     !sensitivityQualityGateRejectEnabled
  || sensitivityQualityScore > optSensitivityQualityFail;
const bool sensitivityLowConfidence =
  static_cast<bool>(optSensitivityQualityGateEnable)
  && sensitivityQualityScore <= optSensitivityQualityWarn;
const scalar feasibilityRestorationClipSlack =
  Foam::max(scalar(0.0), scalar(0.0));
const scalar feasibilityRestorationGrayClipLimit =
  Foam::max(Foam::min(optRogueClipGrayThreshold + feasibilityRestorationClipSlack, scalar(1.0)), scalar(0.0));
const bool feasibilityRestorationClipCompatible =
     !clipForcedReject
  || (
       static_cast<bool>(false)
    && !clipHardForcedReject
    && (
         !hasPreviousClipDiagnostics
      || maxFluidClipFraction <= previousMaxFluidClipFraction + feasibilityRestorationClipSlack
       )
    && maxGrayClipFraction <= feasibilityRestorationGrayClipLimit
    && guardConstraintViolationDelta < 0
     );
const bool feasibilityRestorationEligible =
     static_cast<bool>(false)
  && !optDiagnosticMode
  && optimizationCycleTurbulentOptimizationMode
  && hasPreviousGuardConstraintViolation
  && meritRejectStreak >= Foam::max(label(0), label(1))
  && feasibilityRestorationClipCompatible;
const label feasibilityRestorationRejectExcess =
  Foam::max
  (
    meritRejectStreak - Foam::max(label(0), label(1)),
    label(0)
  );
const scalar feasibilityRestorationObjectiveRelTolAdaptive =
  Foam::min
  (
    scalar(0.60),
    Foam::max
    (
      scalar(0.0),
      scalar(0.05) + scalar(0.02)*feasibilityRestorationRejectExcess
    )
  );
const scalar feasibilityRestorationObjectiveAbsTolAdaptive =
  Foam::max
  (
    scalar(0.0),
    scalar(0.5)*meritAcceptanceTol + scalar(2e-3)*feasibilityRestorationRejectExcess
  );
const bool feasibilityRestorationConstraintImproved =
  feasibilityRestorationEligible
  && (
      guardConstraintViolation
     <= previousGuardConstraintViolation
      - Foam::max
        (
          scalar(0.0),
          scalar(0.1)*optConstraintWorsenTol
        )
     );
const bool feasibilityRestorationObjectiveAcceptable =
  feasibilityRestorationEligible
  && (
       objectiveDelta <= Foam::max(feasibilityRestorationObjectiveAbsTolAdaptive, optObjectiveWorsenTol)
    || (
         hasPreviousObjectiveValue
      && objectiveDelta > 0
      && objectiveDelta/Foam::max(Foam::mag(previousObjectiveValue), SMALL)
         <= Foam::max(feasibilityRestorationObjectiveRelTolAdaptive, scalar(0.0))
       )
     );
const bool feasibilityRestorationAccepted =
  feasibilityRestorationConstraintImproved
  && feasibilityRestorationObjectiveAcceptable;
const bool clipHeavyState =
     maxFluidClipFraction >= clipGuardFluidFractionThreshold
  || maxGrayClipFraction >= optRogueClipGrayThreshold;
const bool restorationStrictlyImproving =
     objectiveDelta <= 0.0
  && guardConstraintViolationDelta < 0.0
  && healthGatesPassed;
const bool feasibilityRestorationAcceptedGuarded =
  feasibilityRestorationAccepted
  && (!clipHeavyState || restorationStrictlyImproving);
(void)feasibilityRestorationAcceptedGuarded;
const bool deadlockConstraintPriorityEligible =
     feasibilityRestorationEligible
  && hasPreviousAcceptedMerit
  && static_cast<bool>(false)
  && meritRejectStreak >= filterDeadlockRejectThreshold
  && meritTrustMoveScale <= scalar(1.05)*meritTrustMoveMinClamped
  && !clipForcedReject
  && !fieldEnvelopeForcedReject
  && healthGatesPassed
  && sensitivityQualityGatePassed
  && (filterAcceptedForMerit || meritRejectStreak >= 2*filterDeadlockRejectThreshold)
  && maxConstraintViolation > optConstraintTol;
(void)deadlockConstraintPriorityEligible;
const bool deadlockConstraintDominatedByVolume =
  volumeViolationPositive >= powerViolationPositive
  && volumeViolationPositive >= thermalEnvelopeViolation;
const bool deadlockConstraintProgress =
     guardConstraintViolation
     <= previousGuardConstraintViolation
       - Foam::max
         (
           scalar(0.0),
           scalar(0.1)*optConstraintWorsenTol
         )
  || volumeConstraintMarginDelta
     <= -Foam::max
        (
          scalar(0.5)*scalar(0.0),
          scalar(0.1)*Foam::max(optConstraintTol, SMALL)
        )
  || (
       deadlockConstraintDominatedByVolume
    && predictedVolumeDelta
       <= -Foam::max
          (
            scalar(0.5)*optConstraintTol,
            scalar(0.0)
          )
     );
(void)deadlockConstraintProgress;
const scalar deadlockObjectiveAbsTol =
  Foam::max
  (
    scalar(5.0)*feasibilityRestorationObjectiveAbsTolAdaptive,
    scalar(2.0)*meritAcceptanceTol
  );
const scalar deadlockObjectiveRelTol =
  Foam::min(scalar(0.75), scalar(2.0)*feasibilityRestorationObjectiveRelTolAdaptive);
const bool deadlockObjectiveAcceptable =
     objectiveDelta <= Foam::max(deadlockObjectiveAbsTol, optObjectiveWorsenTol)
  || (
       hasPreviousObjectiveValue
    && objectiveDelta > 0
    && objectiveDelta/Foam::max(Foam::mag(previousObjectiveValue), SMALL)
       <= Foam::max(deadlockObjectiveRelTol, scalar(0.0))
     );
(void)deadlockObjectiveAcceptable;
// Step-acceptance controller retired: always advance MMA design updates unless
// a hard safety envelope is violated.
const bool globalizationFrameworkEnabled(false);
(void)globalizationFrameworkEnabled;
const bool meritAccepted =
     !fieldEnvelopeForcedReject
  && !clipHardForcedReject
  && healthGatesPassed
  && sensitivityQualityGatePassed;
const bool meritDeteriorating =
  false;
const bool clipRecoverySignal =
  customTurbulenceActive
  && (
       maxFluidClipFraction >= 0.9*clipGuardFluidFractionThreshold
    || maxGrayClipFraction >= 0.85
  );
(void)clipRecoverySignal;
const bool meritRecoverySignal = false;
(void)meritRecoverySignal;
meritRecoveryRequest = false;
const label rogueRejectTriggerClamped =
  Foam::max(optRogueRejectStreakTrigger, label(1));
const label rogueGraceIterationsClamped =
  Foam::max(optRogueGraceIterations, label(0));
const label rogueDirectionSwitchCooldownClamped =
  Foam::max(optRogueDirectionSwitchCooldown, label(1));
const label rogueDirectionSwitchMaxCountClamped =
  Foam::max(optRogueDirectionSwitchMaxCount, label(0));
const scalar rogueClipFluidThresholdClamped =
  Foam::max(Foam::min(optRogueClipFluidThreshold, scalar(1.0)), scalar(0.0));
const scalar rogueClipGrayThresholdClamped =
  Foam::max(Foam::min(optRogueClipGrayThreshold, scalar(1.0)), scalar(0.0));
const scalar rogueDirectionSwitchStepScaleClamped =
  Foam::max(optRogueDirectionSwitchStepScale, scalar(0.0));
const scalar rogueDirectionSwitchConstraintBlendClamped =
  Foam::max(Foam::min(optRogueDirectionSwitchConstraintBlend, scalar(1.0)), scalar(0.0));
const scalar rogueDirectionSwitchNoiseScaleClamped =
  Foam::max(optRogueDirectionSwitchNoiseScale, scalar(0.0));
const scalar rogueDirectionSwitchClipDampFluidFractionClamped =
  Foam::max(Foam::min(optRogueDirectionSwitchClipDampFluidFraction, scalar(1.0)), scalar(0.0));
const scalar rogueDirectionSwitchClipDampGrayFractionClamped =
  Foam::max(Foam::min(optRogueDirectionSwitchClipDampGrayFraction, scalar(1.0)), scalar(0.0));
const scalar rogueDirectionSwitchClipDampScaleClamped =
  Foam::max(Foam::min(optRogueDirectionSwitchClipDampScale, scalar(1.0)), scalar(0.0));
const scalar rogueDirectionSwitchSevereClipFluidFractionClamped =
  Foam::max(Foam::min(optRogueDirectionSwitchSevereClipFluidFraction, scalar(1.0)), scalar(0.0));
const scalar rogueDirectionSwitchSevereClipGrayFractionClamped =
  Foam::max(Foam::min(optRogueDirectionSwitchSevereClipGrayFraction, scalar(1.0)), scalar(0.0));
const scalar rogueDirectionSwitchSevereClipScaleClamped =
  Foam::max(Foam::min(optRogueDirectionSwitchSevereClipScale, scalar(1.0)), scalar(0.0));
const label candidateMeritRejectStreak =
  meritAccepted ? 0 : meritRejectStreak + 1;
const label rogueClipRejectGate =
  Foam::max(rogueRejectTriggerClamped - 1, label(1));
const bool rogueRejectSignal =
  !meritAccepted
  && candidateMeritRejectStreak >= rogueRejectTriggerClamped;
const bool rogueClipSignal =
  customTurbulenceActive
  && candidateMeritRejectStreak >= rogueClipRejectGate
  && maxFluidClipFraction >= rogueClipFluidThresholdClamped
  && maxGrayClipFraction >= rogueClipGrayThresholdClamped;
const bool rogueSignal =
  optimizationCycleTurbulentOptimizationMode
  && !meritAccepted
  && (rogueRejectSignal || rogueClipSignal);
bool rogueDirectionGraceActive = false;
bool rogueDirectionVerified = false;
bool rogueDirectionSwitchTriggered = false;
bool rogueDirectionWatchStarted = false;
bool rogueDirectionWatchGraceAdvanced = false;
bool rogueDirectionWatchCleared = false;
label rogueDirectionWatchStartIterationReport = rogueDirectionWatchStartIteration;
if
(
    optimizationCycleTurbulentOptimizationMode
 && static_cast<bool>(optRogueDirectionSwitchEnable)
 && !optDiagnosticMode
)
{
  if (rogueSignal)
  {
    if (!rogueDirectionWatchActive)
    {
      rogueDirectionWatchActive = true;
      rogueDirectionWatchGraceRemaining = rogueGraceIterationsClamped;
      rogueDirectionWatchStartIteration = opt;
      rogueDirectionWatchStartIterationReport = rogueDirectionWatchStartIteration;
      rogueDirectionWatchStarted = true;
    }
    else if (rogueDirectionWatchGraceRemaining > 0)
    {
      rogueDirectionWatchGraceRemaining -= 1;
      rogueDirectionWatchStartIterationReport = rogueDirectionWatchStartIteration;
      rogueDirectionWatchGraceAdvanced = true;
    }
    rogueDirectionGraceActive = rogueDirectionWatchGraceRemaining > 0;
    rogueDirectionVerified = !rogueDirectionGraceActive;
  }
  else
  {
    rogueDirectionWatchCleared = rogueDirectionWatchActive;
    rogueDirectionWatchActive = false;
    rogueDirectionWatchGraceRemaining = 0;
    rogueDirectionWatchStartIteration = -1;
  }

  const bool rogueDirectionSwitchCooldownSatisfied =
    opt - rogueDirectionLastSwitchIteration >= rogueDirectionSwitchCooldownClamped;
  const bool rogueDirectionSwitchBudgetAvailable =
    rogueDirectionSwitchCount < rogueDirectionSwitchMaxCountClamped;
  rogueDirectionSwitchTriggered =
    rogueDirectionVerified
    && !meritAccepted
    && rogueDirectionSwitchCooldownSatisfied
    && rogueDirectionSwitchBudgetAvailable;
}
else
{
  rogueDirectionWatchCleared = rogueDirectionWatchActive;
  rogueDirectionWatchActive = false;
  rogueDirectionWatchGraceRemaining = 0;
  rogueDirectionWatchStartIteration = -1;
}
const bool gradientCheckActive =
  static_cast<bool>(optGradientCheckEnable)
  && optimizationCycleTurbulentOptimizationMode
  && opt >= optGradientCheckStartIter
  && l2XChange >= optGradientCheckMinStepL2;
const scalar objectiveConsistencyRatio =
  gradientCheckActive
  ? Foam::mag(objectiveDelta - predictedObjectiveDelta)
    / Foam::max(Foam::max(Foam::mag(objectiveDelta), Foam::mag(predictedObjectiveDelta)), SMALL)
  : 0;
const scalar powerConsistencyRatio =
  gradientCheckActive
  ? Foam::mag(powerConstraintMarginDelta - predictedPowerDelta)
    / Foam::max(Foam::max(Foam::mag(powerConstraintMarginDelta), Foam::mag(predictedPowerDelta)), SMALL)
  : 0;
const scalar volumeConsistencyRatio =
  gradientCheckActive
  ? Foam::mag(volumeConstraintMarginDelta - predictedVolumeDelta)
    / Foam::max(Foam::max(Foam::mag(volumeConstraintMarginDelta), Foam::mag(predictedVolumeDelta)), SMALL)
  : 0;
const scalar gradientConsistencyRatio =
  Foam::max(objectiveConsistencyRatio, Foam::max(powerConsistencyRatio, volumeConsistencyRatio));
const bool gradientCheckAcceptedStep = gradientCheckActive && meritAccepted;
const bool gradientCheckWarnTriggered =
  gradientCheckAcceptedStep && gradientConsistencyRatio > optGradientCheckWarnRatio;
const bool gradientCheckFailTriggered =
  gradientCheckAcceptedStep && gradientConsistencyRatio > optGradientCheckFailRatio;
const label candidateGradientCheckFailStreak =
  gradientCheckFailTriggered ? gradientCheckFailStreak + 1 : 0;
const bool gradientCheckStopTriggered =
  static_cast<bool>(optGradientCheckStopOnFailure)
  && candidateGradientCheckFailStreak >= optGradientCheckFailWindow;

closureDerivativeFDStageCheckMeanRatio = 0.0;
closureDerivativeFDStageCheckMaxRatio = 0.0;
closureDerivativeFDStageCheckMaxRatioK = 0.0;
closureDerivativeFDStageCheckMaxRatioOmega = 0.0;
closureDerivativeFDStageCheckMaxStoredMismatch = 0.0;
closureDerivativeFDStageCheckSamplesUsed = 0;
scalar closureDerivativeFDStageCheckWarnFractionObserved = 0.0;
scalar closureDerivativeFDStageCheckFailFractionObserved = 0.0;
bool closureDerivativeFDStageCheckActive = false;

if
(
    customTurbulenceActive
 && optimizationCycleTurbulentOptimizationMode
 && static_cast<bool>(optClosureDerivativeFDStageCheckEnable)
 && optClosureDerivativeFDStageCheckSamples > 0
)
{
  const volTensorField closureFDGradU(fvc::grad(U));
  const volScalarField closureFDS2(2.0*magSqr(symm(closureFDGradU)));
  const scalarField& kInternal = k.internalField();
  const scalarField& omegaInternal = omega.internalField();
  const scalarField& yWallInternal = yWall.internalField();
  const scalar nuValue = nu.value();
  const volScalarField closureFDSqrtS2(Foam::sqrt(closureFDS2));
  const scalarField& sqrtS2Internal = closureFDSqrtS2.internalField();
  const scalarField& nutInternal = nut.internalField();
  const scalarField& dNutDkInternal = dNutDkBlock.internalField();
  const scalarField& dNutDOmegaInternal = dNutDOmegaBlock.internalField();
  const label localCellCount = kInternal.size();
  const label sampleTargetLocal = Foam::max(optClosureDerivativeFDStageCheckSamples, label(1));
  const label sampleStride = Foam::max(localCellCount/sampleTargetLocal, label(1));
  const label sampleOffset = sampleStride > 1 ? (Pstream::myProcNo() % sampleStride) : 0;
  const scalar closureFdRelStep =
    Foam::min(Foam::max(optClosureDerivativeFDRelStep, scalar(1e-7)), scalar(0.25));
  const scalar closureFdExpCoarsePlus = Foam::exp(closureFdRelStep);
  const scalar closureFdExpCoarseMinus = 1.0/closureFdExpCoarsePlus;
  const scalar closureFdHalfRelStep = 0.5*closureFdRelStep;
  const scalar closureFdExpFinePlus = Foam::exp(closureFdHalfRelStep);
  const scalar closureFdExpFineMinus = 1.0/closureFdExpFinePlus;
  const scalar omegaLimiterEpsValue =
    Foam::max(1e-3*Foam::max(omegaMin.value(), SMALL), SMALL);
  const scalar omegaBoundEpsValue =
    Foam::max(1e-3*Foam::max(omegaMin.value(), SMALL), SMALL);
  const scalar kStepFloor = Foam::max(1e-12, SMALL);
  const scalar omegaStepFloor = Foam::max(1e-12, SMALL);

  auto evaluateLocalNut =
  [&](const scalar kValue, const scalar omegaValue, const scalar sqrtS2Value, const scalar yEffValue, const scalar nuValue) -> scalar
  {
    if (!(kValue > 0.0) || !(omegaValue > 0.0))
    {
      return GREAT;
    }
    const scalar boundedY = Foam::max(yEffValue, yMin.value());
    const scalar arg2Term1 =
      2.0*Foam::sqrt(kValue)/(betaStarCoeff.value()*omegaValue*boundedY);
    const scalar arg2Term2 =
      500.0*nuValue/(sqr(boundedY)*omegaValue);
    const scalar arg2Value = Foam::min(Foam::max(arg2Term1, arg2Term2), scalar(100.0));
    const scalar F2Value = Foam::tanh(sqr(arg2Value));
    const scalar omegaLimiterBranchKValue = a1Coeff.value()*omegaValue;
    const scalar omegaLimiterBranchSValue = b1Coeff.value()*F2Value*sqrtS2Value;
    const scalar omegaLimiterBlendDenValue =
      Foam::sqrt
      (
        sqr(omegaLimiterBranchKValue - omegaLimiterBranchSValue)
        + sqr(omegaLimiterEpsValue)
      );
    const scalar omegaLimiterValue =
      0.5
      *(
        omegaLimiterBranchKValue
        + omegaLimiterBranchSValue
        + omegaLimiterBlendDenValue
      );
    const scalar omegaBoundBlendDenValue =
      Foam::sqrt
      (
        sqr(omegaLimiterValue - SMALL)
        + sqr(omegaBoundEpsValue)
      );
    const scalar boundedOmegaLimiterValue =
      0.5*(omegaLimiterValue + SMALL + omegaBoundBlendDenValue);
    return a1Coeff.value()*kValue/boundedOmegaLimiterValue;
  };

  scalar ratioSumLocal = 0.0;
  scalar ratioMaxLocal = 0.0;
  scalar ratioKMaxLocal = 0.0;
  scalar ratioOmegaMaxLocal = 0.0;
  scalar storedMismatchMaxLocal = 0.0;
  label ratioWarnCountLocal = 0;
  label ratioFailCountLocal = 0;
  label ratioSamplesLocal = 0;

  for
  (
    label celli = sampleOffset;
    celli < localCellCount && ratioSamplesLocal < sampleTargetLocal;
    celli += sampleStride
  )
  {
    const scalar kCell = kInternal[celli];
    const scalar omegaCell = omegaInternal[celli];
    const scalar nutCell = nutInternal[celli];
    const scalar sqrtS2Cell = sqrtS2Internal[celli];
    const scalar yWallCell = yWallInternal[celli];
    const scalar nuCell = nuValue;
    const bool nonPhysicalCell =
         !(kCell > 0.0)
      || !(omegaCell > 0.0)
      || !(nutCell > 0.0)
      || !(yWallCell > 0.0)
      || kCell >= GREAT
      || omegaCell >= GREAT
      || nutCell >= GREAT;
    if (nonPhysicalCell)
    {
      continue;
    }

    const scalar kPlusCoarse = kCell*closureFdExpCoarsePlus;
    const scalar kMinusCoarse = kCell*closureFdExpCoarseMinus;
    const scalar kPlusFine = kCell*closureFdExpFinePlus;
    const scalar kMinusFine = kCell*closureFdExpFineMinus;
    const scalar kSpanCoarse = Foam::max(kPlusCoarse - kMinusCoarse, kStepFloor);
    const scalar kSpanFine = Foam::max(kPlusFine - kMinusFine, kStepFloor);

    const scalar omegaPlusCoarse = omegaCell*closureFdExpCoarsePlus;
    const scalar omegaMinusCoarse = omegaCell*closureFdExpCoarseMinus;
    const scalar omegaPlusFine = omegaCell*closureFdExpFinePlus;
    const scalar omegaMinusFine = omegaCell*closureFdExpFineMinus;
    const scalar omegaSpanCoarse =
      Foam::max(omegaPlusCoarse - omegaMinusCoarse, omegaStepFloor);
    const scalar omegaSpanFine =
      Foam::max(omegaPlusFine - omegaMinusFine, omegaStepFloor);

    const scalar nutKPlusCoarse =
      evaluateLocalNut(kPlusCoarse, omegaCell, sqrtS2Cell, yWallCell, nuCell);
    const scalar nutKMinusCoarse =
      evaluateLocalNut(kMinusCoarse, omegaCell, sqrtS2Cell, yWallCell, nuCell);
    const scalar nutKPlusFine =
      evaluateLocalNut(kPlusFine, omegaCell, sqrtS2Cell, yWallCell, nuCell);
    const scalar nutKMinusFine =
      evaluateLocalNut(kMinusFine, omegaCell, sqrtS2Cell, yWallCell, nuCell);

    const scalar nutOmegaPlusCoarse =
      evaluateLocalNut(kCell, omegaPlusCoarse, sqrtS2Cell, yWallCell, nuCell);
    const scalar nutOmegaMinusCoarse =
      evaluateLocalNut(kCell, omegaMinusCoarse, sqrtS2Cell, yWallCell, nuCell);
    const scalar nutOmegaPlusFine =
      evaluateLocalNut(kCell, omegaPlusFine, sqrtS2Cell, yWallCell, nuCell);
    const scalar nutOmegaMinusFine =
      evaluateLocalNut(kCell, omegaMinusFine, sqrtS2Cell, yWallCell, nuCell);

    const scalar dNutDkCoarse = (nutKPlusCoarse - nutKMinusCoarse)/kSpanCoarse;
    const scalar dNutDkFine = (nutKPlusFine - nutKMinusFine)/kSpanFine;
    const scalar dNutDOmegaCoarse =
      (nutOmegaPlusCoarse - nutOmegaMinusCoarse)/omegaSpanCoarse;
    const scalar dNutDOmegaFine =
      (nutOmegaPlusFine - nutOmegaMinusFine)/omegaSpanFine;

    const scalar kRichardsonRatio =
      Foam::mag(dNutDkFine - dNutDkCoarse)
      / Foam::max(Foam::max(Foam::mag(dNutDkFine), Foam::mag(dNutDkCoarse)), SMALL);
    const scalar omegaRichardsonRatio =
      Foam::mag(dNutDOmegaFine - dNutDOmegaCoarse)
      / Foam::max
        (
          Foam::max(Foam::mag(dNutDOmegaFine), Foam::mag(dNutDOmegaCoarse)),
          SMALL
        );
    const scalar kStoredMismatch =
      Foam::mag(dNutDkInternal[celli] - dNutDkFine)
      / Foam::max
        (
          Foam::max(Foam::mag(dNutDkInternal[celli]), Foam::mag(dNutDkFine)),
          SMALL
        );
    const scalar omegaStoredMismatch =
      Foam::mag(dNutDOmegaInternal[celli] - dNutDOmegaFine)
      / Foam::max
        (
          Foam::max(Foam::mag(dNutDOmegaInternal[celli]), Foam::mag(dNutDOmegaFine)),
          SMALL
        );

    const scalar kRatio = Foam::max(kRichardsonRatio, kStoredMismatch);
    const scalar omegaRatio = Foam::max(omegaRichardsonRatio, omegaStoredMismatch);
    const scalar storedMismatch = Foam::max(kStoredMismatch, omegaStoredMismatch);
    const scalar sampleRatio = Foam::max(kRatio, omegaRatio);

    if (!(sampleRatio == sampleRatio) || sampleRatio >= GREAT)
    {
      continue;
    }

    ratioSumLocal += sampleRatio;
    ratioMaxLocal = Foam::max(ratioMaxLocal, sampleRatio);
    ratioKMaxLocal = Foam::max(ratioKMaxLocal, kRatio);
    ratioOmegaMaxLocal = Foam::max(ratioOmegaMaxLocal, omegaRatio);
    storedMismatchMaxLocal = Foam::max(storedMismatchMaxLocal, storedMismatch);
    if (sampleRatio > optClosureDerivativeFDStageCheckWarnRatio)
    {
      ratioWarnCountLocal++;
    }
    if (sampleRatio > optClosureDerivativeFDStageCheckFailRatio)
    {
      ratioFailCountLocal++;
    }
    ratioSamplesLocal++;
  }

  reduce(ratioSumLocal, sumOp<scalar>());
  reduce(ratioMaxLocal, maxOp<scalar>());
  reduce(ratioKMaxLocal, maxOp<scalar>());
  reduce(ratioOmegaMaxLocal, maxOp<scalar>());
  reduce(storedMismatchMaxLocal, maxOp<scalar>());
  reduce(ratioWarnCountLocal, sumOp<label>());
  reduce(ratioFailCountLocal, sumOp<label>());
  reduce(ratioSamplesLocal, sumOp<label>());

  closureDerivativeFDStageCheckSamplesUsed = ratioSamplesLocal;
  if (closureDerivativeFDStageCheckSamplesUsed > 0)
  {
    closureDerivativeFDStageCheckActive = true;
    closureDerivativeFDStageCheckMeanRatio =
      ratioSumLocal/static_cast<scalar>(closureDerivativeFDStageCheckSamplesUsed);
    closureDerivativeFDStageCheckMaxRatio = ratioMaxLocal;
    closureDerivativeFDStageCheckMaxRatioK = ratioKMaxLocal;
    closureDerivativeFDStageCheckMaxRatioOmega = ratioOmegaMaxLocal;
    closureDerivativeFDStageCheckMaxStoredMismatch = storedMismatchMaxLocal;
    const scalar ratioSamplesDen =
      Foam::max(static_cast<scalar>(closureDerivativeFDStageCheckSamplesUsed), scalar(1.0));
    closureDerivativeFDStageCheckWarnFractionObserved =
      static_cast<scalar>(ratioWarnCountLocal)/ratioSamplesDen;
    closureDerivativeFDStageCheckFailFractionObserved =
      static_cast<scalar>(ratioFailCountLocal)/ratioSamplesDen;
  }
}

const scalar closureDerivativeFDStageCheckRatio =
  closureDerivativeFDStageCheckActive
  ? Foam::max(closureDerivativeFDStageCheckMeanRatio, closureDerivativeFDStageCheckMaxRatio)
  : 0.0;
const bool closureDerivativeFDStageCheckWarnTriggered =
  closureDerivativeFDStageCheckActive
  && closureDerivativeFDStageCheckRatio > optClosureDerivativeFDStageCheckWarnRatio
  && closureDerivativeFDStageCheckWarnFractionObserved
     >= Foam::max(Foam::min(optClosureDerivativeFDStageCheckWarnFraction, scalar(1.0)), scalar(0.0));
const bool closureDerivativeFDStageCheckFailTriggered =
  closureDerivativeFDStageCheckActive
  && closureDerivativeFDStageCheckRatio > optClosureDerivativeFDStageCheckFailRatio
  && closureDerivativeFDStageCheckFailFractionObserved
     >= Foam::max(Foam::min(optClosureDerivativeFDStageCheckFailFraction, scalar(1.0)), scalar(0.0));
const label candidateClosureDerivativeFDStageCheckFailStreak =
  closureDerivativeFDStageCheckFailTriggered
  ? closureDerivativeFDStageCheckFailStreak + 1
  : 0;
const bool closureDerivativeFDStageCheckStopTriggered =
  static_cast<bool>(optClosureDerivativeFDStageCheckStopOnFailure)
  && candidateClosureDerivativeFDStageCheckFailStreak
    >= optClosureDerivativeFDStageCheckFailWindow;
const bool closureDerivativeFDStageCheckClipLimited =
    customTurbulenceActive
 && optimizationCycleTurbulentOptimizationMode
 && static_cast<bool>(optClosureDerivativeFDStageCheckEnable)
 && optClosureDerivativeFDStageCheckSamples > 0
 && closureDerivativeFDStageCheckSamplesUsed == 0;

const scalar sensitivityClipPenalty =
  Foam::min
  (
    scalar(1.0),
    maxFluidClipFraction/Foam::max(clipGuardFluidFractionThreshold, scalar(1e-3))
  );
const scalar sensitivityResidualPenalty =
  Foam::min
  (
    scalar(1.0),
    maxLinearFinalResidual/Foam::max(linearSolveWarnResidual, SMALL)
  );
const scalar sensitivityGradientPenalty =
  gradientCheckActive
  ? Foam::min(gradientConsistencyRatio/Foam::max(optGradientCheckFailRatio, SMALL), scalar(1.0))
  : scalar(0.0);
const scalar sensitivityClosurePenalty =
  closureDerivativeFDStageCheckActive
  ? Foam::min
    (
      closureDerivativeFDStageCheckRatio/Foam::max(optClosureDerivativeFDStageCheckFailRatio, SMALL),
      scalar(1.0)
    )
  : scalar(0.0);
const scalar sensitivityFinitePenalty =
  Foam::min
  (
    scalar(1.0),
    static_cast<scalar>(fsensNonFiniteCountLocal + gsensPowerNonFiniteCountLocal + gsensVolNonFiniteCountLocal)
    /Foam::max(N, scalar(1.0))
  );
const scalar sensitivityQualityScoreCurrent =
  Foam::max
  (
    scalar(0.0),
    scalar(1.0)
    - (
        0.35*sensitivityClipPenalty
      + 0.25*sensitivityResidualPenalty
      + 0.20*sensitivityClosurePenalty
      + 0.15*sensitivityGradientPenalty
      + 0.05*sensitivityFinitePenalty
      )
  );
const bool sensitivityQualityWarnCurrent =
  static_cast<bool>(optSensitivityQualityGateEnable)
  && sensitivityQualityScoreCurrent <= optSensitivityQualityWarn;
const bool sensitivityQualityFailCurrent =
  static_cast<bool>(optSensitivityQualityGateEnable)
  && sensitivityQualityScoreCurrent <= optSensitivityQualityFail;
const bool sensitivityQualityHardFailCurrent =
  sensitivityQualityFailCurrent
  && (
       maxFluidClipFraction >= clipHardRejectFluidFractionClamped
    || sensitivityFinitePenalty > 0.0
    || residualGrowthGatePassed == false
     );
const label candidateSensitivityQualityFailStreak =
  sensitivityQualityFailCurrent
  ? sensitivityQualityFailStreak + 1
  : 0;
const bool sensitivityQualityStopTriggered =
  static_cast<bool>(optSensitivityStopOnHardFail)
  && (
       sensitivityQualityHardFailCurrent
    || candidateSensitivityQualityFailStreak >= optSensitivityQualityFailWindow
     );

const bool convergenceCandidate =
  opt >= optMinIterations
  &&
  solverSeverity < 2
  &&
  convergenceConstraintViolation <= optConstraintTol
  &&
  objectiveRelativeDelta <= optObjectiveTol
  &&
  maxAbsXChangeLocal <= optStepTol
  &&
  maxAbsXhChangeLocal <= optStepTol;
const label candidateOptimizerConvergedStreak =
  convergenceCandidate ? optimizerConvergedStreak + 1 : 0;
label continuationClipBackoffStepsApplied(0);
bool rogueDirectionSwitchApplied(false);
bool rogueDirectionSwitchBudgetRearmed(false);
bool rogueDirectionObjectiveUsable(false);
bool rogueDirectionConstraintUsable(false);
scalar rogueDirectionSwitchStepAmplitudeApplied(0.0);
scalar rogueDirectionObjectiveNormApplied(0.0);
scalar rogueDirectionConstraintNormApplied(0.0);
scalar rogueDirectionSwitchClipStepScaleApplied(1.0);

if (!meritAccepted)
{
  for (i=0; i<n; i++)
  {
    x[i] = xOld[i];
    xp[i] = xpOld[i];
    xh[i] = xhOld[i];
  }
  x.correctBoundaryConditions();
  xp.correctBoundaryConditions();
  xh.correctBoundaryConditions();
  alpha=alphaMax*qu*(1-xh)/(qu+xh);
  DT=(ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoCp;
  if (optimizationCycleTurbulentOptimizationMode)
  {
    U.primitiveFieldRef() = UOldInternal;
    p.primitiveFieldRef() = pOldInternal;
    T.primitiveFieldRef() = TOldInternal;
    k.primitiveFieldRef() = kOldInternal;
    omega.primitiveFieldRef() = omegaOldInternal;
    nut.primitiveFieldRef() = nutOldInternal;
    alphat.primitiveFieldRef() = alphatOldInternal;

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();
    T.correctBoundaryConditions();
    k.correctBoundaryConditions();
    omega.correctBoundaryConditions();
    nut.correctBoundaryConditions();
    alphat.correctBoundaryConditions();

    if (customTurbulenceActive)
    {
      kClosureState.primitiveFieldRef() = kOldInternal;
      omegaClosureState.primitiveFieldRef() = omegaOldInternal;
      kClosureState.correctBoundaryConditions();
      omegaClosureState.correctBoundaryConditions();
    }
  }

  meritStepRejected = true;
  meritRejectStreak += 1;
  clipGuardStreak = candidateClipGuardStreak;
  clipGuardNoProgressStreak = candidateClipGuardNoProgressStreak;
  sensitivityQualityScore = sensitivityQualityScoreCurrent;
  sensitivityQualityWarn = sensitivityQualityWarnCurrent;
  sensitivityQualityFail = sensitivityQualityFailCurrent;
  sensitivityQualityHardFail = sensitivityQualityHardFailCurrent;
  sensitivityQualityFailStreak = candidateSensitivityQualityFailStreak;
  sensitivityQualityRequireConservativeStep = sensitivityLowConfidence || sensitivityQualityWarnCurrent;
  sensitivityQualityReuseAcceptedGradient =
      static_cast<bool>(optSensitivityReuseAcceptedGradient)
   && hasAcceptedGradientSnapshot
   && sensitivityQualityFailCurrent;
  meritTrustMoveScale =
    Foam::max(meritTrustMoveScale*rejectMoveScaleClamped, meritTrustMoveMinClamped);
  if (trustRatioDefined && trustRatio < scalar(0.0))
  {
    meritTrustMoveScale =
      Foam::max(scalar(0.8)*meritTrustMoveScale, meritTrustMoveMinClamped);
  }
  const bool deadlockTrustEscapeEligible =
       optimizationCycleTurbulentOptimizationMode
    && !optDiagnosticMode
    && meritRejectStreak >= 3*filterDeadlockRejectThreshold
    && meritTrustMoveScale <= scalar(1.05)*meritTrustMoveMinClamped
    && !clipForcedReject
    && !fieldEnvelopeForcedReject
    && sensitivityQualityGatePassed;
  if (deadlockTrustEscapeEligible)
  {
    meritTrustMoveScale =
      Foam::min(Foam::max(meritTrustMoveScale, scalar(0.10)), scalar(0.25));
  }
  const bool continuationClipBackoffEligible =
       optimizationCycleTurbulentOptimizationMode
    && clipForcedReject
    && guardConstraintViolationDelta > optConstraintWorsenTol
    && continuationIteration > 1;
  const bool continuationFieldEnvelopeBackoffEligible =
       optimizationCycleTurbulentOptimizationMode
    && fieldEnvelopeForcedReject
    && continuationIteration > 1;
  if (continuationClipBackoffEligible || continuationFieldEnvelopeBackoffEligible)
  {
    const bool severeClipLock =
         maxFluidClipFraction >= Foam::max(scalar(2.0)*clipGuardFluidFractionThreshold, scalar(0.20))
      || maxGrayClipFraction >= scalar(0.92);
    const label streakBackoffSteps =
      Foam::min(Foam::max(meritRejectStreak, label(1)), label(6));
    label targetContinuationIteration =
      Foam::max(continuationIteration - streakBackoffSteps, label(1));
    if (continuationIteration > 80)
    {
      targetContinuationIteration = Foam::min(targetContinuationIteration, label(79));
    }
    if (continuationFieldEnvelopeBackoffEligible)
    {
      targetContinuationIteration =
        Foam::max(continuationIteration - Foam::max(streakBackoffSteps, label(8)), label(1));
    }
    if
    (
         severeClipLock
      && meritRejectStreak >= clipGuardWindowDisabled
      && continuationIteration > 35
    )
    {
      targetContinuationIteration = Foam::min(targetContinuationIteration, label(35));
    }
    else if (severeClipLock && meritRejectStreak >= 3 && continuationIteration > 63)
    {
      targetContinuationIteration = Foam::min(targetContinuationIteration, label(63));
    }
    const label clipBackoffSteps =
      Foam::max(continuationIteration - targetContinuationIteration, label(0));
    if (clipBackoffSteps > 0)
    {
      continuationIteration = targetContinuationIteration;
      continuationClipBackoffStepsApplied = clipBackoffSteps;
      meritRejectStopEligibleStreak = 0;
    }
  }
  if (rogueDirectionSwitchTriggered)
  {
    const scalar rogueDirectionDerivativeAbsCap = 1e60;
    scalar objectiveDirectionAbsMaxLocal = 0.0;
    scalar constraintDirectionAbsMaxLocal = 0.0;
    for (i=0; i<n; i++)
    {
      const scalar objectiveValue = Foam::mag(static_cast<scalar>(dfdx[i]));
      const scalar constraintValue =
        Foam::mag(static_cast<scalar>(dgdx[mmaPowerConstraintIndex][i]));
      if (std::isfinite(static_cast<double>(objectiveValue)))
      {
        objectiveDirectionAbsMaxLocal =
          Foam::max(objectiveDirectionAbsMaxLocal, objectiveValue);
      }
      if (std::isfinite(static_cast<double>(constraintValue)))
      {
        constraintDirectionAbsMaxLocal =
          Foam::max(constraintDirectionAbsMaxLocal, constraintValue);
      }
    }
    reduce(objectiveDirectionAbsMaxLocal, maxOp<scalar>());
    reduce(constraintDirectionAbsMaxLocal, maxOp<scalar>());
    rogueDirectionObjectiveUsable =
         objectiveDirectionAbsMaxLocal > SMALL
      && objectiveDirectionAbsMaxLocal < rogueDirectionDerivativeAbsCap
      && std::isfinite(static_cast<double>(objectiveDirectionAbsMaxLocal));
    rogueDirectionConstraintUsable =
         constraintDirectionAbsMaxLocal > SMALL
      && constraintDirectionAbsMaxLocal < rogueDirectionDerivativeAbsCap
      && std::isfinite(static_cast<double>(constraintDirectionAbsMaxLocal));

    scalar objectiveDirectionNorm = 1.0;
    if (rogueDirectionObjectiveUsable)
    {
      scalar objectiveDirectionNormSqScaledLocal = 0.0;
      for (i=0; i<n; i++)
      {
        const scalar objectiveScaled =
          static_cast<scalar>(dfdx[i])/objectiveDirectionAbsMaxLocal;
        if (std::isfinite(static_cast<double>(objectiveScaled)))
        {
          objectiveDirectionNormSqScaledLocal += sqr(objectiveScaled);
        }
      }
      reduce(objectiveDirectionNormSqScaledLocal, sumOp<scalar>());
      objectiveDirectionNorm =
        objectiveDirectionAbsMaxLocal
       *Foam::sqrt(Foam::max(objectiveDirectionNormSqScaledLocal, SMALL));
      if (!std::isfinite(static_cast<double>(objectiveDirectionNorm)))
      {
        rogueDirectionObjectiveUsable = false;
        objectiveDirectionNorm = 1.0;
      }
    }

    scalar constraintDirectionNorm = 1.0;
    if (rogueDirectionConstraintUsable)
    {
      scalar constraintDirectionNormSqScaledLocal = 0.0;
      for (i=0; i<n; i++)
      {
        const scalar constraintScaled =
          static_cast<scalar>(dgdx[mmaPowerConstraintIndex][i])/constraintDirectionAbsMaxLocal;
        if (std::isfinite(static_cast<double>(constraintScaled)))
        {
          constraintDirectionNormSqScaledLocal += sqr(constraintScaled);
        }
      }
      reduce(constraintDirectionNormSqScaledLocal, sumOp<scalar>());
      constraintDirectionNorm =
        constraintDirectionAbsMaxLocal
       *Foam::sqrt(Foam::max(constraintDirectionNormSqScaledLocal, SMALL));
      if (!std::isfinite(static_cast<double>(constraintDirectionNorm)))
      {
        rogueDirectionConstraintUsable = false;
        constraintDirectionNorm = 1.0;
      }
    }

    const bool rogueDirectionHasBasis =
      rogueDirectionObjectiveUsable || rogueDirectionConstraintUsable;
    const scalar directionSwitchMoveFloor = Foam::max(optimizerMoveLimitActive, scalar(0.05));
    scalar rogueDirectionSwitchClipScale = 1.0;
    const bool rogueDirectionSevereClipDampingActive =
         maxFluidClipFraction >= rogueDirectionSwitchSevereClipFluidFractionClamped
      || maxGrayClipFraction >= rogueDirectionSwitchSevereClipGrayFractionClamped;
    const bool rogueDirectionClipDampingActive =
         maxFluidClipFraction >= rogueDirectionSwitchClipDampFluidFractionClamped
      || maxGrayClipFraction >= rogueDirectionSwitchClipDampGrayFractionClamped;
    if (rogueDirectionSevereClipDampingActive)
    {
      rogueDirectionSwitchClipScale =
        Foam::min(rogueDirectionSwitchClipScale, rogueDirectionSwitchSevereClipScaleClamped);
    }
    else if (rogueDirectionClipDampingActive)
    {
      rogueDirectionSwitchClipScale =
        Foam::min(rogueDirectionSwitchClipScale, rogueDirectionSwitchClipDampScaleClamped);
    }
    const scalar directionSwitchStepAmplitude =
      rogueDirectionSwitchClipScale
     *rogueDirectionSwitchStepScaleClamped
     *(rogueDirectionHasBasis ? directionSwitchMoveFloor : scalar(0.25)*directionSwitchMoveFloor);
    rogueDirectionSwitchStepAmplitudeApplied = directionSwitchStepAmplitude;
    rogueDirectionSwitchClipStepScaleApplied = rogueDirectionSwitchClipScale;
    rogueDirectionObjectiveNormApplied =
      rogueDirectionObjectiveUsable ? objectiveDirectionNorm : scalar(0.0);
    rogueDirectionConstraintNormApplied =
      rogueDirectionConstraintUsable ? constraintDirectionNorm : scalar(0.0);

    for (i=0; i<n; i++)
    {
      const scalar objectiveDirection =
        rogueDirectionObjectiveUsable
        ? -static_cast<scalar>(dfdx[i])/objectiveDirectionNorm
        : scalar(0.0);
      const scalar constraintDirection =
        rogueDirectionConstraintUsable
        ? -static_cast<scalar>(dgdx[mmaPowerConstraintIndex][i])/constraintDirectionNorm
        : scalar(0.0);
      const scalar blendedDirection =
        Foam::max
        (
          Foam::min
          (
              (scalar(1.0) - rogueDirectionSwitchConstraintBlendClamped)*objectiveDirection
            + rogueDirectionSwitchConstraintBlendClamped*constraintDirection,
            scalar(1.0)
          ),
          scalar(-1.0)
        );
      const scalar signPattern =
        ((i + opt + Pstream::myProcNo()) % 2 == 0) ? scalar(1.0) : scalar(-1.0);
      const scalar switchedValue =
        xhOld[i]
        + directionSwitchStepAmplitude
         *(
            blendedDirection
          + rogueDirectionSwitchNoiseScaleClamped*signPattern
          );
      const scalar boundedSwitchedValue =
        Foam::max(Foam::min(switchedValue, scalar(1.0)), scalar(0.0));
      xh[i] = boundedSwitchedValue;
      x[i] = boundedSwitchedValue;
      xp[i] = boundedSwitchedValue;
    }

    if (solid_area)
    {
      forAll(cells_solid, cellI)
      {
        const label c = cells_solid[cellI];
        xh[c] = 0.0;
        x[c] = 0.0;
        xp[c] = 0.0;
      }
    }
    if (fluid_area)
    {
      forAll(cells_fluid, cellI)
      {
        const label c = cells_fluid[cellI];
        xh[c] = 1.0;
        x[c] = 1.0;
        xp[c] = 1.0;
      }
    }

    x.correctBoundaryConditions();
    xp.correctBoundaryConditions();
    xh.correctBoundaryConditions();
    alpha=alphaMax*qu*(1-xh)/(qu+xh);
    DT=(ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoCp;

    rogueDirectionSwitchApplied = true;
    rogueDirectionSwitchCount += 1;
    rogueDirectionLastSwitchIteration = opt;
    rogueDirectionWatchActive = false;
    rogueDirectionWatchGraceRemaining = 0;
    rogueDirectionWatchStartIteration = -1;
    stopOptimization = false;
    optimizerStopReason = word("rogueDirectionSwitch");
    meritRejectStopEligibleStreak = 0;
  }
  const bool meritScheduleRebaseEligible =
       hasPreviousAcceptedMerit
    && (
         (
              scheduleIterationCurrent > acceptedMeritScheduleIteration
           && !constraintWorsenScheduleReady
         )
         || constraintWorsenScheduleReadyTransition
       );
  if (meritScheduleRebaseEligible)
  {
    previousAcceptedMerit = meritValue;
    acceptedMeritScheduleIteration = scheduleIterationCurrent;
    previousGuardConstraintViolation = guardConstraintViolation;
    hasPreviousGuardConstraintViolation = true;
    previousMaxConstraintViolation = maxConstraintViolation;
    hasPreviousMaxConstraintViolation = true;
    meritScheduleRebaseCount += 1;
  }
  if (rogueDirectionSwitchApplied)
  {
    stopOptimization = false;
    optimizerStopReason = word("rogueDirectionSwitch");
  }
  else if
  (
       !optDiagnosticMode
    && static_cast<bool>(optProductionContinueAfterSafetyReject)
    && !fieldEnvelopeForcedReject
    && !clipHardForcedReject
    && !sensitivityQualityHardFailCurrent
    && !sensitivityQualityStopTriggered
  )
  {
    stopOptimization = false;
    optimizerStopReason = word("safetyRejectRecovered");
    meritRecoveryRequest = true;
  }
  else
  {
    if (fieldEnvelopeForcedReject)
    {
      stopOptimization = true;
      optimizerStopReason = word("fieldEnvelopeReject");
    }
    else
    {
      if (clipHardForcedReject)
      {
        stopOptimization = true;
        optimizerStopReason = word("clipHardReject");
      }
      else if (clipGuardTriggered && !clipGuardRejectMode)
      {
        stopOptimization = true;
        optimizerStopReason = word("turbulenceClipLock");
      }
      else
      {
        stopOptimization = true;
        optimizerStopReason = word("safetyReject");
      }
    }
  }
  if
  (
      optDiagnosticMode
   && (
        !healthGatesPassed
     || sensitivityQualityHardFailCurrent
     || sensitivityQualityStopTriggered
      )
  )
  {
    stopOptimization = true;
    optimizerStopReason = word("diagnosticHealthFail");
  }
}
else
{
  meritStepRejected = false;
  meritRejectStreak = 0;
  meritRejectStopEligibleStreak = 0;
  meritTrustMoveScale = 1.0;
  const bool rogueDirectionSwitchRearmEligible =
       optimizationCycleTurbulentOptimizationMode
    && static_cast<bool>(optRogueDirectionSwitchEnable)
    && !optDiagnosticMode
    && rogueDirectionSwitchCount > 0
    && !clipForcedReject
    && guardConstraintViolationDelta
       < -Foam::max(optConstraintWorsenTol, clipGuardNoProgressTolReference);
  if (rogueDirectionSwitchRearmEligible)
  {
    rogueDirectionSwitchCount =
      Foam::max(rogueDirectionSwitchCount - 1, label(0));
    rogueDirectionSwitchBudgetRearmed = true;
    rogueDirectionWatchActive = false;
    rogueDirectionWatchGraceRemaining = 0;
    rogueDirectionWatchStartIteration = -1;
  }

  previousObjectiveValue = objectiveValue;
  hasPreviousObjectiveValue = true;
  previousPowerConstraintMargin = powerConstraintMargin;
  hasPreviousPowerConstraintMargin = true;
  previousPowerConstraintFinalMargin = powerConstraintFinalMargin;
  hasPreviousPowerConstraintFinalMargin = true;
  previousVolumeConstraintMargin = volumeConstraintMargin;
  hasPreviousVolumeConstraintMargin = true;
  previousThermalTmaxConstraintMargin = thermalTmaxConstraintMargin;
  hasPreviousThermalTmaxConstraintMargin = static_cast<bool>(optThermalEnvelopeConstraintEnable);
  previousThermalTminConstraintMargin = thermalTminConstraintMargin;
  hasPreviousThermalTminConstraintMargin = static_cast<bool>(optThermalEnvelopeConstraintEnable);
  previousMaxConstraintViolation = maxConstraintViolation;
  hasPreviousMaxConstraintViolation = true;
  previousGuardConstraintViolation = guardConstraintViolation;
  hasPreviousGuardConstraintViolation = true;
  previousAcceptedMerit = meritValue;
  hasPreviousAcceptedMerit = true;
  acceptedMeritScheduleIteration = scheduleIterationCurrent;
  continuationDeadlockRecoveryCount = 0;
  if (filterAcceptanceEnabled)
  {
    const scalar constraintCandidate = guardConstraintViolation;
    const scalar objectiveCandidate = objectiveValue;
    std::vector<scalar> filterConstraintPruned;
    std::vector<scalar> filterObjectivePruned;
    filterConstraintPruned.reserve(filterConstraintArchive.size() + 1);
    filterObjectivePruned.reserve(filterObjectiveArchive.size() + 1);
    for (std::size_t fi = 0; fi < filterConstraintArchive.size(); ++fi)
    {
      const bool dominatedByCandidate =
           filterConstraintArchive[fi] >= constraintCandidate - scalar(0.0)
        && filterObjectiveArchive[fi] >= objectiveCandidate - scalar(0.0);
      if (!dominatedByCandidate)
      {
        filterConstraintPruned.push_back(filterConstraintArchive[fi]);
        filterObjectivePruned.push_back(filterObjectiveArchive[fi]);
      }
    }
    filterConstraintPruned.push_back(constraintCandidate);
    filterObjectivePruned.push_back(objectiveCandidate);
    while (static_cast<label>(filterConstraintPruned.size()) > label(0))
    {
      filterConstraintPruned.erase(filterConstraintPruned.begin());
      filterObjectivePruned.erase(filterObjectivePruned.begin());
    }
    filterConstraintArchive.swap(filterConstraintPruned);
    filterObjectiveArchive.swap(filterObjectivePruned);
  }
  powerConstraintWorsenStreak = candidatePowerConstraintWorsenStreak;
  maxViolationWorsenStreak = candidateMaxViolationWorsenStreak;
  objectiveWorsenStreak = candidateObjectiveWorsenStreak;
  clipGuardStreak = candidateClipGuardStreak;
  clipGuardNoProgressStreak = candidateClipGuardNoProgressStreak;
  gradientCheckFailStreak = candidateGradientCheckFailStreak;
  closureDerivativeFDStageCheckFailStreak =
    candidateClosureDerivativeFDStageCheckFailStreak;
  sensitivityQualityScore = sensitivityQualityScoreCurrent;
  sensitivityQualityWarn = sensitivityQualityWarnCurrent;
  sensitivityQualityFail = sensitivityQualityFailCurrent;
  sensitivityQualityHardFail = sensitivityQualityHardFailCurrent;
  sensitivityQualityFailStreak = candidateSensitivityQualityFailStreak;
  sensitivityQualityRequireConservativeStep = sensitivityLowConfidence || sensitivityQualityWarnCurrent;
  sensitivityQualityReuseAcceptedGradient = false;
  optimizerConvergedStreak = candidateOptimizerConvergedStreak;
  acceptedDfdx = dfdx;
  acceptedDgdx = dgdx;
  hasAcceptedGradientSnapshot = true;
  hasAcceptedStateSnapshot = true;

  if (constraintScheduleHorizonGuardTriggered)
  {
    stopOptimization = true;
    optimizerStopReason = word("constraintScheduleHorizon");
  }
  else if (constraintWorsenGuardTriggered)
  {
    stopOptimization = true;
    optimizerStopReason = word("constraintWorsening");
  }
  else if (constraintWorsenNetRiseTriggered)
  {
    stopOptimization = true;
    optimizerStopReason = word("constraintNetRise");
  }
  else if (constraintWorsenReboundStopTriggered)
  {
    stopOptimization = true;
    optimizerStopReason = word("constraintReboundRise");
  }
  else if (objectiveWorsenGuardTriggered)
  {
    stopOptimization = true;
    optimizerStopReason = word("objectiveBacktrack");
  }
  else if (clipGuardTriggered && !clipGuardRejectMode)
  {
    stopOptimization = true;
    optimizerStopReason = word("turbulenceClipLock");
  }
  else if (gradientCheckStopTriggered)
  {
    stopOptimization = true;
    optimizerStopReason = word("gradientConsistency");
  }
  else if (closureDerivativeFDStageCheckStopTriggered)
  {
    stopOptimization = true;
    optimizerStopReason = word("closureDerivativeConsistency");
  }
  else if (sensitivityQualityStopTriggered)
  {
    stopOptimization = true;
    optimizerStopReason = word("sensitivityQuality");
  }
  else
  {
    stopOptimization = optimizerConvergedStreak >= optStallWindow;
    optimizerStopReason = stopOptimization ? word("converged") : word("running");
  }
}
const bool globalizationRejectStopEligible =
  !static_cast<bool>(optRejectStopRequireConstraintWindow)
  || constraintWorsenCountingActive
  || (
       scheduleIterationCurrent >= optConstraintWorsenMinIterations
    && (
         constraintWorsenScheduleReady
      || constraintScheduleHorizonGuardTriggered
       )
     );
const bool rejectIndicatesDeterioration =
     guardConstraintViolationDelta > 0
  || objectiveWorsening
  || meritDeteriorating;
const bool continuationClipBackoffRecoveryActive =
  continuationClipBackoffStepsApplied > 0;
const bool meritRejectCountsTowardStop =
  meritStepRejected
  && !rogueDirectionSwitchApplied
  && !continuationClipBackoffRecoveryActive
  && (
       !clipForcedReject
    || rejectIndicatesDeterioration
     );
if (meritRejectCountsTowardStop && globalizationRejectStopEligible)
{
  meritRejectStopEligibleStreak += 1;
}
else if (!meritRejectCountsTowardStop || !globalizationRejectStopEligible)
{
  meritRejectStopEligibleStreak = 0;
}
if
(
    static_cast<bool>(false)
 && label(0) > 0
 && meritRejectStopEligibleStreak >= label(0)
)
{
  if (globalizationRejectStopEligible)
  {
    stopOptimization = true;
    optimizerStopReason = word("globalizationRejectStreak");
  }
  else
  {
    stopOptimization = false;
    optimizerStopReason = word("globalizationReject");
  }
}
const bool continuationFrozenDuringIteration =
  continuationFrozenByReject || continuationFrozenByLowTrust || continuationFrozenByClipHold;
const bool scheduleHorizonNoSlackUnderFreeze =
  !constraintWorsenScheduleReady
  && continuationFrozenDuringIteration
  && powerConstraintScheduleIterationsToSettle > 0
  && optimizationIterationsRemaining > 0
  && powerConstraintScheduleIterationsToSettle >= optimizationIterationsRemaining;
const bool scheduleHorizonCatchUpNeeded =
  static_cast<bool>(optConstraintScheduleHorizonGuard)
  && (
       (
            !constraintWorsenScheduleReady
         && !powerConstraintScheduleWillSettleInRun
         && scheduleHorizonLagIterations > 0
       )
    || scheduleHorizonNoSlackUnderFreeze
     )
  && scheduleHorizonConstraintViolation > optConstraintTol;
const bool scheduleHorizonDeadlockOnReject =
  meritStepRejected && scheduleHorizonCatchUpNeeded;
const bool deadlockRecoveryBlockedByGlobalization =
  meritStepRejected
  && (
       fieldEnvelopeForcedReject
    || clipHardForcedReject
    || (clipGuardTriggered && !clipGuardRejectMode)
    || continuationClipBackoffRecoveryActive
    || clipForcedReject
    || objectiveWorsening
    || guardConstraintViolationDelta > optConstraintWorsenTol
     );
label deadlockCatchUpRequest(0);
if (scheduleHorizonCatchUpNeeded)
{
  const bool continuationDeadlockRecoveryEligible =
    scheduleHorizonDeadlockOnReject
    && meritTrustMoveScale <= scalar(1.05)*meritTrustMoveMinClamped
    && static_cast<bool>(optContinuationDeadlockRecoveryEnable)
    && !optDiagnosticMode
    && continuationDeadlockRecoveryCount < optContinuationDeadlockRecoveryMaxOverrides
    && !deadlockRecoveryBlockedByGlobalization;
  if (continuationDeadlockRecoveryEligible)
  {
    const label deadlockRecoveryLagBufferClamped =
      Foam::max(optContinuationDeadlockRecoveryLagBuffer, label(0));
    deadlockCatchUpRequest =
      scheduleHorizonLagIterations + deadlockRecoveryLagBufferClamped;
    const label deadlockCatchUpSteps =
      Foam::min
      (
        deadlockCatchUpRequest,
        Foam::max(optContinuationDeadlockRecoveryMaxCatchUpSteps, label(0))
      );
    continuationAdvanceOverrideOnReject = true;
    continuationAdvanceCatchUpSteps =
      Foam::max(continuationAdvanceCatchUpSteps, deadlockCatchUpSteps);
    continuationDeadlockRecoveryCount += 1;
    stopOptimization = false;
    if (rogueDirectionSwitchApplied)
    {
      optimizerStopReason = word("rogueDirectionSwitch");
    }
    else
    {
      optimizerStopReason = meritStepRejected ? word("globalizationReject") : word("running");
    }
  }
  else if (scheduleHorizonDeadlockOnReject && !deadlockRecoveryBlockedByGlobalization)
  {
    stopOptimization = true;
    optimizerStopReason = word("constraintScheduleHorizonReject");
  }
}
const bool diagnosticHealthEvent =
     fieldEnvelopeForcedReject
  || clipHardForcedReject
  || sensitivityQualityHardFailCurrent
  || sensitivityQualityStopTriggered
  || !healthGatesPassed;
if (optDiagnosticMode && diagnosticHealthEvent)
{
  stopOptimization = true;
  optimizerStopReason = word("diagnosticHealthFail");
}

const label reportedPowerConstraintWorsenStreak =
  meritStepRejected ? candidatePowerConstraintWorsenStreak : powerConstraintWorsenStreak;
const label reportedMaxViolationWorsenStreak =
  meritStepRejected ? candidateMaxViolationWorsenStreak : maxViolationWorsenStreak;
const label reportedObjectiveWorsenStreak =
  meritStepRejected ? candidateObjectiveWorsenStreak : objectiveWorsenStreak;
const label reportedClipGuardStreak =
  meritStepRejected ? candidateClipGuardStreak : clipGuardStreak;
const label reportedClipGuardNoProgressStreak =
  meritStepRejected ? candidateClipGuardNoProgressStreak : clipGuardNoProgressStreak;
const label reportedOptimizerConvergedStreak =
  meritStepRejected ? candidateOptimizerConvergedStreak : optimizerConvergedStreak;
const label reportedGradientCheckFailStreak =
  meritStepRejected ? candidateGradientCheckFailStreak : gradientCheckFailStreak;
const label reportedClosureDerivativeFDStageCheckFailStreak =
  meritStepRejected
  ? candidateClosureDerivativeFDStageCheckFailStreak
  : closureDerivativeFDStageCheckFailStreak;
const label reportedContinuationDeadlockRecoveryCount =
  continuationDeadlockRecoveryCount;
const label reportedContinuationAdvanceCatchUpSteps =
  continuationAdvanceCatchUpSteps;
const label reportedDeadlockCatchUpRequest =
  deadlockCatchUpRequest;
const bool reportedContinuationFrozenByLowTrust =
  continuationFrozenByLowTrust;
const bool reportedContinuationFrozenAny =
  continuationFrozenAny;
const label rogueDirectionSwitchCooldownRemaining =
  Foam::max
  (
    rogueDirectionSwitchCooldownClamped - (opt - rogueDirectionLastSwitchIteration),
    label(0)
  );

previousMaxFluidClipFraction = maxFluidClipFraction;
previousMaxClipFraction = maxClipFraction;
previousGuardConstraintViolationDelta = guardConstraintViolationDelta;
hasPreviousClipDiagnostics = true;
previousHealthGateFailed = !healthGatesPassed;
previousThermalEnvelopeViolation = thermalEnvelopeViolation;
previousMaxLinearFinalResidual = maxLinearFinalResidual;
hasPreviousMaxLinearFinalResidual = true;
const bool firstFailureHealthEvent =
     optimizerStopReason != word("converged")
  && (
       fieldEnvelopeForcedReject
    || clipHardForcedReject
    || sensitivityQualityHardFailCurrent
    || sensitivityQualityStopTriggered
    || !healthGatesPassed
     );
const bool triggerFirstFailureSnapshot =
     static_cast<bool>(optFirstFailureSnapshotEnable)
  && !firstFailureSnapshotWritten
  && firstFailureHealthEvent;
if (triggerFirstFailureSnapshot)
{
  x.write();
  xp.write();
  xh.write();
  U.write();
  p.write();
  T.write();
  k.write();
  omega.write();
  fsensMeanT.write();
  gsensPowerDiss.write();
  gsensVol.write();
  firstFailureSnapshotWritten = true;
}

if (Pstream::master())
{
  const scalar clockTime = runTime.elapsedClockTime();
  const char* designHealth =
    (xNonFiniteCountLocal > 0 || xLocalMin < -1e-8 || xLocalMax > 1.0 + 1e-8 || xhLocalMin < -1e-8 || xhLocalMax > 1.0 + 1e-8)
    ? "WARN"
    : (grayVolumeFraction > 0.50 ? "OBS" : "OK");
  const char* constraintHealth =
    (powerConstraintMargin > 0 || volumeConstraintMargin > 0 || thermalEnvelopeViolation > 0)
    ? "WARN"
    : "OK";
  const char* sensitivityHealth =
    (fsensNonFiniteCountLocal > 0 || gsensPowerNonFiniteCountLocal > 0 || gsensVolNonFiniteCountLocal > 0)
    ? "WARN"
    : "OK";
  const scalar severeTransitionFluidThreshold =
    Foam::max(scalar(0.08), scalar(0.8)*clipGuardFluidFractionThreshold);
  const bool severeTransitionClip =
    customTurbulenceActive
    && maxGrayClipFraction >= 0.85
    && maxFluidClipFraction >= severeTransitionFluidThreshold
    && maxConstraintViolation > optConstraintTol;
  const label turbulenceClipSeverity =
    !customTurbulenceActive ? -1
    : ((maxFluidClipFraction >= turbulenceClipWarnFraction || severeTransitionClip) ? 2
      : (maxFluidClipFraction >= turbulenceClipObsFraction ? 1
        : (maxGrayClipFraction >= turbulenceClipWarnFraction ? 1 : 0)));
  const char* turbulenceClipHealth =
    turbulenceClipSeverity == 2 ? "WARN"
    : (turbulenceClipSeverity == 1 ? "OBS"
      : (turbulenceClipSeverity == 0 ? "OK" : "SKIP"));
  const bool stopSuccessful =
    stopOptimization && optimizerStopReason == word("converged");
  const bool stopFailure =
    stopOptimization && !stopSuccessful;
  const char* stopOutcome =
    !stopOptimization ? "running" : (stopSuccessful ? "converged" : "guard-stop");
  const char* solverHealth =
    (solverSeverity >= 2 || turbulenceClipSeverity >= 2)
    ? "WARN"
    : ((solverSeverity == 1 || turbulenceClipSeverity == 1) ? "OBS" : "OK");
  if (triggerFirstFailureSnapshot)
  {
    std::ofstream firstFailureSnapshotStream(firstFailureSnapshotFile, std::ios::app);
    firstFailureSnapshotStream
      << "Iter=" << opt
      << " reason=" << optimizerStopReason
      << " objective=" << objectiveValue
      << " maxConstraint=" << maxConstraintViolation
      << " guardConstraint=" << guardConstraintViolation
      << " maxLinearFinalResidual=" << maxLinearFinalResidual
      << " residualGrowthRatio=" << linearResidualGrowthRatio
      << " maxFluidClip=" << maxFluidClipFraction
      << " maxGrayClip=" << maxGrayClipFraction
      << " sensQuality=" << sensitivityQualityScoreCurrent
      << " sensFail=" << (sensitivityQualityFailCurrent ? 1 : 0)
      << " sensHardFail=" << (sensitivityQualityHardFailCurrent ? 1 : 0)
      << " thermalTmaxKS=" << thermalTmaxConstraintMargin
      << " thermalTminKS=" << thermalTminConstraintMargin
      << '\n';
  }

  std::ofstream debugLogFile(debugOptimizerLogFile, std::ios::app);
  std::ofstream optimizerDirectionEventLogFile;
  if (static_cast<bool>(optRogueDirectionSwitchEnable))
  {
    optimizerDirectionEventLogFile.open(optimizerDirectionEventFile, std::ios::app);
  }
  const int itemWidth = 18;
  const int valueWidth = 34;
  const int tableWidth = itemWidth + 4*valueWidth;
  auto fmtScalar =
  [&](const scalar value)
  {
    std::ostringstream os;
    os << std::scientific << std::setprecision(6) << value;
    return os.str();
  };
  auto fmtCount =
  [&](const label value)
  {
    std::ostringstream os;
    os << value;
    return os.str();
  };
  auto writeSectionHeader =
  [&](const char* title)
  {
    debugLogFile
      << "\n[" << title << "]\n"
      << "  "
      << std::left
      << std::setw(itemWidth) << "Item"
      << std::setw(valueWidth) << "Value 1"
      << std::setw(valueWidth) << "Value 2"
      << std::setw(valueWidth) << "Value 3"
      << std::setw(valueWidth) << "Value 4"
      << "\n"
      << "  "
      << std::string(tableWidth, '-')
      << "\n";
  };
  auto writeRow =
  [&](const char* item, const std::string& value1, const std::string& value2 = "", const std::string& value3 = "", const std::string& value4 = "")
  {
    debugLogFile
      << "  "
      << std::left
      << std::setw(itemWidth) << item
      << std::setw(valueWidth) << value1
      << std::setw(valueWidth) << value2
      << std::setw(valueWidth) << value3
      << std::setw(valueWidth) << value4
      << "\n";
  };

  if (rogueDirectionWatchStarted && optimizerDirectionEventLogFile.is_open())
  {
    optimizerDirectionEventLogFile
      << "Iter=" << opt
      << " Event=watch_start"
      << " signal=rogue"
      << " rejectSignal=" << (rogueRejectSignal ? 1 : 0)
      << " clipSignal=" << (rogueClipSignal ? 1 : 0)
      << " graceRemaining=" << rogueDirectionWatchGraceRemaining
      << " rejectStreak=" << candidateMeritRejectStreak
      << " fluidClip=" << maxFluidClipFraction
      << " grayClip=" << maxGrayClipFraction
      << '\n';
  }
  if (rogueDirectionWatchGraceAdvanced && optimizerDirectionEventLogFile.is_open())
  {
    optimizerDirectionEventLogFile
      << "Iter=" << opt
      << " Event=watch_grace_countdown"
      << " graceRemaining=" << rogueDirectionWatchGraceRemaining
      << " watchStart=" << rogueDirectionWatchStartIterationReport
      << " rejectStreak=" << candidateMeritRejectStreak
      << '\n';
  }
  if
  (
      rogueDirectionWatchCleared
   && !rogueDirectionSwitchApplied
   && optimizerDirectionEventLogFile.is_open()
  )
  {
    optimizerDirectionEventLogFile
      << "Iter=" << opt
      << " Event=watch_cleared"
      << " reason=signal_recovered"
      << " rejectSignal=" << (rogueRejectSignal ? 1 : 0)
      << " clipSignal=" << (rogueClipSignal ? 1 : 0)
      << '\n';
  }
  if (rogueDirectionSwitchApplied)
  {
    if (optimizerDirectionEventLogFile.is_open())
    {
      optimizerDirectionEventLogFile
        << "Iter=" << opt
        << " Event=switch_applied"
        << " reason=rogue_verified"
        << " switchCount=" << rogueDirectionSwitchCount
        << " maxSwitchCount=" << rogueDirectionSwitchMaxCountClamped
        << " cooldown=" << rogueDirectionSwitchCooldownClamped
        << " stepAmplitude=" << rogueDirectionSwitchStepAmplitudeApplied
        << " clipStepScale=" << rogueDirectionSwitchClipStepScaleApplied
        << " objectiveUsable=" << (rogueDirectionObjectiveUsable ? 1 : 0)
        << " constraintUsable=" << (rogueDirectionConstraintUsable ? 1 : 0)
        << " objectiveNorm=" << rogueDirectionObjectiveNormApplied
        << " constraintNorm=" << rogueDirectionConstraintNormApplied
        << " constraintBlend=" << rogueDirectionSwitchConstraintBlendClamped
        << " noiseScale=" << rogueDirectionSwitchNoiseScaleClamped
        << " rejectStreak=" << meritRejectStreak
        << " fluidClip=" << maxFluidClipFraction
        << " grayClip=" << maxGrayClipFraction
        << '\n';
    }

    Info<< "\n============================================================\n"
        << "Optimizer Direction Switch Event\n"
        << "Iteration: " << opt << '\n'
        << "Switch count: " << rogueDirectionSwitchCount
        << "/" << rogueDirectionSwitchMaxCountClamped << '\n'
        << "Step amplitude: " << rogueDirectionSwitchStepAmplitudeApplied << '\n'
        << "Clip step scale: " << rogueDirectionSwitchClipStepScaleApplied << '\n'
        << "Objective usable: " << (rogueDirectionObjectiveUsable ? "yes" : "no")
        << " | Constraint usable: " << (rogueDirectionConstraintUsable ? "yes" : "no") << '\n'
        << "Blend (constraint): " << rogueDirectionSwitchConstraintBlendClamped << '\n'
        << "Reject streak: " << meritRejectStreak
        << " | fluidClip: " << maxFluidClipFraction
        << " | grayClip: " << maxGrayClipFraction << '\n'
        << "============================================================\n" << endl;
  }
  if (rogueDirectionSwitchBudgetRearmed && optimizerDirectionEventLogFile.is_open())
  {
    optimizerDirectionEventLogFile
      << "Iter=" << opt
      << " Event=switch_budget_rearmed"
      << " reason=accepted_feasibility_progress"
      << " switchCount=" << rogueDirectionSwitchCount
      << " guardConstraintDelta=" << guardConstraintViolationDelta
      << '\n';
  }

  debugLogFile
    << "\n"
    << "====================================================================================================\n"
    << " Optimizer Debug Dashboard\n"
    << "----------------------------------------------------------------------------------------------------\n";

  writeSectionHeader("Overview");
  writeRow("Run", "Iter="+fmtCount(opt), "Clock[s]="+fmtScalar(clockTime), "ObjMode="+std::string(objectiveLabel.c_str()), "ObjVal="+fmtScalar(objectiveValue));
  writeRow("Health", std::string("Design=")+designHealth, std::string("Constr=")+constraintHealth, std::string("Sens=")+sensitivityHealth, std::string("Solver=")+solverHealth);
  writeRow
  (
    "Convergence",
    "objRel="+fmtScalar(objectiveRelativeDelta),
    "maxViol="+fmtScalar(maxConstraintViolation),
    "convViol="+fmtScalar(convergenceConstraintViolation),
    "streak="+fmtCount(reportedOptimizerConvergedStreak)+"/"+fmtCount(optStallWindow)
  );
  writeRow
  (
    "Viol trend",
    "dMax="+fmtScalar(maxConstraintViolationDelta),
    "dPower="+fmtScalar(powerConstraintMarginDelta),
    "dPowerFinal="+fmtScalar(powerConstraintFinalMarginDelta),
    "dVol="+fmtScalar(volumeConstraintMarginDelta)
  );
  writeRow
  (
    "Guard metric",
    "max="+fmtScalar(guardConstraintViolation),
    "dMax="+fmtScalar(guardConstraintViolationDelta),
    "pFinal="+fmtScalar(powerConstraintFinalMargin),
    "rise="+fmtScalar(guardConstraintViolationNetRise)
  );
  writeRow("Viol streak", "dMax+="+fmtCount(reportedMaxViolationWorsenStreak), "dPower+="+fmtCount(reportedPowerConstraintWorsenStreak));
  writeRow
  (
    "Viol guard",
    std::string("en=") + (constraintWorsenGuardEnabled ? "yes" : "no")
      + " cnt=" + (
          constraintWorsenCountingActive ? "yes"
          : (
              !constraintWorsenGuardEnabled ? "no"
              : (
                  opt < optConstraintWorsenMinIterations
                  ? "waitIter"
                  : (constraintWorsenScheduleReady ? "no" : "waitSched")
                )
            )
        ),
    std::string("phase=") + (powerConstraintScheduleSettled ? "settled" : "tightening")
      + " tol=" + fmtScalar(optConstraintWorsenTol),
    "streak="+fmtCount(reportedMaxViolationWorsenStreak)+"/"+fmtCount(optConstraintWorsenWindow),
    "schedIter="+fmtCount(scheduleIterationCurrent)
      + "/" + fmtCount(optConstraintWorsenMinIterations)
      + " den=" + fmtScalar(powerConstraintScheduleDenomCurrent)
  );
  writeRow
  (
    "Sched horizon",
    "toSettle="+fmtCount(powerConstraintScheduleIterationsToSettle),
    "itersLeft="+fmtCount(optimizationIterationsRemaining),
    "settlesInRun="+std::string(powerConstraintScheduleWillSettleInRun ? "yes" : "no"),
    "lag="+fmtCount(scheduleHorizonLagIterations)
      + " tol=" + fmtCount(constraintScheduleHorizonLagTolerance)
      + " trigger=" + std::string(constraintScheduleHorizonGuardTriggered ? "yes" : "no")
  );
  writeRow
  (
    "Guard rise",
    std::string("start=")
      + (hasGuardConstraintViolationAtCountStart ? fmtScalar(guardConstraintViolationAtCountStart) : "n/a"),
    "rise="+fmtScalar(guardConstraintViolationNetRise)+" d="+fmtScalar(guardConstraintViolationDelta),
    "tol="+fmtScalar(effectiveConstraintWorsenNetRiseTol)
      + " age=" + fmtCount(guardConstraintViolationCountAge)
      + "/" + fmtCount(effectiveConstraintWorsenNetRiseMinAge),
    "trigger="+std::string(constraintWorsenNetRiseTriggered ? "yes" : "no")
  );
  writeRow
  (
    "Rebound rise",
    std::string("min=")
      + (hasGuardConstraintViolationAtCountStart ? fmtScalar(guardConstraintViolationMinSinceCountStart) : "n/a"),
    "rise="+fmtScalar(guardConstraintViolationReboundRise)
      + " streak=" + fmtCount(reportedMaxViolationWorsenStreak),
    "win="+fmtCount(effectiveConstraintReboundWindow)
      + " age=" + fmtCount(guardConstraintViolationCountAge)
      + "/" + fmtCount(effectiveConstraintReboundMinAge),
    "tol="+fmtScalar(effectiveConstraintReboundRiseTol)
      + " trig=" + std::string(constraintWorsenReboundStopTriggered ? "yes" : "no")
  );
  writeRow
  (
    "Rebound gate",
    "objW="+std::string(objectiveWorsening ? "yes" : "no"),
    "streak="+fmtCount(reportedMaxViolationWorsenStreak)+">"+fmtCount(effectiveConstraintReboundWindow),
    "netRise="+fmtScalar(guardConstraintViolationNetRise)
      + " tol=" + fmtScalar(constraintReboundNetRiseCorroborationTol),
    "corr="+std::string(constraintWorsenReboundCorroborated ? "yes" : "no")
  );
  writeRow
  (
    "Obj guard",
    "worsen="+std::string(objectiveWorsening ? "yes" : "no")
      + " dObj=" + fmtScalar(objectiveDelta),
    "streak="+fmtCount(reportedObjectiveWorsenStreak)+"/"+fmtCount(optObjectiveWorsenWindow),
    "tol="+fmtScalar(optObjectiveWorsenTol),
    "trigger="+std::string(objectiveWorsenGuardTriggered ? "yes" : "no")
  );
  writeRow
  (
    "Clip guard",
    "fluid="+fmtScalar(maxFluidClipFraction),
    "streak="+fmtCount(reportedClipGuardStreak)+"/"+fmtCount(clipGuardWindowDisabled)
      + " noProg=" + fmtCount(reportedClipGuardNoProgressStreak)
      + "/" + fmtCount(clipGuardWindowDisabled),
    "thr="+fmtScalar(clipGuardFluidFractionThreshold)
      + " noProgTol=" + fmtScalar(clipGuardNoProgressTolClamped),
    "trigger="+std::string(clipGuardTriggered ? "yes" : "no")
      + " rejMode=" + std::string(clipGuardRejectMode ? "yes" : "no")
  );
  writeRow
  (
    "Clip hard",
    "frac="+std::string(clipHardRejectByFraction ? "yes" : "no")
      + " jump=" + std::string(clipHardRejectByJump ? "yes" : "no"),
    "fThr="+fmtScalar(clipHardRejectFluidFractionClamped)
      + " gThr=" + fmtScalar(clipHardRejectGrayFractionClamped),
    "dFluid="+fmtScalar(maxFluidClipFraction - previousMaxFluidClipFraction)
      + " dThr=" + fmtScalar(clipHardRejectFluidJumpClamped),
    "forced="+std::string(clipHardForcedReject ? "yes" : "no")
  );
  writeRow
  (
    "Recovery",
    "active="+std::string(constraintRecoveryActive ? "yes" : "no"),
    "streak="+fmtCount(reportedMaxViolationWorsenStreak)+"/"+fmtCount(optRecoveryWorsenStreak),
    "contScale="+fmtScalar(continuationRampScaleActive)+" delStep="+fmtScalar(delRampStepActive),
    "move="+fmtScalar(optimizerMoveLimitActive)
  );
  writeRow
  (
    "Acceptance",
    "stepAccepted="+std::string(meritAccepted ? "yes" : "no"),
    "stepRejected="+std::string(meritStepRejected ? "yes" : "no"),
    "clipForced="+std::string(clipForcedReject ? "yes" : "no"),
    "hardReject="+std::string(clipHardForcedReject ? "yes" : "no")
  );
  writeRow
  (
    "Envelope",
    "objAbs="+fmtScalar(objectiveAbsValue)
      + " lim=" + fmtScalar(optObjectiveAbsMax)
      + " viol=" + std::string(objectiveEnvelopeViolated ? "yes" : "no"),
    "tAbs="+fmtScalar(thermalFieldAbsMaxObserved)
      + " lim=" + fmtScalar(optThermalFieldAbsMax)
      + " viol=" + std::string(thermalFieldEnvelopeViolated ? "yes" : "no"),
    "forcedReject="+std::string(fieldEnvelopeForcedReject ? "yes" : "no"),
    "clipReject="+std::string(clipForcedReject ? "yes" : "no")
      + " hard=" + std::string(clipHardForcedReject ? "yes" : "no")
  );
  writeRow
  (
    "Step ctrl",
    "rejectStreak="+fmtCount(meritRejectStreak),
    "stepRejected="+std::string(meritStepRejected ? "yes" : "no"),
    "stopElig="+std::string(globalizationRejectStopEligible ? "yes" : "no"),
    "hCatch="+std::string(scheduleHorizonCatchUpNeeded ? "yes" : "no")
  );
  writeRow
  (
    "Sens quality",
    "score="+fmtScalar(sensitivityQualityScoreCurrent)
      + " warn=" + std::string(sensitivityQualityWarnCurrent ? "yes" : "no"),
    "fail="+std::string(sensitivityQualityFailCurrent ? "yes" : "no")
      + " hard=" + std::string(sensitivityQualityHardFailCurrent ? "yes" : "no")
      + " streak=" + fmtCount(candidateSensitivityQualityFailStreak)
      + "/" + fmtCount(optSensitivityQualityFailWindow),
    "reuse="+std::string(sensitivityQualityReuseAcceptedGradient ? "yes" : "no")
      + " lowConf=" + std::string(sensitivityLowConfidence ? "yes" : "no")
      + " gSnap=" + std::string(hasAcceptedGradientSnapshot ? "yes" : "no")
      + " sSnap=" + std::string(hasAcceptedStateSnapshot ? "yes" : "no"),
    "resGrow="+fmtScalar(linearResidualGrowthRatio)
      + " health=" + std::string(healthGatesPassed ? "yes" : "no")
  );
  writeRow
  (
    "Rogue dir",
    "enabled="+std::string(static_cast<bool>(optRogueDirectionSwitchEnable) ? "yes" : "no")
      + " signal=" + std::string(rogueSignal ? "yes" : "no")
      + " rejSig=" + std::string(rogueRejectSignal ? "yes" : "no")
      + " clipSig=" + std::string(rogueClipSignal ? "yes" : "no"),
    "watch="+std::string(rogueDirectionWatchActive ? "yes" : "no")
      + " grace=" + fmtCount(rogueDirectionWatchGraceRemaining)
      + " active=" + std::string(rogueDirectionGraceActive ? "yes" : "no")
      + " verified=" + std::string(rogueDirectionVerified ? "yes" : "no"),
    "switchTrig="+std::string(rogueDirectionSwitchTriggered ? "yes" : "no")
      + " applied=" + std::string(rogueDirectionSwitchApplied ? "yes" : "no")
      + " count=" + fmtCount(rogueDirectionSwitchCount)
      + "/" + fmtCount(rogueDirectionSwitchMaxCountClamped),
    "coolRem="+fmtCount(rogueDirectionSwitchCooldownRemaining)
      + " rearm=" + std::string(rogueDirectionSwitchBudgetRearmed ? "yes" : "no")
      + " usObj=" + std::string(rogueDirectionObjectiveUsable ? "yes" : "no")
      + " usCon=" + std::string(rogueDirectionConstraintUsable ? "yes" : "no")
      + " cScale=" + fmtScalar(rogueDirectionSwitchClipStepScaleApplied)
      + " last=" + fmtCount(rogueDirectionLastSwitchIteration)
      + " start=" + fmtCount(rogueDirectionWatchStartIterationReport)
  );
  writeRow
  (
    "Cont gate",
    "schedIter="+fmtCount(continuationScheduleIterationActive),
    "rawIter="+fmtCount(opt),
    "freezeOnReject="+std::string(static_cast<bool>(optFreezeContinuationOnReject) ? "yes" : "no")
      + " lowTrust=" + std::string(static_cast<bool>(optFreezeContinuationOnLowTrust) ? "yes" : "no")
      + " clipHold=" + std::string(static_cast<bool>(optContinuationClipHoldEnable) ? "yes" : "no"),
    "frozen="+std::string(reportedContinuationFrozenAny ? "yes" : "no")
      + " rej=" + std::string(continuationFrozenByReject ? "yes" : "no")
      + " low=" + std::string(reportedContinuationFrozenByLowTrust ? "yes" : "no")
      + " clip=" + std::string(continuationFrozenByClipHold ? "yes" : "no")
      + " phase=" + std::string(continuationFreezePhaseActive ? "active" : "inactive")
      + " thr=" + fmtScalar(continuationLowTrustThresholdReport)
      + " cUp=" + fmtCount(reportedContinuationAdvanceCatchUpSteps)
  );
  writeRow
  (
    "Grad check",
    "active="+std::string(gradientCheckActive ? "yes" : "no")
      + " ratio=" + fmtScalar(gradientConsistencyRatio),
    "obj="+fmtScalar(objectiveConsistencyRatio)+" pow="+fmtScalar(powerConsistencyRatio),
    "vol="+fmtScalar(volumeConsistencyRatio)
      + " fail=" + fmtCount(reportedGradientCheckFailStreak)
      + "/" + fmtCount(optGradientCheckFailWindow),
    "warn="+std::string(gradientCheckWarnTriggered ? "yes" : "no")
      + " stop=" + std::string(gradientCheckStopTriggered ? "yes" : "no")
  );
  writeRow
  (
    "Closure FD",
    "active="+std::string(closureDerivativeFDStageCheckActive ? "yes" : "no")
      + " ratio=" + fmtScalar(closureDerivativeFDStageCheckRatio),
    "mean="+fmtScalar(closureDerivativeFDStageCheckMeanRatio)
      + " max=" + fmtScalar(closureDerivativeFDStageCheckMaxRatio),
    "kMax="+fmtScalar(closureDerivativeFDStageCheckMaxRatioK)
      + " wMax=" + fmtScalar(closureDerivativeFDStageCheckMaxRatioOmega)
      + " mismatch=" + fmtScalar(closureDerivativeFDStageCheckMaxStoredMismatch),
    "samples="+fmtCount(closureDerivativeFDStageCheckSamplesUsed)
      + " warnFrac=" + fmtScalar(closureDerivativeFDStageCheckWarnFractionObserved)
      + " failFrac=" + fmtScalar(closureDerivativeFDStageCheckFailFractionObserved)
      + " fail=" + fmtCount(reportedClosureDerivativeFDStageCheckFailStreak)
      + "/" + fmtCount(optClosureDerivativeFDStageCheckFailWindow)
      + " warn=" + std::string(closureDerivativeFDStageCheckWarnTriggered ? "yes" : "no")
      + " stop=" + std::string(closureDerivativeFDStageCheckStopTriggered ? "yes" : "no")
      + " clipLimited=" + std::string(closureDerivativeFDStageCheckClipLimited ? "yes" : "no")
  );
  writeRow("Stop", "candidate="+std::string(convergenceCandidate ? "yes" : "no"), "stop="+std::string(stopOptimization ? "yes" : "no"), "reason="+std::string(optimizerStopReason.c_str()));
  writeRow("Outcome", "state="+std::string(stopOutcome), "success="+std::string(stopSuccessful ? "yes" : "no"), "failure="+std::string(stopFailure ? "yes" : "no"));

  writeSectionHeader("Design Diagnostics");
  writeRow("x bounds", "min="+fmtScalar(xLocalMin), "max="+fmtScalar(xLocalMax), "mean="+fmtScalar(xMean));
  writeRow("xp bounds", "min="+fmtScalar(xpLocalMin), "max="+fmtScalar(xpLocalMax));
  writeRow("xh bounds", "min="+fmtScalar(xhLocalMin), "max="+fmtScalar(xhLocalMax), "mean="+fmtScalar(xhMean));
  writeRow("xh phases", "gray="+fmtScalar(grayVolumeFraction), "solid="+fmtScalar(solidVolumeFraction), "fluid="+fmtScalar(fluidVolumeFraction));
  writeRow("x step", "max="+fmtScalar(maxAbsXChangeLocal), "l2="+fmtScalar(l2XChange), "l2/sqrtN="+fmtScalar(normalizedL2XChange));
  writeRow("xh step", "max="+fmtScalar(maxAbsXhChangeLocal), "l2="+fmtScalar(l2XhChange), "l2/sqrtN="+fmtScalar(normalizedL2XhChange));

  writeSectionHeader("Filter Diagnostics");
  writeRow("projection", "del="+fmtScalar(del), "eta5="+fmtScalar(eta5));
  writeRow("drho bounds", "min="+fmtScalar(drhoLocalMin), "max="+fmtScalar(drhoLocalMax));
  writeRow("obj sens", "min="+fmtScalar(fsensLocalMin), "max="+fmtScalar(fsensLocalMax));
  writeRow("power sens", "min="+fmtScalar(gsensPowerLocalMin), "max="+fmtScalar(gsensPowerLocalMax));
  writeRow("vol sens", "min="+fmtScalar(gsensVolLocalMin), "max="+fmtScalar(gsensVolLocalMax));

  writeSectionHeader("Objective And Constraints");
  writeRow("objective", "value="+fmtScalar(objectiveValue), "delta="+fmtScalar(objectiveDelta));
  writeRow("power diss", "actual="+fmtScalar(PowerDissActual), "act/base="+fmtScalar(PowerDissActualOverBase), "constraint="+fmtScalar(PowerDiss));
  writeRow("power margin", "sched="+fmtScalar(powerConstraintMargin), "final="+fmtScalar(powerConstraintFinalMargin));
  writeRow("volume", "value="+fmtScalar(V), "margin="+fmtScalar(volumeConstraintMargin));
  writeRow
  (
    "thermal env",
    "smoothExcess="+fmtScalar(thermalExcessObjective),
    "tmaxKS="+fmtScalar(thermalTmaxConstraintMargin),
    "tminKS="+fmtScalar(thermalTminConstraintMargin),
    "pen="+fmtScalar(optThermalEnvelopePenalty)
  );
  writeRow("margin delta", "power="+fmtScalar(powerConstraintMarginDelta), "powerFinal="+fmtScalar(powerConstraintFinalMarginDelta), "volume="+fmtScalar(volumeConstraintMarginDelta), "maxGuard="+fmtScalar(guardConstraintViolationDelta));
  writeRow("thermal", "meanT="+fmtScalar(MeanT), "maxT="+fmtScalar(maxT));

  writeSectionHeader("MMA Diagnostics");
  writeRow
  (
    "controls",
    "moveLimit="+fmtScalar(optimizerMoveLimitActive),
    "fixed=yes",
    "moveHits="+fmtCount(xMoveLimitHitCountLocal)
  );
  writeRow("bound clamps", "lower="+fmtCount(xLowerBoundCountLocal), "upper="+fmtCount(xUpperBoundCountLocal), "nonFiniteX="+fmtCount(xNonFiniteCountLocal));

  writeSectionHeader("Sensitivity Diagnostics");
  writeRow("objective", "min="+fmtScalar(fsensLocalMin), "max="+fmtScalar(fsensLocalMax), "l2="+fmtScalar(fsensL2), "nonFinite="+fmtCount(fsensNonFiniteCountLocal));
  writeRow("power", "min="+fmtScalar(gsensPowerLocalMin), "max="+fmtScalar(gsensPowerLocalMax), "l2="+fmtScalar(gsensPowerL2), "nonFinite="+fmtCount(gsensPowerNonFiniteCountLocal));
  writeRow("volume", "min="+fmtScalar(gsensVolLocalMin), "max="+fmtScalar(gsensVolLocalMax), "l2="+fmtScalar(gsensVolL2), "nonFinite="+fmtCount(gsensVolNonFiniteCountLocal));

  writeSectionHeader("Solver Diagnostics");
  const label effectivePrimalCorrectors = momentumCorrectors;
  const label effectiveAdjPowerCorrectors =
      optimizationCycleTurbulentOptimizationMode
    ? optimizationCycleFlowAdjointCorrectors
    : (adjointCycleFullAdjointVerification ? adjointCycleAdjointCorrectors : 40);
  const label effectiveAdjHeatCorrectors =
      optimizationCycleTurbulentOptimizationMode
    ? optimizationCycleHeatAdjointCorrectors
    : 40;
  const char* optimizationModeLabel =
      optimizationCycleTurbulentOptimizationMode
    ? "optimizationCycle"
    : (adjointCycleFullAdjointVerification ? "adjointCycle-full"
      : (adjointCycleThermalAdjointOnly ? "adjointCycle-thermal" : "baseline"));
  writeRow
  (
      "config",
      std::string("mode=") + optimizationModeLabel,
      "primal="+fmtCount(effectivePrimalCorrectors),
      "adjPower="+fmtCount(effectiveAdjPowerCorrectors),
      "adjHeat="+fmtCount(effectiveAdjHeatCorrectors)
  );
  writeRow("primal flow", "|U|max="+fmtScalar(UMaxMagLocal), "|U|mean="+fmtScalar(UMeanMag), "pMin="+fmtScalar(pLocalMin), "pMax="+fmtScalar(pLocalMax));
  writeRow("turb closure", "active="+std::string(customTurbulenceActive ? "yes" : "no"), "kMin="+fmtScalar(kGlobalMin), "kMax="+fmtScalar(kGlobalMax), "nutMax="+fmtScalar(nutGlobalMax));
  writeRow("omega field", "omegaMin="+fmtScalar(omegaGlobalMin), "omegaMax="+fmtScalar(omegaGlobalMax));
  writeRow
  (
    "porous turb",
    "aware="+std::string(porousAwareTurbulenceTransport ? "yes" : "no"),
    "prodExp="+fmtScalar(porousTurbulenceProductionExponent)
      + " prodFloor=" + fmtScalar(porousTurbulenceProductionFloor),
    "sinkS="+fmtScalar(porousTurbulenceSuppressionCoeff)
      + " sinkA=" + fmtScalar(porousTurbulenceAlphaSinkCoeff)
  );
  writeRow("clip health", "status="+std::string(turbulenceClipHealth), "fluidFrac="+fmtScalar(maxFluidClipFraction), "grayFrac="+fmtScalar(maxGrayClipFraction), "allFrac="+fmtScalar(maxClipFraction));
  writeRow("clip trigger", "severeFluidThr="+fmtScalar(severeTransitionFluidThreshold), "severeGrayThr=8.500000e-01");
  writeRow("turb clips", "kFloor="+fmtCount(kFloorCountLocal), "omegaFloor="+fmtCount(omegaFloorCountLocal), "nutZero="+fmtCount(nutZeroCountLocal), "nutCap="+fmtCount(nutCapCountLocal));
  writeRow("clip frac", "kFloor="+fmtScalar(kFloorFraction), "omegaFloor="+fmtScalar(omegaFloorFraction), "nutZero="+fmtScalar(nutZeroFraction), "nutCap="+fmtScalar(nutCapFraction));
  writeRow("phase cells", "solid="+fmtCount(solidCellCountLocal), "gray="+fmtCount(grayCellCountLocal), "fluid="+fmtCount(fluidCellCountLocal));
  writeRow("kFloor phase", "solid="+fmtScalar(kFloorSolidFraction), "gray="+fmtScalar(kFloorGrayFraction), "fluid="+fmtScalar(kFloorFluidFraction));
  writeRow("wFloor phase", "solid="+fmtScalar(omegaFloorSolidFraction), "gray="+fmtScalar(omegaFloorGrayFraction), "fluid="+fmtScalar(omegaFloorFluidFraction));
  writeRow("thermal field", "TMin="+fmtScalar(TLocalMin), "TMax="+fmtScalar(TLocalMax));
  writeRow("adjoint power", "|Ua|max="+fmtScalar(UaMaxMagLocal), "|Ua|mean="+fmtScalar(UaMeanMag), "paMin="+fmtScalar(paLocalMin), "paMax="+fmtScalar(paLocalMax));
  writeRow("adjoint heat", "|Ub|max="+fmtScalar(UbMaxMagLocal), "|Ub|mean="+fmtScalar(UbMeanMag), "pbMin="+fmtScalar(pbLocalMin), "pbMax="+fmtScalar(pbLocalMax));

  writeSectionHeader("Interpolation Diagnostics");
  writeRow("interp bounds", "alphaMin="+fmtScalar(alphaLocalMin), "alphaMax="+fmtScalar(alphaLocalMax), "DTMin="+fmtScalar(DTLocalMin), "DTMax="+fmtScalar(DTLocalMax));
  writeRow
  (
    "continuation",
    "alphaMax="+fmtScalar(alphaMax.value()),
    "qu="+fmtScalar(qu),
    "schedIter="+fmtCount(continuationScheduleIterationActive),
    "frozen="+std::string(reportedContinuationFrozenAny ? "yes" : "no")
      + " rej=" + std::string(continuationFrozenByReject ? "yes" : "no")
      + " low=" + std::string(reportedContinuationFrozenByLowTrust ? "yes" : "no")
      + " clip=" + std::string(continuationFrozenByClipHold ? "yes" : "no")
  );
  writeRow("thermal props", "heatSource="+fmtScalar(heatSource.value()), "thermalSrc="+fmtScalar(thermalSource.value()), "rhoCp="+fmtScalar(rhoCp.value()));

  std::ofstream debugJsonFile(debugOptimizerJsonFile, std::ios::app);
  debugJsonFile
    << std::scientific << std::setprecision(8)
    << "{"
    << "\"iteration\":" << opt
    << ",\"clockTime\":" << clockTime
    << ",\"health\":{"
      << "\"design\":\"" << designHealth << "\""
      << ",\"constraints\":\"" << constraintHealth << "\""
      << ",\"sensitivity\":\"" << sensitivityHealth << "\""
      << ",\"solver\":\"" << solverHealth << "\""
      << ",\"turbulenceClip\":\"" << turbulenceClipHealth << "\""
    << "}"
    << ",\"convergence\":{"
      << "\"objectiveRelativeDelta\":" << objectiveRelativeDelta
      << ",\"maxConstraintViolation\":" << maxConstraintViolation
      << ",\"maxConstraintViolationDelta\":" << maxConstraintViolationDelta
      << ",\"convergenceConstraintViolation\":" << convergenceConstraintViolation
      << ",\"guardConstraintViolation\":" << guardConstraintViolation
      << ",\"guardConstraintViolationDelta\":" << guardConstraintViolationDelta
      << ",\"guardConstraintViolationAtCountStart\":" << guardConstraintViolationAtCountStart
      << ",\"guardConstraintViolationMinSinceCountStart\":" << guardConstraintViolationMinSinceCountStart
      << ",\"guardConstraintViolationNetRise\":" << guardConstraintViolationNetRise
      << ",\"guardConstraintViolationReboundRise\":" << guardConstraintViolationReboundRise
      << ",\"maxConstraintViolationWorsening\":" << (maxConstraintViolationWorsening ? 1 : 0)
      << ",\"guardConstraintViolationWorsening\":" << (guardConstraintViolationWorsening ? 1 : 0)
      << ",\"maxConstraintViolationWorsenStreak\":" << reportedMaxViolationWorsenStreak
      << ",\"powerConstraintMarginWorsening\":" << (powerConstraintMarginWorsening ? 1 : 0)
      << ",\"powerConstraintWorsening\":" << (powerConstraintWorsening ? 1 : 0)
      << ",\"powerConstraintFinalMargin\":" << powerConstraintFinalMargin
      << ",\"powerConstraintFinalMarginDelta\":" << powerConstraintFinalMarginDelta
      << ",\"powerConstraintFinalWorsening\":" << (powerConstraintFinalWorsening ? 1 : 0)
      << ",\"powerConstraintWorsenStreak\":" << reportedPowerConstraintWorsenStreak
      << ",\"constraintWorsenGuardEnabled\":" << (constraintWorsenGuardEnabled ? 1 : 0)
      << ",\"constraintWorsenCountingActive\":" << (constraintWorsenCountingActive ? 1 : 0)
      << ",\"constraintWorsenNetRiseTol\":" << effectiveConstraintWorsenNetRiseTol
      << ",\"constraintWorsenNetRiseMinAge\":" << optConstraintWorsenNetRiseMinAge
      << ",\"constraintWorsenNetRiseEffectiveMinAge\":" << effectiveConstraintWorsenNetRiseMinAge
      << ",\"constraintWorsenNetRiseAge\":" << guardConstraintViolationCountAge
      << ",\"constraintWorsenNetRiseMaxAgeInRun\":" << maxGuardConstraintViolationAgeInRun
      << ",\"constraintWorsenNetRiseDelta\":" << guardConstraintViolationDelta
      << ",\"constraintWorsenNetRiseTriggered\":" << (constraintWorsenNetRiseTriggered ? 1 : 0)
      << ",\"constraintWorsenReboundWindow\":" << effectiveConstraintReboundWindow
      << ",\"constraintWorsenReboundMinAge\":" << effectiveConstraintReboundMinAge
      << ",\"constraintWorsenReboundRiseTol\":" << effectiveConstraintReboundRiseTol
      << ",\"constraintWorsenReboundRawTriggered\":" << (constraintWorsenReboundTriggered ? 1 : 0)
      << ",\"constraintWorsenReboundCorroborated\":" << (constraintWorsenReboundCorroborated ? 1 : 0)
      << ",\"constraintWorsenReboundNetRiseCorroborationTol\":" << constraintReboundNetRiseCorroborationTol
      << ",\"constraintWorsenReboundTriggered\":" << (constraintWorsenReboundStopTriggered ? 1 : 0)
      << ",\"constraintWorsenScheduleSettled\":" << (powerConstraintScheduleSettled ? 1 : 0)
      << ",\"constraintWorsenScheduleTightening\":" << (powerConstraintScheduleTightening ? 1 : 0)
      << ",\"constraintWorsenScheduleReady\":" << (constraintWorsenScheduleReady ? 1 : 0)
      << ",\"constraintWorsenScheduleReadyTransition\":" << (constraintWorsenScheduleReadyTransition ? 1 : 0)
      << ",\"constraintWorsenScheduleIterationsToSettle\":" << powerConstraintScheduleIterationsToSettle
      << ",\"optimizationIterationsRemaining\":" << optimizationIterationsRemaining
      << ",\"constraintScheduleHorizonLagIterations\":" << scheduleHorizonLagIterations
      << ",\"constraintScheduleHorizonLagTolerance\":" << constraintScheduleHorizonLagTolerance
      << ",\"constraintWorsenScheduleWillSettleInRun\":" << (powerConstraintScheduleWillSettleInRun ? 1 : 0)
      << ",\"constraintScheduleHorizonGuardEnabled\":" << (static_cast<bool>(optConstraintScheduleHorizonGuard) ? 1 : 0)
      << ",\"constraintScheduleHorizonGuardTriggered\":" << (constraintScheduleHorizonGuardTriggered ? 1 : 0)
      << ",\"powerConstraintScheduleDenom\":" << powerConstraintScheduleDenomCurrent
      << ",\"constraintWorsenGuardTriggered\":" << (constraintWorsenGuardTriggered ? 1 : 0)
      << ",\"objectiveWorsening\":" << (objectiveWorsening ? 1 : 0)
      << ",\"objectiveWorsenTol\":" << optObjectiveWorsenTol
      << ",\"objectiveWorsenStreak\":" << reportedObjectiveWorsenStreak
      << ",\"objectiveWorsenWindow\":" << optObjectiveWorsenWindow
      << ",\"objectiveWorsenGuardTriggered\":" << (objectiveWorsenGuardTriggered ? 1 : 0)
      << ",\"clipGuardEnabled\":" << (clipGuardEnabled ? 1 : 0)
      << ",\"clipGuardActive\":" << (clipGuardActive ? 1 : 0)
      << ",\"clipGuardNoProgressActive\":" << (clipGuardNoProgressActive ? 1 : 0)
      << ",\"clipGuardFluidFractionThreshold\":" << clipGuardFluidFractionThreshold
      << ",\"clipGuardStreak\":" << reportedClipGuardStreak
      << ",\"clipGuardNoProgressStreak\":" << reportedClipGuardNoProgressStreak
      << ",\"clipGuardNoProgressTol\":" << clipGuardNoProgressTolClamped
      << ",\"clipGuardWindow\":" << clipGuardWindowDisabled
      << ",\"clipGuardTriggered\":" << (clipGuardTriggered ? 1 : 0)
      << ",\"clipGuardRejectMode\":" << (clipGuardRejectMode ? 1 : 0)
      << ",\"clipGuardForcedReject\":" << (clipGuardForcedReject ? 1 : 0)
      << ",\"clipHardRejectFluidFraction\":" << clipHardRejectFluidFractionClamped
      << ",\"clipHardRejectGrayFraction\":" << clipHardRejectGrayFractionClamped
      << ",\"clipHardRejectFluidJump\":" << clipHardRejectFluidJumpClamped
      << ",\"clipHardRejectByFraction\":" << (clipHardRejectByFraction ? 1 : 0)
      << ",\"clipHardRejectByJump\":" << (clipHardRejectByJump ? 1 : 0)
      << ",\"clipHardForcedReject\":" << (clipHardForcedReject ? 1 : 0)
      << ",\"clipForcedReject\":" << (clipForcedReject ? 1 : 0)
      << ",\"objectiveAbsValue\":" << objectiveAbsValue
      << ",\"objectiveAbsMax\":" << optObjectiveAbsMax
      << ",\"objectiveEnvelopeViolated\":" << (objectiveEnvelopeViolated ? 1 : 0)
      << ",\"thermalFieldAbsMaxObserved\":" << thermalFieldAbsMaxObserved
      << ",\"thermalFieldAbsMax\":" << optThermalFieldAbsMax
      << ",\"thermalFieldEnvelopeViolated\":" << (thermalFieldEnvelopeViolated ? 1 : 0)
      << ",\"fieldEnvelopeForcedReject\":" << (fieldEnvelopeForcedReject ? 1 : 0)
      << ",\"constraintNonDeteriorationAccepted\":" << (constraintNonDeteriorationAccepted ? 1 : 0)
      << ",\"clipJumpHealthGatePassed\":" << (clipJumpHealthGatePassed ? 1 : 0)
      << ",\"residualGrowthGatePassed\":" << (residualGrowthGatePassed ? 1 : 0)
      << ",\"healthGatesPassed\":" << (healthGatesPassed ? 1 : 0)
      << ",\"maxLinearFinalResidual\":" << maxLinearFinalResidual
      << ",\"linearResidualGrowthRatio\":" << linearResidualGrowthRatio
      << ",\"sensitivityQualityScore\":" << sensitivityQualityScoreCurrent
      << ",\"sensitivityQualityWarn\":" << (sensitivityQualityWarnCurrent ? 1 : 0)
      << ",\"sensitivityQualityFail\":" << (sensitivityQualityFailCurrent ? 1 : 0)
      << ",\"sensitivityQualityHardFail\":" << (sensitivityQualityHardFailCurrent ? 1 : 0)
      << ",\"sensitivityQualityStopTriggered\":" << (sensitivityQualityStopTriggered ? 1 : 0)
      << ",\"acceptanceEnabled\":0"
      << ",\"stepAccepted\":" << (meritAccepted ? 1 : 0)
      << ",\"stepRejected\":" << (meritStepRejected ? 1 : 0)
      << ",\"rejectStreak\":" << meritRejectStreak
      << ",\"globalizationRejectStopEligible\":" << (globalizationRejectStopEligible ? 1 : 0)
      << ",\"rejectCountsTowardStop\":" << (meritRejectCountsTowardStop ? 1 : 0)
      << ",\"scheduleHorizonCatchUpNeeded\":" << (scheduleHorizonCatchUpNeeded ? 1 : 0)
      << ",\"scheduleHorizonNoSlackUnderFreeze\":" << (scheduleHorizonNoSlackUnderFreeze ? 1 : 0)
      << ",\"scheduleHorizonDeadlockOnReject\":" << (scheduleHorizonDeadlockOnReject ? 1 : 0)
      << ",\"deadlockRecoveryBlockedByGlobalization\":"
      << (deadlockRecoveryBlockedByGlobalization ? 1 : 0)
      << ",\"continuationDeadlockRecoveryEnabled\":" << (static_cast<bool>(optContinuationDeadlockRecoveryEnable) ? 1 : 0)
      << ",\"continuationDeadlockRecoveryCount\":" << reportedContinuationDeadlockRecoveryCount
      << ",\"continuationDeadlockRecoveryMaxOverrides\":" << optContinuationDeadlockRecoveryMaxOverrides
      << ",\"continuationDeadlockRecoveryMaxCatchUpSteps\":" << optContinuationDeadlockRecoveryMaxCatchUpSteps
      << ",\"continuationDeadlockRecoveryLagBuffer\":" << optContinuationDeadlockRecoveryLagBuffer
      << ",\"continuationDeadlockCatchUpRequest\":" << reportedDeadlockCatchUpRequest
      << ",\"continuationAdvanceOverrideOnReject\":" << (continuationAdvanceOverrideOnReject ? 1 : 0)
      << ",\"continuationAdvanceCatchUpSteps\":" << reportedContinuationAdvanceCatchUpSteps
      << ",\"continuationClipBackoffStepsApplied\":" << continuationClipBackoffStepsApplied
      << ",\"trustMoveScale\":" << meritTrustMoveScale
      << ",\"trustMoveMin\":" << meritTrustMoveMinClamped
      << ",\"recoveryRequest\":" << (meritRecoveryRequest ? 1 : 0)
      << ",\"rogueDirectionSwitchEnabled\":" << (static_cast<bool>(optRogueDirectionSwitchEnable) ? 1 : 0)
      << ",\"rogueRejectSignal\":" << (rogueRejectSignal ? 1 : 0)
      << ",\"rogueClipSignal\":" << (rogueClipSignal ? 1 : 0)
      << ",\"rogueSignal\":" << (rogueSignal ? 1 : 0)
      << ",\"rogueDirectionWatchActive\":" << (rogueDirectionWatchActive ? 1 : 0)
      << ",\"rogueDirectionGraceRemaining\":" << rogueDirectionWatchGraceRemaining
      << ",\"rogueDirectionGraceActive\":" << (rogueDirectionGraceActive ? 1 : 0)
      << ",\"rogueDirectionVerified\":" << (rogueDirectionVerified ? 1 : 0)
      << ",\"rogueDirectionWatchStartIteration\":" << rogueDirectionWatchStartIteration
      << ",\"rogueDirectionWatchStartIterationReport\":" << rogueDirectionWatchStartIterationReport
      << ",\"rogueDirectionSwitchTriggered\":" << (rogueDirectionSwitchTriggered ? 1 : 0)
      << ",\"rogueDirectionSwitchApplied\":" << (rogueDirectionSwitchApplied ? 1 : 0)
      << ",\"rogueDirectionObjectiveUsable\":" << (rogueDirectionObjectiveUsable ? 1 : 0)
      << ",\"rogueDirectionConstraintUsable\":" << (rogueDirectionConstraintUsable ? 1 : 0)
      << ",\"rogueDirectionSwitchBudgetRearmed\":" << (rogueDirectionSwitchBudgetRearmed ? 1 : 0)
      << ",\"rogueDirectionSwitchCount\":" << rogueDirectionSwitchCount
      << ",\"rogueDirectionSwitchMaxCount\":" << rogueDirectionSwitchMaxCountClamped
      << ",\"rogueDirectionSwitchCooldown\":" << rogueDirectionSwitchCooldownClamped
      << ",\"rogueDirectionSwitchCooldownRemaining\":" << rogueDirectionSwitchCooldownRemaining
      << ",\"rogueDirectionLastSwitchIteration\":" << rogueDirectionLastSwitchIteration
      << ",\"rogueDirectionSwitchClipDampFluidFraction\":"
      << rogueDirectionSwitchClipDampFluidFractionClamped
      << ",\"rogueDirectionSwitchClipDampGrayFraction\":"
      << rogueDirectionSwitchClipDampGrayFractionClamped
      << ",\"rogueDirectionSwitchClipDampScale\":"
      << rogueDirectionSwitchClipDampScaleClamped
      << ",\"rogueDirectionSwitchSevereClipFluidFraction\":"
      << rogueDirectionSwitchSevereClipFluidFractionClamped
      << ",\"rogueDirectionSwitchSevereClipGrayFraction\":"
      << rogueDirectionSwitchSevereClipGrayFractionClamped
      << ",\"rogueDirectionSwitchSevereClipScale\":"
      << rogueDirectionSwitchSevereClipScaleClamped
      << ",\"rogueDirectionSwitchClipStepScaleApplied\":"
      << rogueDirectionSwitchClipStepScaleApplied
      << ",\"rogueDirectionSwitchStepAmplitudeApplied\":"
      << rogueDirectionSwitchStepAmplitudeApplied
      << ",\"freezeContinuationOnRejectEnabled\":" << (static_cast<bool>(optFreezeContinuationOnReject) ? 1 : 0)
      << ",\"freezeContinuationOnLowTrustEnabled\":" << (static_cast<bool>(optFreezeContinuationOnLowTrust) ? 1 : 0)
      << ",\"continuationLowTrustThreshold\":" << continuationLowTrustThresholdReport
      << ",\"rejectStopRequireConstraintWindow\":" << (static_cast<bool>(optRejectStopRequireConstraintWindow) ? 1 : 0)
      << ",\"earlyRecoveryOnRejectEnabled\":" << (static_cast<bool>(optEarlyRecoveryOnGlobalizationReject) ? 1 : 0)
      << ",\"recoveryEnabled\":" << (static_cast<bool>(optRecoveryEnable) ? 1 : 0)
      << ",\"recoveryActive\":" << (constraintRecoveryActive ? 1 : 0)
      << ",\"recoveryWorsenStreak\":" << optRecoveryWorsenStreak
      << ",\"recoveryMoveLimitScale\":" << optRecoveryMoveLimitScale
      << ",\"recoveryContinuationRampScale\":" << optRecoveryContinuationRampScale
      << ",\"recoveryContinuationRampScaleActive\":" << continuationRampScaleActive
      << ",\"recoveryDelRampStepActive\":" << delRampStepActive
      << ",\"continuationScheduleIteration\":" << continuationScheduleIterationActive
      << ",\"continuationFrozenByReject\":" << (continuationFrozenByReject ? 1 : 0)
      << ",\"continuationFrozenByLowTrust\":" << (reportedContinuationFrozenByLowTrust ? 1 : 0)
      << ",\"continuationFrozenByClipHold\":" << (continuationFrozenByClipHold ? 1 : 0)
      << ",\"continuationFrozen\":" << (reportedContinuationFrozenAny ? 1 : 0)
      << ",\"continuationFreezePhaseActive\":" << (continuationFreezePhaseActive ? 1 : 0)
      << ",\"gradientCheckEnabled\":" << (static_cast<bool>(optGradientCheckEnable) ? 1 : 0)
      << ",\"gradientCheckActive\":" << (gradientCheckActive ? 1 : 0)
      << ",\"gradientCheckMinStepL2\":" << optGradientCheckMinStepL2
      << ",\"gradientCheckWarnRatio\":" << optGradientCheckWarnRatio
      << ",\"gradientCheckFailRatio\":" << optGradientCheckFailRatio
      << ",\"gradientCheckWarnTriggered\":" << (gradientCheckWarnTriggered ? 1 : 0)
      << ",\"gradientCheckFailTriggered\":" << (gradientCheckFailTriggered ? 1 : 0)
      << ",\"gradientCheckFailStreak\":" << reportedGradientCheckFailStreak
      << ",\"gradientCheckFailWindow\":" << optGradientCheckFailWindow
      << ",\"gradientCheckStopTriggered\":" << (gradientCheckStopTriggered ? 1 : 0)
      << ",\"closureDerivativeFDEnabled\":" << (static_cast<bool>(optClosureDerivativeFDEnable) ? 1 : 0)
      << ",\"closureDerivativeFDRelStep\":" << optClosureDerivativeFDRelStep
      << ",\"closureDerivativeFDStageCheckEnabled\":" << (static_cast<bool>(optClosureDerivativeFDStageCheckEnable) ? 1 : 0)
      << ",\"closureDerivativeFDStageCheckActive\":" << (closureDerivativeFDStageCheckActive ? 1 : 0)
      << ",\"closureDerivativeFDStageCheckSamples\":" << optClosureDerivativeFDStageCheckSamples
      << ",\"closureDerivativeFDStageCheckSamplesUsed\":" << closureDerivativeFDStageCheckSamplesUsed
      << ",\"closureDerivativeFDStageCheckWarnRatio\":" << optClosureDerivativeFDStageCheckWarnRatio
      << ",\"closureDerivativeFDStageCheckWarnFractionRequired\":"
      << optClosureDerivativeFDStageCheckWarnFraction
      << ",\"closureDerivativeFDStageCheckFailRatio\":" << optClosureDerivativeFDStageCheckFailRatio
      << ",\"closureDerivativeFDStageCheckFailFractionRequired\":"
      << optClosureDerivativeFDStageCheckFailFraction
      << ",\"closureDerivativeFDStageCheckWarnTriggered\":" << (closureDerivativeFDStageCheckWarnTriggered ? 1 : 0)
      << ",\"closureDerivativeFDStageCheckFailTriggered\":" << (closureDerivativeFDStageCheckFailTriggered ? 1 : 0)
      << ",\"closureDerivativeFDStageCheckFailStreak\":" << reportedClosureDerivativeFDStageCheckFailStreak
      << ",\"closureDerivativeFDStageCheckFailWindow\":" << optClosureDerivativeFDStageCheckFailWindow
      << ",\"closureDerivativeFDStageCheckStopTriggered\":" << (closureDerivativeFDStageCheckStopTriggered ? 1 : 0)
      << ",\"closureDerivativeFDStageCheckClipLimited\":" << (closureDerivativeFDStageCheckClipLimited ? 1 : 0)
      << ",\"closureDerivativeFDStageCheckRatio\":" << closureDerivativeFDStageCheckRatio
      << ",\"closureDerivativeFDStageCheckMeanRatio\":" << closureDerivativeFDStageCheckMeanRatio
      << ",\"closureDerivativeFDStageCheckMaxRatio\":" << closureDerivativeFDStageCheckMaxRatio
      << ",\"closureDerivativeFDStageCheckMaxRatioK\":" << closureDerivativeFDStageCheckMaxRatioK
      << ",\"closureDerivativeFDStageCheckMaxRatioOmega\":" << closureDerivativeFDStageCheckMaxRatioOmega
      << ",\"closureDerivativeFDStageCheckMaxStoredMismatch\":" << closureDerivativeFDStageCheckMaxStoredMismatch
      << ",\"closureDerivativeFDStageCheckWarnFractionObserved\":"
      << closureDerivativeFDStageCheckWarnFractionObserved
      << ",\"closureDerivativeFDStageCheckFailFractionObserved\":"
      << closureDerivativeFDStageCheckFailFractionObserved
      << ",\"gradientConsistencyRatio\":" << gradientConsistencyRatio
      << ",\"gradientObjectiveRatio\":" << objectiveConsistencyRatio
      << ",\"gradientPowerRatio\":" << powerConsistencyRatio
      << ",\"gradientVolumeRatio\":" << volumeConsistencyRatio
      << ",\"predictedObjectiveDelta\":" << predictedObjectiveDelta
      << ",\"predictedPowerDelta\":" << predictedPowerDelta
      << ",\"predictedVolumeDelta\":" << predictedVolumeDelta
      << ",\"predictedThermalTmaxDelta\":" << predictedThermalTmaxDelta
      << ",\"predictedThermalTminDelta\":" << predictedThermalTminDelta
      << ",\"constraintWorsenMinIterations\":" << optConstraintWorsenMinIterations
      << ",\"constraintWorsenScheduleIteration\":" << scheduleIterationCurrent
      << ",\"constraintWorsenWindow\":" << optConstraintWorsenWindow
      << ",\"constraintWorsenTol\":" << optConstraintWorsenTol
      << ",\"xStepMax\":" << maxAbsXChangeLocal
      << ",\"xhStepMax\":" << maxAbsXhChangeLocal
      << ",\"xStepL2OverSqrtN\":" << normalizedL2XChange
      << ",\"xhStepL2OverSqrtN\":" << normalizedL2XhChange
      << ",\"streak\":" << reportedOptimizerConvergedStreak
      << ",\"stallWindow\":" << optStallWindow
      << ",\"candidate\":" << (convergenceCandidate ? 1 : 0)
      << ",\"stop\":" << (stopOptimization ? 1 : 0)
      << ",\"stopSuccessful\":" << (stopSuccessful ? 1 : 0)
      << ",\"stopFailure\":" << (stopFailure ? 1 : 0)
      << ",\"stopOutcome\":\"" << stopOutcome << "\""
      << ",\"reason\":\"" << optimizerStopReason << "\""
    << "}"
    << ",\"design\":{"
      << "\"xMin\":" << xLocalMin
      << ",\"xMax\":" << xLocalMax
      << ",\"xMean\":" << xMean
      << ",\"xpMin\":" << xpLocalMin
      << ",\"xpMax\":" << xpLocalMax
      << ",\"xhMin\":" << xhLocalMin
      << ",\"xhMax\":" << xhLocalMax
      << ",\"xhMean\":" << xhMean
      << ",\"xhGrayVolumeFraction\":" << grayVolumeFraction
      << ",\"xhSolidVolumeFraction\":" << solidVolumeFraction
      << ",\"xhFluidVolumeFraction\":" << fluidVolumeFraction
      << ",\"xStepMax\":" << maxAbsXChangeLocal
      << ",\"xStepL2\":" << l2XChange
      << ",\"xhStepMax\":" << maxAbsXhChangeLocal
      << ",\"xhStepL2\":" << l2XhChange
    << "}"
    << ",\"filter\":{"
      << "\"del\":" << del
      << ",\"eta5\":" << eta5
      << ",\"drhoMin\":" << drhoLocalMin
      << ",\"drhoMax\":" << drhoLocalMax
      << ",\"objectiveSensitivityFilteredMin\":" << fsensLocalMin
      << ",\"objectiveSensitivityFilteredMax\":" << fsensLocalMax
      << ",\"powerConstraintSensitivityFilteredMin\":" << gsensPowerLocalMin
      << ",\"powerConstraintSensitivityFilteredMax\":" << gsensPowerLocalMax
      << ",\"volumeConstraintSensitivityFilteredMin\":" << gsensVolLocalMin
      << ",\"volumeConstraintSensitivityFilteredMax\":" << gsensVolLocalMax
    << "}"
    << ",\"objective\":{"
      << "\"name\":\"" << objectiveLabel << "\""
      << ",\"value\":" << objectiveValue
      << ",\"delta\":" << objectiveDelta
      << ",\"powerDissipationActual\":" << PowerDissActual
      << ",\"powerDissipationActualOverPowerDiss0\":" << PowerDissActualOverBase
      << ",\"powerDissipationConstraintValue\":" << PowerDiss
      << ",\"powerDissipationConstraintMargin\":" << powerConstraintMargin
      << ",\"powerDissipationConstraintMarginDelta\":" << powerConstraintMarginDelta
      << ",\"powerDissipationFinalConstraintMargin\":" << powerConstraintFinalMargin
      << ",\"powerDissipationFinalConstraintMarginDelta\":" << powerConstraintFinalMarginDelta
      << ",\"volumeConstraintValue\":" << V
      << ",\"volumeConstraintMargin\":" << volumeConstraintMargin
      << ",\"volumeConstraintMarginDelta\":" << volumeConstraintMarginDelta
      << ",\"thermalExcessObjective\":" << thermalExcessObjective
      << ",\"thermalTmaxKSConstraintMargin\":" << thermalTmaxConstraintMargin
      << ",\"thermalTmaxKSConstraintMarginDelta\":" << thermalTmaxConstraintMarginDelta
      << ",\"thermalTminKSConstraintMargin\":" << thermalTminConstraintMargin
      << ",\"thermalTminKSConstraintMarginDelta\":" << thermalTminConstraintMarginDelta
      << ",\"meanT\":" << MeanT
      << ",\"maxT\":" << maxT
    << "}"
    << ",\"mma\":{"
      << "\"iteration\":" << opt
      << ",\"moveLimit\":" << optimizerMoveLimitActive
      << ",\"moveLimitMode\":\"fixed\""
      << ",\"asymInc\":" << mma_inc
      << ",\"asymDec\":" << mma_dec
      << ",\"recoveryActive\":" << (constraintRecoveryActive ? 1 : 0)
      << ",\"trustMoveScale\":" << meritTrustMoveScale
      << ",\"stepRejected\":" << (meritStepRejected ? 1 : 0)
      << ",\"moveLimitHitCount\":" << xMoveLimitHitCountLocal
      << ",\"lowerBoundCount\":" << xLowerBoundCountLocal
      << ",\"upperBoundCount\":" << xUpperBoundCountLocal
      << ",\"designNonFiniteCount\":" << xNonFiniteCountLocal
    << "}"
    << ",\"sensitivity\":{"
      << "\"objectiveMin\":" << fsensLocalMin
      << ",\"objectiveMax\":" << fsensLocalMax
      << ",\"objectiveL2\":" << fsensL2
      << ",\"objectiveNonFiniteCount\":" << fsensNonFiniteCountLocal
      << ",\"powerConstraintMin\":" << gsensPowerLocalMin
      << ",\"powerConstraintMax\":" << gsensPowerLocalMax
      << ",\"powerConstraintL2\":" << gsensPowerL2
      << ",\"powerConstraintNonFiniteCount\":" << gsensPowerNonFiniteCountLocal
      << ",\"volumeConstraintMin\":" << gsensVolLocalMin
      << ",\"volumeConstraintMax\":" << gsensVolLocalMax
      << ",\"volumeConstraintL2\":" << gsensVolL2
      << ",\"volumeConstraintNonFiniteCount\":" << gsensVolNonFiniteCountLocal
      << ",\"qualityScore\":" << sensitivityQualityScoreCurrent
      << ",\"qualityWarn\":" << (sensitivityQualityWarnCurrent ? 1 : 0)
      << ",\"qualityFail\":" << (sensitivityQualityFailCurrent ? 1 : 0)
      << ",\"qualityHardFail\":" << (sensitivityQualityHardFailCurrent ? 1 : 0)
      << ",\"qualityFailStreak\":" << candidateSensitivityQualityFailStreak
      << ",\"reuseAcceptedGradient\":" << (sensitivityQualityReuseAcceptedGradient ? 1 : 0)
      << ",\"hasAcceptedGradientSnapshot\":" << (hasAcceptedGradientSnapshot ? 1 : 0)
      << ",\"hasAcceptedStateSnapshot\":" << (hasAcceptedStateSnapshot ? 1 : 0)
    << "}"
    << ",\"solver\":{"
      << "\"primalFlowSweeps\":" << effectivePrimalCorrectors
      << ",\"adjointFlowPowerSweeps\":" << effectiveAdjPowerCorrectors
      << ",\"adjointFlowHeatSweeps\":" << effectiveAdjHeatCorrectors
      << ",\"uMaxMag\":" << UMaxMagLocal
      << ",\"uMeanMag\":" << UMeanMag
      << ",\"pMin\":" << pLocalMin
      << ",\"pMax\":" << pLocalMax
      << ",\"tMin\":" << TLocalMin
      << ",\"tMax\":" << TLocalMax
      << ",\"kMin\":" << kGlobalMin
      << ",\"kMax\":" << kGlobalMax
      << ",\"omegaMin\":" << omegaGlobalMin
      << ",\"omegaMax\":" << omegaGlobalMax
      << ",\"nutMax\":" << nutGlobalMax
      << ",\"porousAwareTurbulenceTransport\":"
      << (porousAwareTurbulenceTransport ? 1 : 0)
      << ",\"porousTurbulenceProductionExponent\":"
      << porousTurbulenceProductionExponent
      << ",\"porousTurbulenceProductionFloor\":"
      << porousTurbulenceProductionFloor
      << ",\"porousTurbulenceSuppressionCoeff\":"
      << porousTurbulenceSuppressionCoeff
      << ",\"porousTurbulenceAlphaSinkCoeff\":"
      << porousTurbulenceAlphaSinkCoeff
      << ",\"kFloorCount\":" << kFloorCountLocal
      << ",\"omegaFloorCount\":" << omegaFloorCountLocal
      << ",\"nutZeroCount\":" << nutZeroCountLocal
      << ",\"nutCapCount\":" << nutCapCountLocal
      << ",\"solidCellCount\":" << solidCellCountLocal
      << ",\"grayCellCount\":" << grayCellCountLocal
      << ",\"fluidCellCount\":" << fluidCellCountLocal
      << ",\"kFloorSolidCount\":" << kFloorSolidCountLocal
      << ",\"kFloorGrayCount\":" << kFloorGrayCountLocal
      << ",\"kFloorFluidCount\":" << kFloorFluidCountLocal
      << ",\"omegaFloorSolidCount\":" << omegaFloorSolidCountLocal
      << ",\"omegaFloorGrayCount\":" << omegaFloorGrayCountLocal
      << ",\"omegaFloorFluidCount\":" << omegaFloorFluidCountLocal
      << ",\"kFloorFraction\":" << kFloorFraction
      << ",\"omegaFloorFraction\":" << omegaFloorFraction
      << ",\"nutZeroFraction\":" << nutZeroFraction
      << ",\"nutCapFraction\":" << nutCapFraction
      << ",\"kFloorSolidFraction\":" << kFloorSolidFraction
      << ",\"kFloorGrayFraction\":" << kFloorGrayFraction
      << ",\"kFloorFluidFraction\":" << kFloorFluidFraction
      << ",\"omegaFloorSolidFraction\":" << omegaFloorSolidFraction
      << ",\"omegaFloorGrayFraction\":" << omegaFloorGrayFraction
      << ",\"omegaFloorFluidFraction\":" << omegaFloorFluidFraction
      << ",\"clipHealth\":\"" << turbulenceClipHealth << "\""
      << ",\"clipMaxFraction\":" << maxClipFraction
      << ",\"clipFluidFraction\":" << maxFluidClipFraction
      << ",\"clipGrayFraction\":" << maxGrayClipFraction
      << ",\"clipTransitionFraction\":" << maxTransitionClipFraction
      << ",\"clipObsThreshold\":" << turbulenceClipObsFraction
      << ",\"clipWarnThreshold\":" << turbulenceClipWarnFraction
      << ",\"clipSevereTransitionFluidThreshold\":" << severeTransitionFluidThreshold
      << ",\"clipSevereTransition\":" << (severeTransitionClip ? 1 : 0)
      << ",\"uaMaxMag\":" << UaMaxMagLocal
      << ",\"uaMeanMag\":" << UaMeanMag
      << ",\"paMin\":" << paLocalMin
      << ",\"paMax\":" << paLocalMax
      << ",\"ubMaxMag\":" << UbMaxMagLocal
      << ",\"ubMeanMag\":" << UbMeanMag
      << ",\"pbMin\":" << pbLocalMin
      << ",\"pbMax\":" << pbLocalMax
    << "}"
    << ",\"interpolation\":{"
      << "\"alphaMin\":" << alphaLocalMin
      << ",\"alphaMax\":" << alphaLocalMax
      << ",\"dtMin\":" << DTLocalMin
      << ",\"dtMax\":" << DTLocalMax
      << ",\"alphaMaxSchedule\":" << alphaMax.value()
      << ",\"qu\":" << qu
      << ",\"scheduleIteration\":" << continuationScheduleIterationActive
      << ",\"rawIteration\":" << opt
      << ",\"frozenByReject\":" << (continuationFrozenByReject ? 1 : 0)
      << ",\"frozenByLowTrust\":" << (reportedContinuationFrozenByLowTrust ? 1 : 0)
      << ",\"frozenByClipHold\":" << (continuationFrozenByClipHold ? 1 : 0)
      << ",\"frozen\":" << (reportedContinuationFrozenAny ? 1 : 0)
      << ",\"heatSource\":" << heatSource.value()
      << ",\"thermalSource\":" << thermalSource.value()
      << ",\"rhoCp\":" << rhoCp.value()
    << "}"
    << "}\n";
}

previousConstraintWorsenScheduleReady = constraintWorsenScheduleReady;
hasPreviousConstraintWorsenScheduleReady = true;
