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
const scalar convergenceConstraintViolation = maxConstraintViolation;
const scalar guardConstraintViolationNetRise = guardConstraintViolationDelta;
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
const bool clipHardSafetyStop = false;
const bool clipSafetyStop = false;
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
const scalar predictedObjectiveReduction = -predictedObjectiveDelta;
const scalar actualObjectiveReduction =
  hasPreviousObjectiveValue ? previousObjectiveValue - objectiveValue : 0.0;
const bool trustRatioDefined =
  hasPreviousObjectiveValue && predictedObjectiveReduction > SMALL;
const scalar trustRatio =
  trustRatioDefined
  ? actualObjectiveReduction/Foam::max(predictedObjectiveReduction, SMALL)
  : 1.0;
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
const bool fieldEnvelopeStop =
  optimizationCycleTurbulentOptimizationMode
  && (objectiveEnvelopeViolated || thermalFieldEnvelopeViolated);
const scalar residualGrowthAcceptLimit =
  optDiagnosticMode
  ? scalar(1.2)
  : scalar(2.0);
const bool residualGrowthWithinLimit =
     !hasPreviousMaxLinearFinalResidual
  || maxLinearFinalResidual <= linearSolveWarnResidual
  || linearResidualGrowthRatio <= residualGrowthAcceptLimit;
const bool healthChecksPassed =
  !fieldEnvelopeStop;
const scalar objectiveConsistencyRatio =
  Foam::mag(objectiveDelta - predictedObjectiveDelta)
  /Foam::max(Foam::max(Foam::mag(objectiveDelta), Foam::mag(predictedObjectiveDelta)), SMALL);
const scalar powerConsistencyRatio =
  Foam::mag(powerConstraintMarginDelta - predictedPowerDelta)
  /Foam::max(Foam::max(Foam::mag(powerConstraintMarginDelta), Foam::mag(predictedPowerDelta)), SMALL);
const scalar volumeConsistencyRatio =
  Foam::mag(volumeConstraintMarginDelta - predictedVolumeDelta)
  /Foam::max(Foam::max(Foam::mag(volumeConstraintMarginDelta), Foam::mag(predictedVolumeDelta)), SMALL);
const scalar gradientConsistencyRatio =
  Foam::max(objectiveConsistencyRatio, Foam::max(powerConsistencyRatio, volumeConsistencyRatio));

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
optimizerConvergedStreak = candidateOptimizerConvergedStreak;

if (fieldEnvelopeStop)
{
  stopOptimization = true;
  optimizerStopReason = word("fieldEnvelopeViolation");
}
else if (clipHardSafetyStop)
{
  stopOptimization = true;
  optimizerStopReason = word("clipHardViolation");
}
else if (!healthChecksPassed)
{
  stopOptimization = true;
  optimizerStopReason = word("healthViolation");
}
else
{
  stopOptimization = optimizerConvergedStreak >= optStallWindow;
  optimizerStopReason = stopOptimization ? word("converged") : word("running");
}

if (optDiagnosticMode && !residualGrowthWithinLimit)
{
  stopOptimization = true;
  optimizerStopReason = word("diagnosticResidualGrowth");
}

const label reportedOptimizerConvergedStreak = optimizerConvergedStreak;
previousTrustRatio = trustRatio;
hasPreviousTrustRatio = trustRatioDefined;
previousMaxLinearFinalResidual = maxLinearFinalResidual;
hasPreviousMaxLinearFinalResidual = true;
const bool firstFailureHealthEvent =
     stopOptimization
  && optimizerStopReason != word("converged");
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

const std::string flowAdjStopReason(adjointCycleAdjointStopReason.c_str());
const std::string heatAdjStopReason(optimizationCycleHeatAdjointStopReason.c_str());
const bool flowMagnitudeGuardHit =
  flowAdjStopReason.find("magnitude-guard") != std::string::npos;
const bool heatMagnitudeGuardHit =
  heatAdjStopReason.find("magnitude-guard") != std::string::npos;
const bool flowNoImproveHit =
  flowAdjStopReason.find("no-improve") != std::string::npos;
const bool heatNoImproveHit =
  heatAdjStopReason.find("no-improve") != std::string::npos;
const bool flowGrowthHit =
  flowAdjStopReason.find("growth-diverged") != std::string::npos;
const bool flowRollbackHit = adjointCycleAdjointRolledBack;
const bool heatRollbackHit = optimizationCycleHeatAdjointRolledBack;
const bool proModeCatastrophicAdjointEvent =
    flowRollbackHit
 && heatRollbackHit
 && (flowMagnitudeGuardHit || heatMagnitudeGuardHit || adjointCycleAdjointFlowDiverged);
if (turbulenceProModeActive && proModeAdjointStabilityPackActive)
{
  if (((opt - 1) % proModeRollbackDiagnosticWindow) == 0)
  {
    proModeFlowRollbackWindowCount = 0;
    proModeHeatRollbackWindowCount = 0;
    proModeFlowMagnitudeGuardWindowCount = 0;
    proModeHeatMagnitudeGuardWindowCount = 0;
    proModeFlowNoImproveWindowCount = 0;
    proModeHeatNoImproveWindowCount = 0;
    proModeFlowGrowthWindowCount = 0;
  }
  proModeFlowRollbackTotal += flowRollbackHit ? 1 : 0;
  proModeHeatRollbackTotal += heatRollbackHit ? 1 : 0;
  proModeFlowMagnitudeGuardTotal += flowMagnitudeGuardHit ? 1 : 0;
  proModeHeatMagnitudeGuardTotal += heatMagnitudeGuardHit ? 1 : 0;
  proModeFlowNoImproveTotal += flowNoImproveHit ? 1 : 0;
  proModeHeatNoImproveTotal += heatNoImproveHit ? 1 : 0;
  proModeFlowGrowthTotal += flowGrowthHit ? 1 : 0;

  proModeFlowRollbackWindowCount += flowRollbackHit ? 1 : 0;
  proModeHeatRollbackWindowCount += heatRollbackHit ? 1 : 0;
  proModeFlowMagnitudeGuardWindowCount += flowMagnitudeGuardHit ? 1 : 0;
  proModeHeatMagnitudeGuardWindowCount += heatMagnitudeGuardHit ? 1 : 0;
  proModeFlowNoImproveWindowCount += flowNoImproveHit ? 1 : 0;
  proModeHeatNoImproveWindowCount += heatNoImproveHit ? 1 : 0;
  proModeFlowGrowthWindowCount += flowGrowthHit ? 1 : 0;

  if (flowRollbackHit || heatRollbackHit)
  {
    proModeAdjointRollbackFreeStreak = 0;
  }
  else
  {
    ++proModeAdjointRollbackFreeStreak;
  }
  proModeForceAdjointResetNextIter = proModeCatastrophicAdjointEvent;
}
else
{
  proModeAdjointRollbackFreeStreak = 0;
  proModeForceAdjointResetNextIter = false;
}
const label rollbackWindowStartIter =
  ((opt - 1)/proModeRollbackDiagnosticWindow)*proModeRollbackDiagnosticWindow + 1;
const label rollbackWindowSpan = opt - rollbackWindowStartIter + 1;

if (Pstream::master())
{
  const scalar clockTime = runTime.elapsedClockTime();
  const char* designHealth =
    (xNonFiniteCountLocal > 0 || xLocalMin < -1e-8 || xLocalMax > 1.0 + 1e-8 || xhLocalMin < -1e-8 || xhLocalMax > 1.0 + 1e-8)
    ? "WARN"
    : (grayVolumeFraction > 0.50 ? "OBS" : "OK");
  const char* constraintHealth =
    (powerViolationPositive > 0 || volumeViolationPositive > 0 || thermalEnvelopeViolation > 0)
    ? "WARN"
    : "OK";
  const char* sensitivityHealth =
    (fsensNonFiniteCountLocal > 0 || gsensPowerNonFiniteCountLocal > 0 || gsensVolNonFiniteCountLocal > 0)
    ? "WARN"
    : "OK";
  const scalar severeTransitionFluidThreshold =
    scalar(0.08);
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
  scalar rollbackPenalty = 0.0;
  if (adjointCycleAdjointRolledBack)
  {
    rollbackPenalty += 20.0;
  }
  if (optimizationCycleHeatAdjointRolledBack)
  {
    rollbackPenalty += 25.0;
  }
  if (adjointCycleAdjointFlowDiverged)
  {
    rollbackPenalty += 15.0;
  }
  if (flowMagnitudeGuardHit)
  {
    rollbackPenalty += 10.0;
  }
  if (heatMagnitudeGuardHit)
  {
    rollbackPenalty += 10.0;
  }
  rollbackPenalty = Foam::min(rollbackPenalty, scalar(45.0));
  const scalar constraintScale = Foam::max(optConstraintTol, scalar(1e-8));
  const scalar constraintStress =
    std::log10
    (
      scalar(1.0)
    + Foam::max(guardConstraintViolation, scalar(0.0))/constraintScale
    );
  const scalar constraintPenalty =
    Foam::min(scalar(25.0), scalar(5.0)*constraintStress);
  const scalar residualScale =
    Foam::max(linearSolveWarnResidual, scalar(1e-14));
  const scalar residualStress =
    std::log10
    (
      scalar(1.0)
    + Foam::max(maxLinearFinalResidual, scalar(0.0))/residualScale
    );
  scalar linearPenalty =
    Foam::min(scalar(15.0), scalar(6.0)*residualStress);
  if (!residualGrowthWithinLimit)
  {
    linearPenalty = Foam::min(linearPenalty + scalar(5.0), scalar(20.0));
  }
  const scalar clipPenalty =
    customTurbulenceActive
    ? Foam::min(scalar(15.0), scalar(80.0)*Foam::max(maxClipFraction, scalar(0.0)))
    : scalar(0.0);
  const scalar consistencyPenalty =
    (opt > 10)
    ? Foam::min
      (
        scalar(10.0),
        scalar(4.0)*Foam::max(gradientConsistencyRatio - scalar(1.0), scalar(0.0))
      )
    : scalar(0.0);
  scalar trustPenalty(0.0);
  if (trustRatioDefined)
  {
    const scalar trustDeviation =
      Foam::mag(trustRatio - scalar(1.0))/Foam::max(Foam::mag(trustRatio), scalar(1.0));
    trustPenalty = Foam::min(scalar(10.0), scalar(10.0)*trustDeviation);
  }
  const scalar healthPenaltyTotal =
      rollbackPenalty
    + constraintPenalty
    + linearPenalty
    + clipPenalty
    + consistencyPenalty
    + trustPenalty;
  const scalar optimizerHealthScore =
    Foam::max(scalar(0.0), scalar(100.0) - healthPenaltyTotal);
  const char* optimizerHealthBand =
    optimizerHealthScore >= 85.0 ? "good"
    : (optimizerHealthScore >= 70.0 ? "watch"
      : (optimizerHealthScore >= 50.0 ? "poor" : "critical"));
  std::ofstream historyFile(optimizationHistoryFile, std::ios::app);
  if (opt > 1 && ((opt - 1) % 25 == 0))
  {
    historyFile
      << std::left
      << std::setw(8) << "Iter"
      << std::setw(16) << "PDActLog"
      << std::setw(16) << "PDCon"
      << std::setw(16) << "PD0Sug"
      << std::setw(16) << "V"
      << std::setw(16) << "MeanT"
      << std::setw(16) << "MaxT"
      << '\n';
  }
  historyFile
    << std::right
    << std::setw(8) << opt
    << std::scientific << std::setprecision(8)
    << std::setw(16) << PDActLog
    << std::setw(16) << PowerDiss
    << std::setw(16) << PD0Sug
    << std::setw(16) << V
    << std::setw(16) << MeanT
    << std::setw(16) << maxT
    << '\n';

  std::ofstream healthFile(healthMonitorFile, std::ios::app);
  if (opt > 1 && ((opt - 1) % 25 == 0))
  {
    healthFile
      << std::left
      << std::setw(8) << "Iter"
      << std::setw(14) << "HealthScore"
      << std::setw(14) << "PenTotal"
      << std::setw(14) << "RollPen"
      << std::setw(14) << "LinPen"
      << std::setw(14) << "ConPen"
      << std::setw(14) << "ClipPen"
      << std::setw(14) << "ConsPen"
      << std::setw(14) << "TrustPen"
      << std::setw(10) << "FlowRB"
      << std::setw(10) << "HeatRB"
      << '\n';
  }
  healthFile
    << std::right
    << std::setw(8) << opt
    << std::fixed << std::setprecision(2)
    << std::setw(14) << optimizerHealthScore
    << std::scientific << std::setprecision(6)
    << std::setw(14) << healthPenaltyTotal
    << std::setw(14) << rollbackPenalty
    << std::setw(14) << linearPenalty
    << std::setw(14) << constraintPenalty
    << std::setw(14) << clipPenalty
    << std::setw(14) << consistencyPenalty
    << std::setw(14) << trustPenalty
    << std::setw(10) << static_cast<int>(adjointCycleAdjointRolledBack)
    << std::setw(10) << static_cast<int>(optimizationCycleHeatAdjointRolledBack)
    << '\n';
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
      << " trustRatioDefined=" << (trustRatioDefined ? 1 : 0)
      << " trustRatio=" << trustRatio
      << " predictedObjectiveReduction=" << predictedObjectiveReduction
      << " actualObjectiveReduction=" << actualObjectiveReduction
      << " gradientConsistencyRatio=" << gradientConsistencyRatio
      << " maxFluidClip=" << maxFluidClipFraction
      << " maxGrayClip=" << maxGrayClipFraction
      << " thermalTmaxKS=" << thermalTmaxConstraintMargin
      << " thermalTminKS=" << thermalTminConstraintMargin
      << '\n';
  }

  std::ofstream debugLogFile(debugOptimizerLogFile, std::ios::app);
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
    "Health score",
    "score="+fmtScalar(optimizerHealthScore),
    "band="+std::string(optimizerHealthBand),
    "rollback="+fmtScalar(rollbackPenalty),
    "other="+fmtScalar(healthPenaltyTotal - rollbackPenalty)
  );
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
  writeRow
  (
    "Model agree",
    "predObjRed="+fmtScalar(predictedObjectiveReduction),
    "actObjRed="+fmtScalar(actualObjectiveReduction),
    "trust="+fmtScalar(trustRatio),
    "defined="+std::string(trustRatioDefined ? "yes" : "no")
  );
  writeRow
  (
    "Grad quality",
    "cons="+fmtScalar(gradientConsistencyRatio),
    "obj="+fmtScalar(objectiveConsistencyRatio),
    "pow="+fmtScalar(powerConsistencyRatio),
    "vol="+fmtScalar(volumeConsistencyRatio)
  );
  writeRow
  (
    "Linear health",
    "maxFinal="+fmtScalar(maxLinearFinalResidual),
    "growth="+fmtScalar(linearResidualGrowthRatio),
    "growthCheck="+std::string(residualGrowthWithinLimit ? "yes" : "no"),
    "healthCheck="+std::string(healthChecksPassed ? "yes" : "no")
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
    "envelopeStop="+std::string(fieldEnvelopeStop ? "yes" : "no"),
    "clipStop="+std::string(clipSafetyStop ? "yes" : "no")
      + " hard=" + std::string(clipHardSafetyStop ? "yes" : "no")
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
    "mode="+std::string(static_cast<bool>(optAdaptiveMoveLimitEnable) ? "adaptive" : "fixed"),
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
    : (turbulenceProModeActive
      ? proModeFlowAdjointCorrectors
      : (adjointCycleFullAdjointVerification ? adjointCycleAdjointCorrectors : 40));
  const label effectiveAdjHeatCorrectors =
      optimizationCycleTurbulentOptimizationMode
    ? optimizationCycleHeatAdjointCorrectors
    : (turbulenceProModeActive ? proModeHeatAdjointCorrectors : 40);
  const char* optimizationModeLabel =
      optimizationCycleTurbulentOptimizationMode
    ? "optimizationCycle"
    : (turbulenceProModeActive ? "turbulencePro"
    : (adjointCycleFullAdjointVerification ? "adjointCycle-full"
      : (adjointCycleThermalAdjointOnly ? "adjointCycle-thermal" : "baseline")));
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
  if (turbulenceProModeActive && proModeAdjointStabilityPackActive)
  {
    writeRow
    (
      "adj rb events",
      "flowRB="+fmtCount(flowRollbackHit ? 1 : 0),
      "heatRB="+fmtCount(heatRollbackHit ? 1 : 0),
      "flowMagG="+fmtCount(flowMagnitudeGuardHit ? 1 : 0),
      "heatMagG="+fmtCount(heatMagnitudeGuardHit ? 1 : 0)
    );
    writeRow
    (
      "adj rb reasons",
      "flowNoImp="+fmtCount(flowNoImproveHit ? 1 : 0),
      "heatNoImp="+fmtCount(heatNoImproveHit ? 1 : 0),
      "flowGrowth="+fmtCount(flowGrowthHit ? 1 : 0),
      "flowDiv="+fmtCount(adjointCycleAdjointFlowDiverged ? 1 : 0)
    );
    writeRow
    (
      "adj rb action",
      "catastrophic="+std::string(proModeCatastrophicAdjointEvent ? "yes" : "no"),
      "forceResetNext="+std::string(proModeForceAdjointResetNextIter ? "yes" : "no")
    );
    writeRow
    (
      "adj rb window",
      "iter="+fmtCount(rollbackWindowStartIter)+"-"+fmtCount(opt),
      "flow="+fmtCount(proModeFlowRollbackWindowCount)+"/"+fmtCount(rollbackWindowSpan)
        + " (" + fmtScalar(static_cast<scalar>(proModeFlowRollbackWindowCount)/Foam::max(static_cast<scalar>(rollbackWindowSpan), scalar(1.0))) + ")",
      "heat="+fmtCount(proModeHeatRollbackWindowCount)+"/"+fmtCount(rollbackWindowSpan)
        + " (" + fmtScalar(static_cast<scalar>(proModeHeatRollbackWindowCount)/Foam::max(static_cast<scalar>(rollbackWindowSpan), scalar(1.0))) + ")",
      "freeStreak="+fmtCount(proModeAdjointRollbackFreeStreak)
    );
    writeRow
    (
      "adj rb total",
      "flowRB="+fmtCount(proModeFlowRollbackTotal),
      "heatRB="+fmtCount(proModeHeatRollbackTotal),
      "flowMagG="+fmtCount(proModeFlowMagnitudeGuardTotal),
      "heatMagG="+fmtCount(proModeHeatMagnitudeGuardTotal)
    );
  }

  writeSectionHeader("Interpolation Diagnostics");
  writeRow("interp bounds", "alphaMin="+fmtScalar(alphaLocalMin), "alphaMax="+fmtScalar(alphaLocalMax), "DTMin="+fmtScalar(DTLocalMin), "DTMax="+fmtScalar(DTLocalMax));
  writeRow
  (
    "continuation",
    "alphaMax="+fmtScalar(alphaMax.value()),
    "qu="+fmtScalar(qu),
    "schedIter="+fmtCount(continuationScheduleIterationActive)
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
      << ",\"score\":" << optimizerHealthScore
      << ",\"band\":\"" << optimizerHealthBand << "\""
      << ",\"rollbackPenalty\":" << rollbackPenalty
      << ",\"constraintPenalty\":" << constraintPenalty
      << ",\"linearPenalty\":" << linearPenalty
      << ",\"clipPenalty\":" << clipPenalty
      << ",\"consistencyPenalty\":" << consistencyPenalty
      << ",\"trustPenalty\":" << trustPenalty
    << "}"
    << ",\"convergence\":{"
      << "\"objectiveRelativeDelta\":" << objectiveRelativeDelta
      << ",\"maxConstraintViolation\":" << maxConstraintViolation
      << ",\"maxConstraintViolationDelta\":" << maxConstraintViolationDelta
      << ",\"convergenceConstraintViolation\":" << convergenceConstraintViolation
      << ",\"guardConstraintViolation\":" << guardConstraintViolation
      << ",\"guardConstraintViolationDelta\":" << guardConstraintViolationDelta
      << ",\"powerConstraintFinalMargin\":" << powerConstraintFinalMargin
      << ",\"powerConstraintFinalMarginDelta\":" << powerConstraintFinalMarginDelta
      << ",\"volumeConstraintMarginDelta\":" << volumeConstraintMarginDelta
      << ",\"objectiveAbsValue\":" << objectiveAbsValue
      << ",\"objectiveAbsMax\":" << optObjectiveAbsMax
      << ",\"objectiveEnvelopeViolated\":" << (objectiveEnvelopeViolated ? 1 : 0)
      << ",\"thermalFieldAbsMaxObserved\":" << thermalFieldAbsMaxObserved
      << ",\"thermalFieldAbsMax\":" << optThermalFieldAbsMax
      << ",\"thermalFieldEnvelopeViolated\":" << (thermalFieldEnvelopeViolated ? 1 : 0)
      << ",\"fieldEnvelopeStop\":" << (fieldEnvelopeStop ? 1 : 0)
      << ",\"maxLinearFinalResidual\":" << maxLinearFinalResidual
      << ",\"linearResidualGrowthRatio\":" << linearResidualGrowthRatio
      << ",\"trustRatioDefined\":" << (trustRatioDefined ? 1 : 0)
      << ",\"trustRatio\":" << trustRatio
      << ",\"continuationScheduleIteration\":" << continuationScheduleIterationActive
      << ",\"closureDerivativeFDEnabled\":" << (static_cast<bool>(optClosureDerivativeFDEnable) ? 1 : 0)
      << ",\"closureDerivativeFDRelStep\":" << optClosureDerivativeFDRelStep
      << ",\"gradientConsistencyRatio\":" << gradientConsistencyRatio
      << ",\"gradientObjectiveRatio\":" << objectiveConsistencyRatio
      << ",\"gradientPowerRatio\":" << powerConsistencyRatio
      << ",\"gradientVolumeRatio\":" << volumeConsistencyRatio
      << ",\"predictedObjectiveDelta\":" << predictedObjectiveDelta
      << ",\"predictedPowerDelta\":" << predictedPowerDelta
      << ",\"predictedVolumeDelta\":" << predictedVolumeDelta
      << ",\"predictedThermalTmaxDelta\":" << predictedThermalTmaxDelta
      << ",\"predictedThermalTminDelta\":" << predictedThermalTminDelta
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
      << ",\"moveLimitMode\":\""
      << (static_cast<bool>(optAdaptiveMoveLimitEnable) ? "adaptive" : "fixed")
      << "\""
      << ",\"asymInc\":" << mma_inc
      << ",\"asymDec\":" << mma_dec
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
      << ",\"heatSource\":" << heatSource.value()
      << ",\"thermalSource\":" << thermalSource.value()
      << ",\"rhoCp\":" << rhoCp.value()
    << "}"
    << "}\n";
}
