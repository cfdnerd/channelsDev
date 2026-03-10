Info<<"cost function"<<endl;
MeanT=fvc::domainIntegrate(T).value()/area;
maxT=gMax(T);
thermalExcessObjective = MeanT;
thermalTmaxKSConstraint = 0.0;
thermalTminKSConstraint = 0.0;

const scalar thermalScaleValue =
  Foam::max(thermalObjectiveSoftplusScale.value(), SMALL);
const scalar thermalRefValue = thermalObjectiveTref.value();
const scalar thermalTmaxLimitValue = thermalTmaxLimit.value();
const scalar thermalTminLimitValue = thermalTminLimit.value();

if
(
    !turbulenceProModeActive
 && (
      (objFunction == 3 && static_cast<bool>(optThermalBoundedObjectiveEnable))
   || static_cast<bool>(optThermalEnvelopeConstraintEnable)
    )
)
{
  scalar localExcessSum = 0.0;
  scalar localThetaMaxPeak = -GREAT;
  scalar localThetaMinPeak = -GREAT;
  for (i = 0; i < n; i++)
  {
    const scalar Ti = T[i];
    const scalar cellV = mesh.V()[i];
    const scalar z = (Ti - thermalRefValue)/thermalScaleValue;
    scalar softplusValue = 0.0;
    if (z > 40.0)
    {
      softplusValue = z;
    }
    else if (z < -40.0)
    {
      softplusValue = std::exp(z);
    }
    else
    {
      softplusValue = std::log(1.0 + std::exp(z));
    }
    localExcessSum += thermalScaleValue*softplusValue*cellV;

    const scalar thetaMax = (Ti - thermalTmaxLimitValue)/thermalScaleValue;
    const scalar thetaMin = (thermalTminLimitValue - Ti)/thermalScaleValue;
    localThetaMaxPeak = Foam::max(localThetaMaxPeak, thetaMax);
    localThetaMinPeak = Foam::max(localThetaMinPeak, thetaMin);
  }
  reduce(localExcessSum, sumOp<scalar>());
  reduce(localThetaMaxPeak, maxOp<scalar>());
  reduce(localThetaMinPeak, maxOp<scalar>());

  thermalExcessObjective = localExcessSum/Foam::max(area, SMALL);

  scalar localKsMaxExpSum = 0.0;
  scalar localKsMinExpSum = 0.0;
  for (i = 0; i < n; i++)
  {
    const scalar Ti = T[i];
    const scalar cellV = mesh.V()[i];
    const scalar thetaMax = (Ti - thermalTmaxLimitValue)/thermalScaleValue;
    const scalar thetaMin = (thermalTminLimitValue - Ti)/thermalScaleValue;
    localKsMaxExpSum += std::exp(thermalKSBeta*(thetaMax - localThetaMaxPeak))*cellV;
    localKsMinExpSum += std::exp(thermalKSBeta*(thetaMin - localThetaMinPeak))*cellV;
  }
  reduce(localKsMaxExpSum, sumOp<scalar>());
  reduce(localKsMinExpSum, sumOp<scalar>());

  const scalar ksMaxLogMean =
    std::log
    (
      Foam::max(localKsMaxExpSum/Foam::max(area, SMALL), SMALL)
    );
  const scalar ksMinLogMean =
    std::log
    (
      Foam::max(localKsMinExpSum/Foam::max(area, SMALL), SMALL)
    );
  thermalTmaxKSConstraint =
    localThetaMaxPeak + ksMaxLogMean/Foam::max(thermalKSBeta, SMALL);
  thermalTminKSConstraint =
    localThetaMinPeak + ksMinLogMean/Foam::max(thermalKSBeta, SMALL);
}
V=fvc::domainIntegrate(xh).value()/area-(voluse+0*Foam::max(0,0.5-opt*0.015));  
grayPenaltyObjective = 0.0;
if (static_cast<bool>(optGrayPenaltyEnable) && grayPenaltyWeightActive > SMALL)
{
  grayPenaltyObjective =
    grayPenaltyWeightActive
   *fvc::domainIntegrate(xh*(scalar(1.0) - xh)).value()
   /Foam::max(area, SMALL);
}

PowerDiss=0;
for (i=0; i<nObjPatch; i++)
{
     PowerDiss =PowerDiss - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
}
PowerDiss=PowerDiss/PowerDiss0;
if(GeoDim==2)
{
  PowerDiss=PowerDiss/std::pow(mesh.V()[0],1.0/3);
}
reduce(PowerDiss, sumOp<scalar>());

const scalar meshFactor = std::pow(mesh.V()[0],1.0/3);
PowerDissActualOverBase = PowerDiss;
if(GeoDim==2)
{
  PowerDissActualOverBase = PowerDissActualOverBase*meshFactor;
}
PowerDissActual = PowerDissActualOverBase*PowerDiss0;
const scalar PDAct = PowerDiss*PowerDiss0*meshFactor;
PDActLog = PDAct/meshFactor;
PD0Sug = PowerDiss0;
if (PowerDiss > PowerDissRelax && PowerDissMax > SMALL)
{
  PD0Sug = PowerDiss0*(PowerDiss/Foam::max(PowerDissMax, SMALL));
}
const bool pd0SugDiffObserved = Foam::mag(PD0Sug - PowerDiss0) > SMALL;
if (pd0SugDiffObserved)
{
  PD0SugMaxCum = Foam::max(PD0SugMaxCum, PD0Sug);
  PD0SugSumCum += PD0Sug;
  PD0SugCountCum++;
}

if (Pstream::master())
{ 
   if (pd0SugDiffObserved && !PD0SugLogInitialized)
   {
     std::ofstream PD0SugLog(suggestedPowerDiss0File, std::ios::out);
     PD0SugLog
       << "# MTO_ThermalFluid suggested PowerDiss0 summary\n"
       << "# Cols: Iter PD0SugMax_1_to_Iter PD0SugAvg_1_to_Iter\n"
       << "# Stats exclude unchanged values where PD0Sug == PowerDiss0\n"
       << std::left
       << std::setw(8) << "Iter"
       << std::setw(20) << "PD0SugMax_1_to_Iter"
       << std::setw(20) << "PD0SugAvg_1_to_Iter"
       << '\n';
     PD0SugLogInitialized = true;
   }
   if (PD0SugLogInitialized && PD0SugCountCum > 0 && opt > 0 && (opt % 100 == 0))
   {
     const scalar PD0SugAvgCum =
       PD0SugSumCum/Foam::max(static_cast<scalar>(PD0SugCountCum), SMALL);
     std::ofstream PD0SugLog(suggestedPowerDiss0File, std::ios::app);
     PD0SugLog
       << std::right
       << std::setw(8) << opt
       << std::scientific << std::setprecision(8)
       << std::setw(20) << PD0SugMaxCum
       << std::setw(20) << PD0SugAvgCum
       << '\n';
   }
}
