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

if ((objFunction == 3 && static_cast<bool>(optThermalBoundedObjectiveEnable))
 || static_cast<bool>(optThermalEnvelopeConstraintEnable))
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

if (Pstream::master())
{ 
   std::ofstream historyFile(optimizationHistoryFile, std::ios::app);
   if (opt > 1 && ((opt - 1) % 25 == 0))
   {
     historyFile
       << std::left
       << std::setw(8) << "Iter"
       << std::setw(16) << "PDActLog"
       << std::setw(16) << "PDCon"
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
     << std::setw(16) << V
     << std::setw(16) << MeanT
     << std::setw(16) << maxT
     << '\n';
}
