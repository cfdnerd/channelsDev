Info<< "Reading optimization parameters\n" << endl;
IOdictionary optProperties
(
   IOobject
   (
       "optProperties",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE
   )
);
IOdictionary tuneOptParameters
(
   IOobject
   (
       "tuneOptParameters",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE
   )
);
//optimization parameters
scalar raa0(readScalar(optProperties.lookup("raa0")));
scalar mma_init(readScalar(optProperties.lookup("mma_init")));
scalar mma_dec(readScalar(optProperties.lookup("mma_dec")));
scalar mma_inc(readScalar(optProperties.lookup("mma_inc")));
scalar movlim(readScalar(optProperties.lookup("movlim")));
scalar voluse(readScalar(optProperties.lookup("voluse")));
scalar GeoDim(readScalar(optProperties.lookup("GeoDim")));
scalar PowerDiss0(readScalar(optProperties.lookup("PowerDiss0")));
scalar PowerDissMax(readScalar(optProperties.lookup("PowerDissMax")));
scalar PowerDissRelax(readScalar(optProperties.lookup("PowerDissRelax")));
scalar filterR(readScalar(optProperties.lookup("filterR")));
dimensionedScalar b("b", dimensionSet(0,-2,0,0,0,0,0),1.0);//filter radius of the PDE filter
dimensionedScalar b1("b1", dimensionSet(0,2,-2,-2,0,0,0),1.0);
Switch solid_area(optProperties.lookup("solid_area"));
Switch fluid_area(optProperties.lookup("fluid_area"));
Switch test_area(optProperties.lookup("test_area"));
scalar Pnorm(readScalar(optProperties.lookup("Pnorm")));
scalar objFunction(readScalar(optProperties.lookup("objFunction")));
const dictionary& convergenceControl = tuneOptParameters.subDict("convergenceControl");
label optMinIterations(convergenceControl.lookupOrDefault<label>("optMinIterations", 50));
label optStallWindow(convergenceControl.lookupOrDefault<label>("optStallWindow", 5));
scalar optObjectiveTol(convergenceControl.lookupOrDefault<scalar>("optObjectiveTol", 1e-3));
scalar optConstraintTol(convergenceControl.lookupOrDefault<scalar>("optConstraintTol", 1e-3));
scalar optStepTol(convergenceControl.lookupOrDefault<scalar>("optStepTol", 5e-3));
scalar optPowerConstraintWeight
(
  convergenceControl.lookupOrDefault<scalar>("optPowerConstraintWeight", 1.0)
);
scalar optVolumeConstraintWeight
(
  convergenceControl.lookupOrDefault<scalar>("optVolumeConstraintWeight", 1.0)
);
label optConstraintWorsenMinIterations
(
  convergenceControl.lookupOrDefault<label>("optConstraintWorsenMinIterations", optMinIterations)
);
label optConstraintWorsenWindow
(
  convergenceControl.lookupOrDefault<label>("optConstraintWorsenWindow", 0)
);
scalar optConstraintWorsenTol
(
  convergenceControl.lookupOrDefault<scalar>("optConstraintWorsenTol", 0.0)
);
scalar optConstraintWorsenNetRiseTol
(
  convergenceControl.lookupOrDefault<scalar>("optConstraintWorsenNetRiseTol", 0.25)
);
label optConstraintWorsenNetRiseMinAge
(
  convergenceControl.lookupOrDefault<label>("optConstraintWorsenNetRiseMinAge", optConstraintWorsenWindow)
);
label optConstraintReboundWindow
(
  convergenceControl.lookupOrDefault<label>("optConstraintReboundWindow", 0)
);
label optConstraintReboundMinAge
(
  convergenceControl.lookupOrDefault<label>("optConstraintReboundMinAge", optConstraintWorsenNetRiseMinAge)
);
scalar optConstraintReboundRiseTol
(
  convergenceControl.lookupOrDefault<scalar>("optConstraintReboundRiseTol", optConstraintWorsenNetRiseTol)
);
scalar optObjectiveWorsenTol
(
  convergenceControl.lookupOrDefault<scalar>("optObjectiveWorsenTol", 0.0)
);
label optObjectiveWorsenWindow
(
  convergenceControl.lookupOrDefault<label>("optObjectiveWorsenWindow", 0)
);
Switch optRecoveryEnable
(
  convergenceControl.lookupOrDefault<Switch>("optRecoveryEnable", false)
);
label optRecoveryWorsenStreak
(
  convergenceControl.lookupOrDefault<label>("optRecoveryWorsenStreak", 3)
);
scalar optRecoveryMoveLimitScale
(
  convergenceControl.lookupOrDefault<scalar>("optRecoveryMoveLimitScale", 0.5)
);
scalar optRecoveryContinuationRampScale
(
  convergenceControl.lookupOrDefault<scalar>("optRecoveryContinuationRampScale", 0.35)
);
scalar optObjectiveAbsMax
(
  convergenceControl.lookupOrDefault<scalar>("optObjectiveAbsMax", 1e6)
);
scalar optThermalFieldAbsMax
(
  convergenceControl.lookupOrDefault<scalar>("optThermalFieldAbsMax", 1e6)
);
Switch optFreezeContinuationOnReject(false);
Switch optFreezeContinuationOnLowTrust(false);
scalar optContinuationLowTrustThreshold(0.0);
label optContinuationLowTrustRejectStreak(0);
(void)optContinuationLowTrustRejectStreak;
Switch optRejectStopRequireConstraintWindow(false);
Switch optEarlyRecoveryOnGlobalizationReject(false);
// Retired dictionary tuning family: continuation deadlock/catch-up heuristics.
Switch optContinuationDeadlockRecoveryEnable(false);
Switch optContinuationClipHoldEnable(false);
scalar optContinuationClipHoldFluidFraction(0.20);
scalar optContinuationClipHoldNoProgressTol(5e-4);
label optContinuationDeadlockRecoveryMaxOverrides(0);
label optContinuationDeadlockRecoveryMaxCatchUpSteps(0);
label optContinuationDeadlockRecoveryLagBuffer(0);
// Retired dictionary tuning family: rogue-direction switching heuristics.
Switch optRogueDirectionSwitchEnable(false);
label optRogueGraceIterations(6);
label optRogueRejectStreakTrigger(3);
scalar optRogueClipFluidThreshold(0.20);
scalar optRogueClipGrayThreshold(0.90);
label optRogueDirectionSwitchCooldown(8);
label optRogueDirectionSwitchMaxCount(0);
scalar optRogueDirectionSwitchStepScale(0.6);
scalar optRogueDirectionSwitchConstraintBlend(0.7);
scalar optRogueDirectionSwitchNoiseScale(0.05);
scalar optRogueDirectionSwitchClipDampFluidFraction(0.60);
scalar optRogueDirectionSwitchClipDampGrayFraction(0.95);
scalar optRogueDirectionSwitchClipDampScale(0.50);
scalar optRogueDirectionSwitchSevereClipFluidFraction(0.80);
scalar optRogueDirectionSwitchSevereClipGrayFraction(0.98);
scalar optRogueDirectionSwitchSevereClipScale(0.25);
Switch optClosureDerivativeFDEnable
(
  convergenceControl.lookupOrDefault<Switch>("optClosureDerivativeFDEnable", true)
);
scalar optClosureDerivativeFDRelStep
(
  convergenceControl.lookupOrDefault<scalar>("optClosureDerivativeFDRelStep", 1e-4)
);
Switch optConstraintScheduleHorizonGuard
(
  convergenceControl.lookupOrDefault<Switch>("optConstraintScheduleHorizonGuard", true)
);
label optConstraintScheduleHorizonLagTolerance
(
  convergenceControl.lookupOrDefault<label>("optConstraintScheduleHorizonLagTolerance", 1)
);
Switch optRequireTransportedClosureAdjoints
(
  convergenceControl.lookupOrDefault<Switch>("optRequireTransportedClosureAdjoints", true)
);
Switch optRequireDualConsistentClosureBC
(
  convergenceControl.lookupOrDefault<Switch>("optRequireDualConsistentClosureBC", true)
);
word optRunMode
(
  convergenceControl.lookupOrDefault<word>("optRunMode", "production")
);
const bool optDiagnosticMode = optRunMode == word("diagnostic");
Switch optFirstFailureSnapshotEnable
(
  convergenceControl.lookupOrDefault<Switch>("optFirstFailureSnapshotEnable", true)
);
Switch optAdjointRollbackGuardsEnable
(
  convergenceControl.lookupOrDefault<Switch>("optAdjointRollbackGuardsEnable", true)
);
Switch optContinuationEventDriven
(
  convergenceControl.lookupOrDefault<Switch>("optContinuationEventDriven", true)
);
Switch optDisableLateStageContinuationAcceleration
(
  convergenceControl.lookupOrDefault<Switch>("optDisableLateStageContinuationAcceleration", true)
);
Switch optThermalBoundedObjectiveEnable
(
  convergenceControl.lookupOrDefault<Switch>("optThermalBoundedObjectiveEnable", true)
);
word optThermalAdjointMode
(
  convergenceControl.lookupOrDefault<word>("optThermalAdjointMode", "ta")
);
if
(
    optThermalAdjointMode != word("ta")
 && optThermalAdjointMode != word("tbOriginal")
 && optThermalAdjointMode != word("tbLaminar")
)
{
  FatalErrorInFunction
    << "Unsupported tuneOptParameters.convergenceControl.optThermalAdjointMode='"
    << optThermalAdjointMode
    << "'. Supported values: ta, tbOriginal, tbLaminar."
    << exit(FatalError);
}
scalar optThermalExcessTref
(
  convergenceControl.lookupOrDefault<scalar>("optThermalExcessTref", 300.0)
);
scalar optThermalExcessSoftplusScale
(
  convergenceControl.lookupOrDefault<scalar>("optThermalExcessSoftplusScale", 5.0)
);
Switch optThermalEnvelopeConstraintEnable
(
  convergenceControl.lookupOrDefault<Switch>("optThermalEnvelopeConstraintEnable", true)
);
Switch optThermalEnvelopeMMAConstraintEnable
(
  convergenceControl.lookupOrDefault<Switch>("optThermalEnvelopeMMAConstraintEnable", false)
);
scalar optThermalTmaxLimit
(
  convergenceControl.lookupOrDefault<scalar>("optThermalTmaxLimit", 450.0)
);
scalar optThermalTminLimit
(
  convergenceControl.lookupOrDefault<scalar>("optThermalTminLimit", 250.0)
);
scalar optThermalKSBeta
(
  convergenceControl.lookupOrDefault<scalar>("optThermalKSBeta", 24.0)
);
scalar optThermalEnvelopePenalty
(
  convergenceControl.lookupOrDefault<scalar>("optThermalEnvelopePenalty", 0.10)
);
Switch optProductionContinueAfterSafetyReject
(
  convergenceControl.lookupOrDefault<Switch>("optProductionContinueAfterSafetyReject", true)
);
scalar optSensitivityGradientL2Eps
(
  convergenceControl.lookupOrDefault<scalar>("optSensitivityGradientL2Eps", 1e-12)
);
Switch optDualConsistentStateMapEnforce
(
  convergenceControl.lookupOrDefault<Switch>("optDualConsistentStateMapEnforce", true)
);

if (optRunMode != word("production") && optRunMode != word("diagnostic"))
{
  FatalErrorInFunction
    << "tuneOptParameters.convergenceControl.optRunMode must be 'production' or 'diagnostic'. "
    << "Found " << optRunMode
    << exit(FatalError);
}
if (optMinIterations < 1 || optStallWindow < 1)
{
  FatalErrorInFunction
    << "tuneOptParameters.convergenceControl.optMinIterations and optStallWindow must be >= 1. "
    << "Found optMinIterations=" << optMinIterations
    << ", optStallWindow=" << optStallWindow
    << exit(FatalError);
}
if (optObjectiveTol <= 0.0 || optConstraintTol <= 0.0 || optStepTol <= 0.0)
{
  FatalErrorInFunction
    << "tuneOptParameters.convergenceControl.optObjectiveTol, optConstraintTol, and optStepTol "
    << "must be strictly positive. Found optObjectiveTol=" << optObjectiveTol
    << ", optConstraintTol=" << optConstraintTol
    << ", optStepTol=" << optStepTol
    << exit(FatalError);
}
if (optPowerConstraintWeight <= 0.0)
{
  FatalErrorInFunction
    << "tuneOptParameters.convergenceControl.optPowerConstraintWeight must be > 0. Found "
    << optPowerConstraintWeight
    << exit(FatalError);
}
if (optVolumeConstraintWeight <= 0.0)
{
  FatalErrorInFunction
    << "tuneOptParameters.convergenceControl.optVolumeConstraintWeight must be > 0. Found "
    << optVolumeConstraintWeight
    << exit(FatalError);
}
if (optThermalExcessSoftplusScale <= 0.0 || optThermalKSBeta <= 0.0)
{
  FatalErrorInFunction
    << "tuneOptParameters thermal smoothing parameters must be > 0. Found "
    << "optThermalExcessSoftplusScale=" << optThermalExcessSoftplusScale
    << ", optThermalKSBeta=" << optThermalKSBeta
    << exit(FatalError);
}
if (optThermalTmaxLimit <= optThermalTminLimit)
{
  FatalErrorInFunction
    << "tuneOptParameters thermal limits must satisfy optThermalTmaxLimit > optThermalTminLimit. "
    << "Found optThermalTmaxLimit=" << optThermalTmaxLimit
    << ", optThermalTminLimit=" << optThermalTminLimit
    << exit(FatalError);
}
if (optThermalEnvelopePenalty < 0.0)
{
  FatalErrorInFunction
    << "tuneOptParameters.convergenceControl.optThermalEnvelopePenalty must be >= 0. Found "
    << optThermalEnvelopePenalty
    << exit(FatalError);
}
if (optSensitivityGradientL2Eps <= 0.0)
{
  FatalErrorInFunction
    << "optSensitivityGradientL2Eps must be > 0. Found "
    << optSensitivityGradientL2Eps
    << exit(FatalError);
}

Info<< "Reading pseudo-density field\n" << endl;
volScalarField xp
(
    IOobject
    (
        "xp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField xh
(
    IOobject
    (
        "xh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField x(xp);
volScalarField drho(x);

// find the cells in the non-design area
labelList cells_solid,cells_fluid,cells_test;
if(solid_area)
{
  word zoneName="solid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_solid=mesh.cellZones()[zoneID];
  forAll(cells_solid, i)
  {
     x[cells_solid[i]]=0;
  }
}
if(fluid_area)
{
  word zoneName="fluid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_fluid=mesh.cellZones()[zoneID];
  forAll(cells_fluid, i)
  {
     x[cells_fluid[i]]=1.0;
  }
}
if(test_area)
{
  word zoneName="test_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_test=mesh.cellZones()[zoneID];
  forAll(cells_test, i)
  {
     x[cells_test[i]]=1.0;
     xh[cells_test[i]]=1.0;
  }
}
volScalarField gsensVol
(
    IOobject
    (
        "gsensVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    xh,
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshVol(xh);

// find the cells on the boundary
dictionary conFunc = mesh.solutionDict().subDict("constraintFunctionDict_flow");
int nObjPatch =conFunc.lookupOrDefault<scalar>("numberConstraintPatches",2);
wordList conPatchNames = conFunc.lookup("constraintPatchesNames");
label conPatchList[nObjPatch];//patch list
int iLoop;
for(iLoop=0; iLoop<nObjPatch; iLoop++)
{
    conPatchList[iLoop]=mesh.boundaryMesh().findPatchID(conPatchNames[iLoop]);
}   
