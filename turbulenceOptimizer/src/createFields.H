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
scalar optObjectiveAbsMax
(
  convergenceControl.lookupOrDefault<scalar>("optObjectiveAbsMax", 1e6)
);
scalar optThermalFieldAbsMax
(
  convergenceControl.lookupOrDefault<scalar>("optThermalFieldAbsMax", 1e6)
);
Switch optClosureDerivativeFDEnable
(
  convergenceControl.lookupOrDefault<Switch>("optClosureDerivativeFDEnable", true)
);
scalar optClosureDerivativeFDRelStep
(
  convergenceControl.lookupOrDefault<scalar>("optClosureDerivativeFDRelStep", 1e-4)
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
 && optThermalAdjointMode != word("tb")
)
{
  FatalErrorInFunction
    << "Unsupported tuneOptParameters.convergenceControl.optThermalAdjointMode='"
    << optThermalAdjointMode
    << "'. Supported values: ta, tb."
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
const Switch optThermalEnvelopeConstraintsEnable
(
  convergenceControl.lookupOrDefault<Switch>
  (
    "optThermalEnvelopeConstraintsEnable",
    convergenceControl.lookupOrDefault<Switch>
    (
      "optThermalEnvelopeConstraintEnable",
      convergenceControl.lookupOrDefault<Switch>("optThermalEnvelopeMMAConstraintEnable", false)
    )
  )
);
Switch optThermalEnvelopeConstraintEnable(optThermalEnvelopeConstraintsEnable);
Switch optThermalEnvelopeMMAConstraintEnable(optThermalEnvelopeConstraintsEnable);
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
scalar optSensitivityGradientL2Eps
(
  convergenceControl.lookupOrDefault<scalar>("optSensitivityGradientL2Eps", 1e-12)
);
Switch optDualConsistentStateMapEnforce
(
  convergenceControl.lookupOrDefault<Switch>("optDualConsistentStateMapEnforce", true)
);
Switch optClosureLogRegularizationEnable
(
  convergenceControl.lookupOrDefault<Switch>("optClosureLogRegularizationEnable", true)
);
Switch optClosureSmoothDenominatorEnable
(
  convergenceControl.lookupOrDefault<Switch>("optClosureSmoothDenominatorEnable", true)
);
Switch optKOmegaPositivityPreservingLinearizationEnable
(
  convergenceControl.lookupOrDefault<Switch>
  (
    "optKOmegaPositivityPreservingLinearizationEnable",
    true
  )
);
Switch optDivisionAuditEnable
(
  convergenceControl.lookupOrDefault<Switch>("optDivisionAuditEnable", true)
);
Switch optProjectionSmoothingEnable
(
  convergenceControl.lookupOrDefault<Switch>("optProjectionSmoothingEnable", true)
);
scalar optProjectionBetaScale
(
  convergenceControl.lookupOrDefault<scalar>("optProjectionBetaScale", 1.0)
);
Switch optSensitivityVolumeWeightingEnable
(
  convergenceControl.lookupOrDefault<Switch>("optSensitivityVolumeWeightingEnable", true)
);
Switch optAdaptiveMoveLimitEnable
(
  convergenceControl.lookupOrDefault<Switch>("optAdaptiveMoveLimitEnable", true)
);
scalar optAdaptiveMoveLimitMin
(
  convergenceControl.lookupOrDefault<scalar>("optAdaptiveMoveLimitMin", 0.05)
);
scalar optAdaptiveMoveLimitMax
(
  convergenceControl.lookupOrDefault<scalar>("optAdaptiveMoveLimitMax", 0.30)
);
scalar optAdaptiveMoveLimitTrustRatioShrink
(
  convergenceControl.lookupOrDefault<scalar>("optAdaptiveMoveLimitTrustRatioShrink", 0.25)
);
scalar optAdaptiveMoveLimitTrustRatioGrow
(
  convergenceControl.lookupOrDefault<scalar>("optAdaptiveMoveLimitTrustRatioGrow", 0.75)
);
scalar optAdaptiveMoveLimitShrinkFactor
(
  convergenceControl.lookupOrDefault<scalar>("optAdaptiveMoveLimitShrinkFactor", 0.70)
);
scalar optAdaptiveMoveLimitGrowFactor
(
  convergenceControl.lookupOrDefault<scalar>("optAdaptiveMoveLimitGrowFactor", 1.15)
);
scalar optAdaptiveMoveLimitNegativeTrustDamping
(
  convergenceControl.lookupOrDefault<scalar>("optAdaptiveMoveLimitNegativeTrustDamping", 0.50)
);
Switch optWallFunctionYPlusDiagnosticsEnable
(
  convergenceControl.lookupOrDefault<Switch>("optWallFunctionYPlusDiagnosticsEnable", true)
);
scalar optWallFunctionYPlusMin
(
  convergenceControl.lookupOrDefault<scalar>("optWallFunctionYPlusMin", 30.0)
);
scalar optWallFunctionYPlusMax
(
  convergenceControl.lookupOrDefault<scalar>("optWallFunctionYPlusMax", 300.0)
);
scalar optWallFunctionYPlusInBandWarnFraction
(
  convergenceControl.lookupOrDefault<scalar>("optWallFunctionYPlusInBandWarnFraction", 0.50)
);
Switch optConvergedGrayTrustGateEnable
(
  convergenceControl.lookupOrDefault<Switch>("optConvergedGrayTrustGateEnable", true)
);
label optConvergedGrayTrustGateMinIterations
(
  convergenceControl.lookupOrDefault<label>
  (
    "optConvergedGrayTrustGateMinIterations",
    optMinIterations
  )
);
scalar optConvergedGrayFractionThreshold
(
  convergenceControl.lookupOrDefault<scalar>("optConvergedGrayFractionThreshold", 0.50)
);
scalar optConvergedTrustPenaltySaturationThreshold
(
  convergenceControl.lookupOrDefault<scalar>("optConvergedTrustPenaltySaturationThreshold", 9.90)
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
if (optProjectionBetaScale <= 0.0)
{
  FatalErrorInFunction
    << "optProjectionBetaScale must be > 0. Found "
    << optProjectionBetaScale
    << exit(FatalError);
}
if (optAdaptiveMoveLimitMin <= 0.0 || optAdaptiveMoveLimitMax <= 0.0)
{
  FatalErrorInFunction
    << "optAdaptiveMoveLimitMin and optAdaptiveMoveLimitMax must be > 0. Found "
    << "min=" << optAdaptiveMoveLimitMin
    << ", max=" << optAdaptiveMoveLimitMax
    << exit(FatalError);
}
if (optWallFunctionYPlusMin <= 0.0 || optWallFunctionYPlusMax <= optWallFunctionYPlusMin)
{
  FatalErrorInFunction
    << "Wall-function yPlus band must satisfy 0 < optWallFunctionYPlusMin < optWallFunctionYPlusMax. "
    << "Found min=" << optWallFunctionYPlusMin
    << ", max=" << optWallFunctionYPlusMax
    << exit(FatalError);
}
if
(
    optWallFunctionYPlusInBandWarnFraction < 0.0
 || optWallFunctionYPlusInBandWarnFraction > 1.0
)
{
  FatalErrorInFunction
    << "optWallFunctionYPlusInBandWarnFraction must be in [0, 1]. Found "
    << optWallFunctionYPlusInBandWarnFraction
    << exit(FatalError);
}
if (optConvergedGrayTrustGateMinIterations < 1)
{
  FatalErrorInFunction
    << "optConvergedGrayTrustGateMinIterations must be >= 1. Found "
    << optConvergedGrayTrustGateMinIterations
    << exit(FatalError);
}
if
(
    optConvergedGrayFractionThreshold < 0.0
 || optConvergedGrayFractionThreshold > 1.0
)
{
  FatalErrorInFunction
    << "optConvergedGrayFractionThreshold must be in [0, 1]. Found "
    << optConvergedGrayFractionThreshold
    << exit(FatalError);
}
if
(
    optConvergedTrustPenaltySaturationThreshold < 0.0
 || optConvergedTrustPenaltySaturationThreshold > 10.0
)
{
  FatalErrorInFunction
    << "optConvergedTrustPenaltySaturationThreshold must be in [0, 10]. Found "
    << optConvergedTrustPenaltySaturationThreshold
    << exit(FatalError);
}
if (optAdaptiveMoveLimitMax < optAdaptiveMoveLimitMin)
{
  FatalErrorInFunction
    << "optAdaptiveMoveLimitMax must be >= optAdaptiveMoveLimitMin. Found "
    << "min=" << optAdaptiveMoveLimitMin
    << ", max=" << optAdaptiveMoveLimitMax
    << exit(FatalError);
}
if
(
    optAdaptiveMoveLimitTrustRatioShrink < 0.0
 || optAdaptiveMoveLimitTrustRatioGrow < 0.0
 || optAdaptiveMoveLimitTrustRatioGrow < optAdaptiveMoveLimitTrustRatioShrink
)
{
  FatalErrorInFunction
    << "adaptive move-limit trust-ratio thresholds must satisfy 0 <= shrink <= grow. Found "
    << "shrink=" << optAdaptiveMoveLimitTrustRatioShrink
    << ", grow=" << optAdaptiveMoveLimitTrustRatioGrow
    << exit(FatalError);
}
if
(
    optAdaptiveMoveLimitShrinkFactor <= 0.0
 || optAdaptiveMoveLimitShrinkFactor > 1.0
 || optAdaptiveMoveLimitGrowFactor < 1.0
)
{
  FatalErrorInFunction
    << "adaptive move-limit factors must satisfy 0 < shrink <= 1 and grow >= 1. Found "
    << "shrinkFactor=" << optAdaptiveMoveLimitShrinkFactor
    << ", growFactor=" << optAdaptiveMoveLimitGrowFactor
    << exit(FatalError);
}
if (optAdaptiveMoveLimitNegativeTrustDamping <= 0.0 || optAdaptiveMoveLimitNegativeTrustDamping > 1.0)
{
  FatalErrorInFunction
    << "optAdaptiveMoveLimitNegativeTrustDamping must satisfy 0 < value <= 1. Found "
    << optAdaptiveMoveLimitNegativeTrustDamping
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
