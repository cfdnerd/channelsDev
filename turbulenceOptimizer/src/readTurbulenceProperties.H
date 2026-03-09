Info<< "Reading custom turbulence properties\n" << endl;
IOdictionary customTurbulenceProperties
(
    IOobject
    (
        "customTurbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Switch customTurbulenceActive(customTurbulenceProperties.lookup("enable"));
Switch turbulenceProMode
(
    customTurbulenceProperties.lookupOrDefault<Switch>("turbulenceProMode", false)
);
Switch proModeClosureAdjointSolves
(
    customTurbulenceProperties.lookupOrDefault<Switch>("proModeClosureAdjointSolves", true)
);
const bool turbulenceProModeActive =
    customTurbulenceActive
 && static_cast<bool>(turbulenceProMode);
const bool proModeClosureAdjointSolvesActive =
    turbulenceProModeActive
 && static_cast<bool>(proModeClosureAdjointSolves);
Switch adjointCycleThermalAdjointOnly
(
    customTurbulenceProperties.lookupOrDefault<Switch>("adjointCycleThermalAdjointOnly", false)
);
Switch adjointCycleFullAdjointVerification
(
    customTurbulenceProperties.lookupOrDefault<Switch>("adjointCycleFullAdjointVerification", false)
);
Switch optimizationCycleGradientClosure
(
    customTurbulenceProperties.lookupOrDefault<Switch>("optimizationCycleGradientClosure", false)
);
Switch enableTurbulentOptimization
(
    customTurbulenceProperties.lookupOrDefault<Switch>("enableTurbulentOptimization", false)
);
const bool optimizationCycleTurbulentOptimizationMode =
    customTurbulenceActive
 && !turbulenceProModeActive
 && optimizationCycleGradientClosure
 && enableTurbulentOptimization;
const bool turbulenceAdjointClosureMode =
    optimizationCycleTurbulentOptimizationMode
 || turbulenceProModeActive;
word inletTurbulenceMode
(
    customTurbulenceProperties.lookupOrDefault<word>("inletTurbulenceMode", "intensityLengthScale")
);
scalar inletVelocityMagnitude
(
    customTurbulenceProperties.lookupOrDefault<scalar>("Uin", 1.0)
);
scalar inletTurbulenceIntensity
(
    customTurbulenceProperties.lookupOrDefault<scalar>("turbulenceIntensity", 0.05)
);
dimensionedScalar inletTurbulenceLengthScale
(
    customTurbulenceProperties.lookupOrDefault<dimensionedScalar>
    (
        "turbulenceLengthScale",
        dimensionedScalar("turbulenceLengthScale", dimLength, 3.5e-4)
    )
);
scalar inletConsistencyRelTol
(
    customTurbulenceProperties.lookupOrDefault<scalar>("inletConsistencyRelTol", 0.20)
);
word customTurbulenceModelName(customTurbulenceProperties.lookup("model"));
dimensionedScalar alphaK1Coeff(customTurbulenceProperties.lookup("alphaK1"));
dimensionedScalar alphaK2Coeff(customTurbulenceProperties.lookup("alphaK2"));
dimensionedScalar alphaOmega1Coeff(customTurbulenceProperties.lookup("alphaOmega1"));
dimensionedScalar alphaOmega2Coeff(customTurbulenceProperties.lookup("alphaOmega2"));
dimensionedScalar gamma1Coeff(customTurbulenceProperties.lookup("gamma1"));
dimensionedScalar gamma2Coeff(customTurbulenceProperties.lookup("gamma2"));
dimensionedScalar beta1Coeff(customTurbulenceProperties.lookup("beta1"));
dimensionedScalar beta2Coeff(customTurbulenceProperties.lookup("beta2"));
dimensionedScalar betaStarCoeff(customTurbulenceProperties.lookup("betaStar"));
dimensionedScalar a1Coeff(customTurbulenceProperties.lookup("a1"));
dimensionedScalar b1Coeff(customTurbulenceProperties.lookup("b1"));
dimensionedScalar c1Coeff(customTurbulenceProperties.lookup("c1"));
dimensionedScalar Prt(customTurbulenceProperties.lookup("Prt"));
dimensionedScalar kMin(customTurbulenceProperties.lookup("kMin"));
dimensionedScalar omegaMin(customTurbulenceProperties.lookup("omegaMin"));
dimensionedScalar omegaMax(customTurbulenceProperties.lookup("omegaMax"));
scalar omegaWallFactor
(
    customTurbulenceProperties.lookupOrDefault<scalar>("omegaWallFactor", 60.0)
);
dimensionedScalar nutMax(customTurbulenceProperties.lookup("nutMax"));
dimensionedScalar yMin(customTurbulenceProperties.lookup("yMin"));
dimensionedScalar nutFloor(customTurbulenceProperties.lookup("nutFloor"));
scalar turbulenceDampingCoeff(readScalar(customTurbulenceProperties.lookup("turbulenceDampingCoeff")));
scalar turbulencePrimalTransportDiffusionScale
(
    customTurbulenceProperties.lookupOrDefault<scalar>
    (
        "turbulencePrimalTransportDiffusionScale",
        1.0
    )
);
scalar turbulenceAdjointTransportDiffusionScale
(
    customTurbulenceProperties.lookupOrDefault<scalar>
    (
        "turbulenceAdjointTransportDiffusionScale",
        1.0
    )
);
if (kMin.value() <= SMALL)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.kMin must be strictly positive. Found "
        << kMin.value()
        << exit(FatalError);
}
if (omegaMin.value() <= SMALL)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.omegaMin must be strictly positive. Found "
        << omegaMin.value()
        << exit(FatalError);
}
if (omegaMax.value() <= omegaMin.value())
{
    FatalErrorInFunction
        << "customTurbulenceProperties.omegaMax must be greater than omegaMin. "
        << "Found omegaMax=" << omegaMax.value()
        << ", omegaMin=" << omegaMin.value()
        << exit(FatalError);
}
if (nutFloor.value() <= SMALL)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.nutFloor must be strictly positive. Found "
        << nutFloor.value()
        << exit(FatalError);
}
if (nutMax.value() <= nutFloor.value())
{
    FatalErrorInFunction
        << "customTurbulenceProperties.nutMax must be greater than nutFloor. "
        << "Found nutMax=" << nutMax.value()
        << ", nutFloor=" << nutFloor.value()
        << exit(FatalError);
}
if (yMin.value() <= SMALL)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.yMin must be strictly positive. Found "
        << yMin.value()
        << exit(FatalError);
}
if (inletVelocityMagnitude <= 0.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.Uin must be strictly positive. Found "
        << inletVelocityMagnitude
        << exit(FatalError);
}
if (inletConsistencyRelTol < 0.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.inletConsistencyRelTol must be non-negative. Found "
        << inletConsistencyRelTol
        << exit(FatalError);
}
if (omegaWallFactor <= 0.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.omegaWallFactor must be strictly positive. Found "
        << omegaWallFactor
        << exit(FatalError);
}
if (turbulenceDampingCoeff < 0.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.turbulenceDampingCoeff must be non-negative. Found "
        << turbulenceDampingCoeff
        << exit(FatalError);
}
if (turbulencePrimalTransportDiffusionScale < 1.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.turbulencePrimalTransportDiffusionScale must be >= 1. "
        << "Found " << turbulencePrimalTransportDiffusionScale
        << exit(FatalError);
}
if (turbulenceAdjointTransportDiffusionScale < 1.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.turbulenceAdjointTransportDiffusionScale must be >= 1. "
        << "Found " << turbulenceAdjointTransportDiffusionScale
        << exit(FatalError);
}
Switch porousAwareTurbulenceTransport
(
    customTurbulenceProperties.lookupOrDefault<Switch>("porousAwareTurbulenceTransport", true)
);
scalar porousTurbulenceProductionExponent
(
    customTurbulenceProperties.lookupOrDefault<scalar>("porousTurbulenceProductionExponent", 2.0)
);
scalar porousTurbulenceProductionFloor
(
    customTurbulenceProperties.lookupOrDefault<scalar>("porousTurbulenceProductionFloor", 5e-3)
);
scalar porousTurbulenceSuppressionCoeff
(
    customTurbulenceProperties.lookupOrDefault<scalar>("porousTurbulenceSuppressionCoeff", 1.0)
);
scalar porousTurbulenceSuppressionExponent
(
    customTurbulenceProperties.lookupOrDefault<scalar>("porousTurbulenceSuppressionExponent", 4.0)
);
scalar porousTurbulenceAlphaSinkCoeff
(
    customTurbulenceProperties.lookupOrDefault<scalar>("porousTurbulenceAlphaSinkCoeff", 0.5)
);
label momentumCorrectors
(
    customTurbulenceProperties.lookupOrDefault<label>("momentumCorrectors", 60)
);
label turbulenceCorrectors
(
    customTurbulenceProperties.lookupOrDefault<label>("turbulenceCorrectors", 1)
);
label minTurbulenceCorrectors
(
    customTurbulenceProperties.lookupOrDefault<label>("minTurbulenceCorrectors", 2)
);
scalar turbulenceKNonlinearTol
(
    customTurbulenceProperties.lookupOrDefault<scalar>("kNonlinearTol", 1e-4)
);
scalar turbulenceOmegaNonlinearTol
(
    customTurbulenceProperties.lookupOrDefault<scalar>("omegaNonlinearTol", 1e-4)
);
scalar turbulenceNutNonlinearTol
(
    customTurbulenceProperties.lookupOrDefault<scalar>("nutNonlinearTol", 1e-4)
);
// Retired dictionary tuning knobs: fixed Adjoint-Cycle/Optimization-Cycle legacy fallback caps.
label adjointCycleAdjointCorrectors(8);
label optimizationCycleAdjointCorrectors(12);
label optimizationCycleFlowAdjointCorrectors
(
    customTurbulenceProperties.lookupOrDefault<label>
    (
        "optimizationCycleFlowAdjointCorrectors",
        optimizationCycleAdjointCorrectors
    )
);
label optimizationCycleHeatAdjointCorrectors
(
    customTurbulenceProperties.lookupOrDefault<label>
    (
        "optimizationCycleHeatAdjointCorrectors",
        optimizationCycleAdjointCorrectors
    )
);
Switch optimizationCycleTransportedClosureAdjoints
(
    customTurbulenceProperties.lookupOrDefault<Switch>("optimizationCycleTransportedClosureAdjoints", false)
);
word optimizationCycleClosureAdjointFormulation
(
    customTurbulenceProperties.lookupOrDefault<word>("optimizationCycleClosureAdjointFormulation", "proxyNormalized")
);
Switch optimizationCycleClosureAdjointDualConsistentBC
(
    customTurbulenceProperties.lookupOrDefault<Switch>("optimizationCycleClosureAdjointDualConsistentBC", false)
);
scalar optimizationCycleClosureAdjointInvTau
(
    customTurbulenceProperties.lookupOrDefault<scalar>("optimizationCycleClosureAdjointInvTau", 1.0)
);
scalar optimizationCycleClosureCrossCoupling
(
    customTurbulenceProperties.lookupOrDefault<scalar>("optimizationCycleClosureCrossCoupling", 0.35)
);
Switch optimizationCycleClosureConsistencyUseBoundedState
(
    customTurbulenceProperties.lookupOrDefault<Switch>("optimizationCycleClosureConsistencyUseBoundedState", true)
);
Switch optimizationCycleUseLogClosureState
(
    customTurbulenceProperties.lookupOrDefault<Switch>("optimizationCycleUseLogClosureState", true)
);
// Retired dictionary tuning knobs: fixed staged-adjoint guard policy.
scalar adjointCycleAdjointUNonlinearTol(5e-3);
scalar adjointCycleAdjointPNonlinearTol(5e-3);
scalar adjointCycleAdjointMaxMagnitude(1e8);
scalar optimizationCycleAdjointMaxMagnitude(1e10);
Switch adjointCycleResetAdjointState(true);
Switch optimizationCycleResetAdjointState(true);
Switch adjointCycleUseNuEff(true);
Switch adjointCycleIncludeThermalMomentumSource(true);
Switch adjointCycleUseLegacyPressureRelax(false);
Switch adjointCycleDivergenceCheck(true);
Switch adjointCycleStagnationCheck(true);
label adjointCycleDivergenceMinCorrectors
(
    Foam::min
    (
        label(5),
        adjointCycleAdjointCorrectors
    )
);
scalar adjointCycleDivergenceGrowthFactor(1.5);
label adjointCycleStagnationPatience(1);
Switch optimizationCycleNoImproveCheck(true);
label optimizationCycleNoImproveMinCorrectors
(
    Foam::min
    (
        label(2),
        Foam::max(optimizationCycleFlowAdjointCorrectors, optimizationCycleHeatAdjointCorrectors)
    )
);
label optimizationCycleNoImprovePatience(1);

if (customTurbulenceModelName != "kOmegaSST")
{
    FatalErrorInFunction
        << "customTurbulenceProperties.model must be kOmegaSST. Found "
        << customTurbulenceModelName
        << exit(FatalError);
}
if
(
    optimizationCycleClosureAdjointFormulation != "proxyNormalized"
 && optimizationCycleClosureAdjointFormulation != "transportedResidualDraft"
 && optimizationCycleClosureAdjointFormulation != "transportedCrossCoupledV3"
)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.optimizationCycleClosureAdjointFormulation must "
        << "be proxyNormalized, transportedResidualDraft, or "
        << "transportedCrossCoupledV3. Found "
        << optimizationCycleClosureAdjointFormulation
        << exit(FatalError);
}
if (porousTurbulenceProductionExponent < 1.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.porousTurbulenceProductionExponent must be >= 1. "
        << "Found " << porousTurbulenceProductionExponent
        << exit(FatalError);
}
if (porousTurbulenceProductionFloor < 0.0 || porousTurbulenceProductionFloor > 1.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.porousTurbulenceProductionFloor must be in [0, 1]. "
        << "Found " << porousTurbulenceProductionFloor
        << exit(FatalError);
}
if (porousTurbulenceSuppressionCoeff < 0.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.porousTurbulenceSuppressionCoeff must be >= 0. "
        << "Found " << porousTurbulenceSuppressionCoeff
        << exit(FatalError);
}
if (porousTurbulenceSuppressionExponent < 1.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.porousTurbulenceSuppressionExponent must be >= 1. "
        << "Found " << porousTurbulenceSuppressionExponent
        << exit(FatalError);
}
if (porousTurbulenceAlphaSinkCoeff < 0.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.porousTurbulenceAlphaSinkCoeff must be >= 0. "
        << "Found " << porousTurbulenceAlphaSinkCoeff
        << exit(FatalError);
}
if (momentumCorrectors < 1)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.momentumCorrectors must be >= 1. Found "
        << momentumCorrectors
        << exit(FatalError);
}
if (turbulenceCorrectors < 1)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.turbulenceCorrectors must be >= 1. Found "
        << turbulenceCorrectors
        << exit(FatalError);
}
if (minTurbulenceCorrectors < 1 || minTurbulenceCorrectors > turbulenceCorrectors)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.minTurbulenceCorrectors must satisfy "
        << "1 <= minTurbulenceCorrectors <= turbulenceCorrectors. Found "
        << "minTurbulenceCorrectors=" << minTurbulenceCorrectors
        << ", turbulenceCorrectors=" << turbulenceCorrectors
        << exit(FatalError);
}
if
(
    turbulenceKNonlinearTol <= 0.0
 || turbulenceOmegaNonlinearTol <= 0.0
 || turbulenceNutNonlinearTol <= 0.0
)
{
    FatalErrorInFunction
        << "kNonlinearTol, omegaNonlinearTol, and nutNonlinearTol must all be > 0. Found "
        << "kNonlinearTol=" << turbulenceKNonlinearTol
        << ", omegaNonlinearTol=" << turbulenceOmegaNonlinearTol
        << ", nutNonlinearTol=" << turbulenceNutNonlinearTol
        << exit(FatalError);
}
if (optimizationCycleFlowAdjointCorrectors < 1 || optimizationCycleHeatAdjointCorrectors < 1)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.optimizationCycleFlowAdjointCorrectors and "
        << "optimizationCycleHeatAdjointCorrectors must both be >= 1. Found flow="
        << optimizationCycleFlowAdjointCorrectors
        << ", heat=" << optimizationCycleHeatAdjointCorrectors
        << exit(FatalError);
}
if (optimizationCycleClosureAdjointInvTau <= SMALL)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.optimizationCycleClosureAdjointInvTau must be > 0. Found "
        << optimizationCycleClosureAdjointInvTau
        << exit(FatalError);
}
if (optimizationCycleClosureCrossCoupling < 0.0 || optimizationCycleClosureCrossCoupling > 1.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.optimizationCycleClosureCrossCoupling must be in [0, 1]. Found "
        << optimizationCycleClosureCrossCoupling
        << exit(FatalError);
}
if
(
    (adjointCycleThermalAdjointOnly && adjointCycleFullAdjointVerification)
 || (adjointCycleThermalAdjointOnly && optimizationCycleTurbulentOptimizationMode)
 || (adjointCycleFullAdjointVerification && optimizationCycleTurbulentOptimizationMode)
 || (adjointCycleThermalAdjointOnly && turbulenceProModeActive)
 || (adjointCycleFullAdjointVerification && turbulenceProModeActive)
 || (optimizationCycleTurbulentOptimizationMode && turbulenceProModeActive)
)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.adjointCycleThermalAdjointOnly, "
        << "customTurbulenceProperties.adjointCycleFullAdjointVerification, and "
        << "the Optimization-Cycle turbulent optimization mode, and "
        << "customTurbulenceProperties.turbulenceProMode are "
        << "mutually exclusive run modes."
        << exit(FatalError);
}
if (customTurbulenceActive && optimizationCycleGradientClosure != enableTurbulentOptimization)
{
    FatalErrorInFunction
        << "Optimization-Cycle turbulent optimization is guarded by a strict two-key "
        << "runtime gate. Set both "
        << "customTurbulenceProperties.optimizationCycleGradientClosure and "
        << "customTurbulenceProperties.enableTurbulentOptimization to true "
        << "to enable it, or set both to false to keep it disabled."
        << exit(FatalError);
}
if (customTurbulenceActive && turbulenceProModeActive && (optimizationCycleGradientClosure || enableTurbulentOptimization))
{
    FatalErrorInFunction
        << "customTurbulenceProperties.turbulenceProMode=true requires "
        << "optimizationCycleGradientClosure=false and enableTurbulentOptimization=false. "
        << "Disable the strict Optimization-Cycle gate keys when running Pro mode."
        << exit(FatalError);
}
if (customTurbulenceActive && adjointCycleThermalAdjointOnly)
{
    Info<< "Custom kOmegaSST Adjoint-Cycle thermal-adjoint verification mode is "
        << "enabled. The explicit Ta block will be solved with alphaEff, while "
        << "the turbulent momentum adjoint, k-adjoint, omega-adjoint, and "
        << "turbulent sensitivities remain disabled."
        << endl;
}
else if (customTurbulenceActive && adjointCycleFullAdjointVerification)
{
    Info<< "Custom kOmegaSST Adjoint-Cycle full-adjoint verification mode is "
        << "enabled. The explicit Ta and Ua/pa verification blocks will be "
        << "solved, and ka/omegaa proxy fields will be reconstructed from the "
        << "same Stage-B residual blocks. Turbulent sensitivities remain "
        << "disabled."
        << endl;
    Info<< "Adjoint-Cycle verification toggles: reset="
        << (adjointCycleResetAdjointState ? "true" : "false")
        << " rollbackGuards=" << (optAdjointRollbackGuardsEnable ? "true" : "false")
        << " viscosity=" << (adjointCycleUseNuEff ? "nuEff" : "nu")
        << " thermalSource="
        << (adjointCycleIncludeThermalMomentumSource ? "true" : "false")
        << " legacyPressureRelax="
        << (adjointCycleUseLegacyPressureRelax ? "true" : "false")
        << " divergenceGuard="
        << (adjointCycleDivergenceCheck ? "true" : "false")
        << " growthFactor=" << adjointCycleDivergenceGrowthFactor
        << " minPass=" << adjointCycleDivergenceMinCorrectors
        << endl;
}
else if (customTurbulenceActive && optimizationCycleTurbulentOptimizationMode)
{
    Info<< "Custom kOmegaSST Optimization-Cycle turbulent optimization mode is enabled "
        << "through the strict runtime gate. The solver will assemble "
        << "turbulence-aware sensitivities and include explicit residual "
        << "closure derivatives in the optimizer gradients. "
        << "Adjoint reset="
        << (optimizationCycleResetAdjointState ? "true" : "false") << ". "
        << "Adjoint correctors per loop: flow="
        << optimizationCycleFlowAdjointCorrectors
        << ", heat=" << optimizationCycleHeatAdjointCorrectors << ". "
        << "Adjoint magnitude guard: |Uadj|max <= "
        << optimizationCycleAdjointMaxMagnitude << ". "
        << "Closure adjoints="
        << (optimizationCycleTransportedClosureAdjoints ? "transported-PDE" : "proxy-reconstruction")
        << " (invTau=" << optimizationCycleClosureAdjointInvTau
        << ", crossCoupling=" << optimizationCycleClosureCrossCoupling
        << ", formulation=" << optimizationCycleClosureAdjointFormulation
        << ", dualBC=" << (optimizationCycleClosureAdjointDualConsistentBC ? "true" : "false")
        << ", clipWeighting=false"
        << "). "
        << "Porous-aware transport="
        << (porousAwareTurbulenceTransport ? "true" : "false")
        << " (prodExp=" << porousTurbulenceProductionExponent
        << ", prodFloor=" << porousTurbulenceProductionFloor
        << ", sinkExp=" << porousTurbulenceSuppressionExponent
        << ", sinkS=" << porousTurbulenceSuppressionCoeff
        << ", sinkA=" << porousTurbulenceAlphaSinkCoeff << "). "
        << "Bounded closure state="
        << (optimizationCycleClosureConsistencyUseBoundedState ? "true" : "false")
        << " logState=" << (optimizationCycleUseLogClosureState ? "true" : "false")
        << ". "
        << "No-improve rollback guard="
        << (optimizationCycleNoImproveCheck ? "true" : "false")
        << " globalRollbackGuards=" << (optAdjointRollbackGuardsEnable ? "true" : "false")
        << " (minPass=" << optimizationCycleNoImproveMinCorrectors
        << ", patience=" << optimizationCycleNoImprovePatience << "). "
        << "Thermal adjoint mode=" << optThermalAdjointMode
        << " (ta=explicit Ta, tb=legacy turbulent Tb)."
        << endl;
}
else if (customTurbulenceActive && turbulenceProModeActive)
{
    Info<< "Custom kOmegaSST turbulence Pro mode is enabled. "
        << "Primal turbulence fields (k/omega) are solved to resolve U/p/T, "
        << "and closure adjoints (ka/omegaa) are solved/reconstructed for the "
        << "adjoint state, while optimizer updates use laminar-style "
        << "cost/sensitivity/filter/update logic without turbulent sensitivity terms."
        << endl;
    Info<< "Pro mode closure adjoints="
        << (optimizationCycleTransportedClosureAdjoints ? "transported-PDE" : "proxy-reconstruction")
        << " (invTau=" << optimizationCycleClosureAdjointInvTau
        << ", crossCoupling=" << optimizationCycleClosureCrossCoupling
        << ", formulation=" << optimizationCycleClosureAdjointFormulation
        << ", dualBC=" << (optimizationCycleClosureAdjointDualConsistentBC ? "true" : "false")
        << "). "
        << "Pro closure-adjoint solves="
        << (proModeClosureAdjointSolvesActive ? "enabled" : "disabled")
        << ". "
        << "Adjoint correctors per loop: flow=40, heat=40"
        << " (laminar-aligned Pro mode)"
        << "."
        << endl;
}
else if (customTurbulenceActive)
{
    FatalErrorInFunction
        << "Custom kOmegaSST is currently implemented only through the Stage-B "
        << "primal and staged verification paths, plus the strict Optimization-Cycle "
        << "mode and turbulence Pro mode. Select exactly one of adjointCycleThermalAdjointOnly, "
        << "adjointCycleFullAdjointVerification, or set "
        << "both optimizationCycleGradientClosure and enableTurbulentOptimization to true, "
        << "or set turbulenceProMode=true (with both gate keys false)."
        << exit(FatalError);
}
if (customTurbulenceActive)
{
    Info<< "Turbulence transport diffusion scales: primal(k/omega)="
        << turbulencePrimalTransportDiffusionScale
        << ", closureAdjoint(ka/omegaa)="
        << turbulenceAdjointTransportDiffusionScale
        << endl;
}

if
(
    inletTurbulenceMode != "intensityLengthScale"
 && inletTurbulenceMode != "directFields"
)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.inletTurbulenceMode must be "
        << "intensityLengthScale or directFields. Found "
        << inletTurbulenceMode
        << exit(FatalError);
}

wallDist wallDistance(mesh);
const volScalarField& yWall = wallDistance.y();

Info<< "Reading field k\n" << endl;
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField ka
(
    IOobject
    (
        "ka",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroKa", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Reading field omega\n" << endl;
volScalarField omega
(
    IOobject
    (
        "omega",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField omegaa
(
    IOobject
    (
        "omegaa",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroOmegaa", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField kClosureState
(
    IOobject
    (
        "kClosureState",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    k
);
volScalarField omegaClosureState
(
    IOobject
    (
        "omegaClosureState",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    omega
);
volScalarField logKClosureState
(
    IOobject
    (
        "logKClosureState",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroLogKClosureState", dimless, 0.0)
);
volScalarField logOmegaClosureState
(
    IOobject
    (
        "logOmegaClosureState",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroLogOmegaClosureState", dimless, 0.0)
);

const label inletPatchId = mesh.boundaryMesh().findPatchID("inlet");
scalar turbulenceInletKTarget(0.0);
scalar turbulenceInletOmegaTarget(0.0);
if (customTurbulenceActive)
{
    if (inletTurbulenceMode == "intensityLengthScale")
    {
        if (inletTurbulenceIntensity <= 0.0 || inletTurbulenceIntensity > 1.0)
        {
            FatalErrorInFunction
                << "customTurbulenceProperties.turbulenceIntensity must satisfy "
                << "0 < I <= 1. Found I=" << inletTurbulenceIntensity
                << exit(FatalError);
        }
        if (inletTurbulenceLengthScale.value() <= SMALL)
        {
            FatalErrorInFunction
                << "customTurbulenceProperties.turbulenceLengthScale must be "
                << "strictly positive. Found Lt="
                << inletTurbulenceLengthScale.value()
                << exit(FatalError);
        }
        const scalar characteristicLength = inletTurbulenceLengthScale.value();
        const scalar cmuQuarter = Foam::pow(betaStarCoeff.value(), 0.25);

        turbulenceInletKTarget =
            1.5*sqr(inletVelocityMagnitude*inletTurbulenceIntensity);
        turbulenceInletOmegaTarget =
            Foam::sqrt(Foam::max(turbulenceInletKTarget, SMALL))
            /(Foam::max(cmuQuarter*characteristicLength, SMALL));
    }
    else if (inletPatchId >= 0)
    {
        const fvPatchScalarField& inletKPatch = k.boundaryField()[inletPatchId];
        const fvPatchScalarField& inletOmegaPatch = omega.boundaryField()[inletPatchId];
        scalar inletKSum = 0.0;
        scalar inletOmegaSum = 0.0;
        scalar inletFaceCount = 0.0;

        forAll(inletKPatch, facei)
        {
            inletKSum += inletKPatch[facei];
            inletOmegaSum += inletOmegaPatch[facei];
            inletFaceCount += 1.0;
        }

        reduce(inletKSum, sumOp<scalar>());
        reduce(inletOmegaSum, sumOp<scalar>());
        reduce(inletFaceCount, sumOp<scalar>());

        turbulenceInletKTarget = inletKSum/Foam::max(inletFaceCount, SMALL);
        turbulenceInletOmegaTarget =
            inletOmegaSum/Foam::max(inletFaceCount, SMALL);
    }
}

volScalarField nut
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroNut", dimensionSet(0,2,-1,0,0,0,0), 0.0)
);

volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroAlphat", DT.dimensions(), 0.0)
);

volScalarField nuEff
(
    IOobject
    (
        "nuEff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("baseNuEff", nu.dimensions(), nu.value())
);

volScalarField alphaEff
(
    IOobject
    (
        "alphaEff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    DT
);

volScalarField F1SST
(
    IOobject
    (
        "F1SST",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("F1Init", dimless, 1.0)
);

volScalarField F2SST
(
    IOobject
    (
        "F2SST",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("F2Init", dimless, 1.0)
);

volScalarField alphaKField
(
    IOobject
    (
        "alphaKField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    alphaK2Coeff
);

volScalarField alphaOmegaField
(
    IOobject
    (
        "alphaOmegaField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    alphaOmega2Coeff
);

volScalarField betaField
(
    IOobject
    (
        "betaField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    beta2Coeff
);

volScalarField gammaField
(
    IOobject
    (
        "gammaField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    gamma2Coeff
);

volScalarField PkField
(
    IOobject
    (
        "PkField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroPk", dimensionSet(0,2,-3,0,0,0,0), 0.0)
);

volScalarField omegaCrossDiffField
(
    IOobject
    (
        "omegaCrossDiffField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroOmegaCrossDiff", dimensionSet(0,0,-2,0,0,0,0), 0.0)
);

volScalarField momentumResistanceBlock
(
    IOobject
    (
        "momentumResistanceBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha
);

volScalarField momentumDiffusivityBlock
(
    IOobject
    (
        "momentumDiffusivityBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    nuEff
);

volScalarField thermalDiffusivityBlock
(
    IOobject
    (
        "thermalDiffusivityBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaEff
);

volScalarField kDiffusionBlock
(
    IOobject
    (
        "kDiffusionBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    nu + alphaKField*nut
);

volScalarField omegaDiffusionBlock
(
    IOobject
    (
        "omegaDiffusionBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    nu + alphaOmegaField*nut
);

volScalarField kDestructionBlock
(
    IOobject
    (
        "kDestructionBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    betaStarCoeff*omega
);

volScalarField omegaProductionBlock
(
    IOobject
    (
        "omegaProductionBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    gammaField*PkField/Foam::max(nut, nutFloor)
);

volScalarField omegaDestructionBlock
(
    IOobject
    (
        "omegaDestructionBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    betaField*omega
);

volScalarField turbulenceDampingBlock
(
    IOobject
    (
        "turbulenceDampingBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    turbulenceDampingCoeff*alpha
);

volScalarField dNutDkBlock
(
    IOobject
    (
        "dNutDkBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroDNutDk", dimensionSet(0,0,1,0,0,0,0), 0.0)
);

volScalarField dNutDOmegaBlock
(
    IOobject
    (
        "dNutDOmegaBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroDNutDOmega", dimensionSet(0,2,0,0,0,0,0), 0.0)
);

volScalarField dAlphaEffDkBlock
(
    IOobject
    (
        "dAlphaEffDkBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroDAlphaEffDk", dimensionSet(0,0,1,0,0,0,0), 0.0)
);

volScalarField dAlphaEffDOmegaBlock
(
    IOobject
    (
        "dAlphaEffDOmegaBlock",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroDAlphaEffDOmega", dimensionSet(0,2,0,0,0,0,0), 0.0)
);
