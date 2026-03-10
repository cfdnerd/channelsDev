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
Switch proModeAdjointStabilityPack
(
    customTurbulenceProperties.lookupOrDefault<Switch>("proModeAdjointStabilityPack", true)
);
word proModeAdjointClosureModel
(
    customTurbulenceProperties.lookupOrDefault<word>("proModeAdjointClosureModel", "sst")
);
word wallTreatmentMode
(
    customTurbulenceProperties.lookupOrDefault<word>("wallTreatmentMode", "legacyFixedKOmega")
);
const bool wallTreatmentModeLegacyFixedKOmega =
    wallTreatmentMode == word("legacyFixedKOmega");
const bool wallTreatmentModeWallFunctionCompatible =
    wallTreatmentMode == word("wallFunctionCompatible");
const bool turbulenceProModeActive =
    customTurbulenceActive
 && static_cast<bool>(turbulenceProMode);
const bool proModeClosureAdjointSolvesActive =
    turbulenceProModeActive
 && static_cast<bool>(proModeClosureAdjointSolves);
const bool proModeAdjointStabilityPackActive =
    turbulenceProModeActive
 && static_cast<bool>(proModeAdjointStabilityPack);
const bool adjointRollbackGuardsEffective =
    static_cast<bool>(optAdjointRollbackGuardsEnable)
 || proModeAdjointStabilityPackActive;
const bool proModeSAAdjointClosureSelected =
    turbulenceProModeActive
 && proModeAdjointClosureModel == word("sa");
const bool proModeSAAdjointClosureActive =
    proModeClosureAdjointSolvesActive
 && proModeSAAdjointClosureSelected;
// Pro-mode stability pack (hardcoded): activated by one dict toggle.
const label proModeFlowAdjointCorrectors =
    proModeAdjointStabilityPackActive ? label(24) : label(40);
const label proModeHeatAdjointCorrectors =
    proModeAdjointStabilityPackActive ? label(24) : label(40);
const bool proModeResetAdjointState = false;
const bool proModeAdjointNoImproveCheck = proModeAdjointStabilityPackActive;
const label proModeAdjointNoImproveMinCorrectors(10);
const label proModeAdjointNoImprovePatience(8);
const scalar proModeAdjointMaxMagnitude(1e8);
const bool proModeAdjointUseNuEff = proModeAdjointStabilityPackActive;
const bool proModeAdjointUseLegacyPressureRelax = proModeAdjointStabilityPackActive;
const label proModeAdjointWarmupIterations(120);
const label proModeAdjointRampRollbackFreeWindow(8);
const scalar proModePseudoInvTauFloorRatio(0.20);
const scalar proModeUaPseudoInvTauBase
(
    proModeSAAdjointClosureSelected ? 16.0 : 12.0
);
const scalar proModeUbPseudoInvTauBase
(
    proModeSAAdjointClosureSelected ? 30.0 : 20.0
);
const scalar proModeClosurePseudoInvTauScale(0.75);
const scalar proModeUaEqnRelaxMin(0.08);
const scalar proModeUaEqnRelaxMax(0.20);
const scalar proModeUbEqnRelaxMin(0.05);
const scalar proModeUbEqnRelaxMax(0.15);
const scalar proModeClosureEqnRelaxMin(0.10);
const scalar proModeClosureEqnRelaxMax(0.20);
const scalar proModePrimalTransportDiffusionMinScale(2.0);
const scalar proModeAdjointTransportDiffusionMinScale(3.0);
const bool proModeTransportedClosureAdjoints = true;
const word proModeClosureAdjointFormulation("transportedResidualDraft");
const bool proModeClosureAdjointDualConsistentBC(true);
const scalar proModeClosureAdjointInvTau(3.0);
const scalar proModeClosureCrossCoupling(0.15);
const bool proModeRobustLinearSolverActive = proModeAdjointStabilityPackActive;
const bool proModeRobustScalarSolverActive = proModeRobustLinearSolverActive;
dictionary proModeRobustScalarSolverControls;
dictionary proModeRobustNonsymmetricSolverControls;
dictionary proModeRobustSymmetricSolverControls;
if (proModeRobustLinearSolverActive)
{
    // Generic scalar fallback robust for both symmetric/non-symmetric scalar systems.
    proModeRobustScalarSolverControls.add("solver", word("smoothSolver"));
    proModeRobustScalarSolverControls.add("smoother", word("symGaussSeidel"));
    proModeRobustScalarSolverControls.add("tolerance", 1e-10);
    proModeRobustScalarSolverControls.add("relTol", 0.0);
    proModeRobustScalarSolverControls.add("maxIter", 20000);

    // Non-symmetric systems: transport-dominated adjoint equations (Ua/Ub/ka/omegaa).
    proModeRobustNonsymmetricSolverControls.add("solver", word("PBiCGStab"));
    proModeRobustNonsymmetricSolverControls.add("preconditioner", word("DILU"));
    proModeRobustNonsymmetricSolverControls.add("tolerance", 1e-10);
    proModeRobustNonsymmetricSolverControls.add("relTol", 0.0);
    proModeRobustNonsymmetricSolverControls.add("maxIter", 20000);

    // Symmetric systems: pressure-correction adjoint Poisson solves (pa/pb).
    // Use multigrid here to reduce high-iteration/stalled adjoint pressure passes
    // observed with Krylov+DIC on ill-conditioned late-stage designs.
    proModeRobustSymmetricSolverControls.add("solver", word("GAMG"));
    proModeRobustSymmetricSolverControls.add("tolerance", 1e-10);
    proModeRobustSymmetricSolverControls.add("relTol", 0.0);
    proModeRobustSymmetricSolverControls.add("smoother", word("GaussSeidel"));
    proModeRobustSymmetricSolverControls.add("cacheAgglomeration", true);
    proModeRobustSymmetricSolverControls.add("agglomerator", word("faceAreaPair"));
    proModeRobustSymmetricSolverControls.add("mergeLevels", 1);
    proModeRobustSymmetricSolverControls.add("nPreSweeps", 0);
    proModeRobustSymmetricSolverControls.add("nPostSweeps", 2);
    proModeRobustSymmetricSolverControls.add("nFinestSweeps", 2);
    proModeRobustSymmetricSolverControls.add("nCellsInCoarsestLevel", 20);
    proModeRobustSymmetricSolverControls.add("maxIter", 20000);
}
Switch adjointCycleThermalAdjointOnly
(
    customTurbulenceProperties.lookupOrDefault<Switch>("adjointCycleThermalAdjointOnly", false)
);
Switch adjointCycleFullAdjointVerification
(
    customTurbulenceProperties.lookupOrDefault<Switch>("adjointCycleFullAdjointVerification", false)
);
const bool optimizationCycleGradientClosure =
    customTurbulenceActive && !turbulenceProModeActive;
const bool enableTurbulentOptimization =
    customTurbulenceActive && !turbulenceProModeActive;
if
(
    customTurbulenceActive
 && (
        customTurbulenceProperties.found("optimizationCycleGradientClosure")
     || customTurbulenceProperties.found("enableTurbulentOptimization")
    )
)
{
    Info<< "Ignoring customTurbulenceProperties.optimizationCycleGradientClosure and "
        << "customTurbulenceProperties.enableTurbulentOptimization. "
        << "Both are now derived from turbulenceProMode (false->true,true | true->false,false)."
        << endl;
}
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
const scalar turbulencePrimalTransportDiffusionScaleActive =
    (turbulenceProModeActive && proModeAdjointStabilityPackActive)
  ? Foam::max(turbulencePrimalTransportDiffusionScale, proModePrimalTransportDiffusionMinScale)
  : turbulencePrimalTransportDiffusionScale;
const scalar turbulenceAdjointTransportDiffusionScaleActive =
    (turbulenceProModeActive && proModeAdjointStabilityPackActive)
  ? Foam::max(turbulenceAdjointTransportDiffusionScale, proModeAdjointTransportDiffusionMinScale)
  : turbulenceAdjointTransportDiffusionScale;
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
    customTurbulenceProperties.lookupOrDefault<Switch>("optimizationCycleTransportedClosureAdjoints", true)
);
word optimizationCycleClosureAdjointFormulation
(
    customTurbulenceProperties.lookupOrDefault<word>("optimizationCycleClosureAdjointFormulation", "transportedResidualDraft")
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
const bool closureAdjointTransportedActive =
    turbulenceProModeActive
 ? (
      proModeAdjointStabilityPackActive
      ? static_cast<bool>(proModeTransportedClosureAdjoints)
      : static_cast<bool>(optimizationCycleTransportedClosureAdjoints)
   )
 : static_cast<bool>(optimizationCycleTransportedClosureAdjoints);
const word closureAdjointFormulationActive =
    turbulenceProModeActive
 ? (
      proModeAdjointStabilityPackActive
      ? proModeClosureAdjointFormulation
      : optimizationCycleClosureAdjointFormulation
   )
 : optimizationCycleClosureAdjointFormulation;
const bool closureAdjointDualConsistentBCActive =
    turbulenceProModeActive
 ? (
      proModeAdjointStabilityPackActive
      ? static_cast<bool>(proModeClosureAdjointDualConsistentBC)
      : static_cast<bool>(optimizationCycleClosureAdjointDualConsistentBC)
   )
 : static_cast<bool>(optimizationCycleClosureAdjointDualConsistentBC);
const scalar closureAdjointInvTauActive =
    turbulenceProModeActive
 ? (
      proModeAdjointStabilityPackActive
      ? proModeClosureAdjointInvTau
      : optimizationCycleClosureAdjointInvTau
   )
 : optimizationCycleClosureAdjointInvTau;
const scalar closureAdjointCrossCouplingActive =
    turbulenceProModeActive
 ? (
      proModeAdjointStabilityPackActive
      ? proModeClosureCrossCoupling
      : optimizationCycleClosureCrossCoupling
   )
 : optimizationCycleClosureCrossCoupling;
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
    customTurbulenceActive
 && wallTreatmentMode != "legacyFixedKOmega"
 && wallTreatmentMode != "wallFunctionCompatible"
)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.wallTreatmentMode must be "
        << "legacyFixedKOmega or wallFunctionCompatible. Found "
        << wallTreatmentMode
        << exit(FatalError);
}
if
(
    customTurbulenceActive
 && optimizationCycleClosureAdjointFormulation != "proxyNormalized"
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
if
(
    customTurbulenceActive
 && proModeClosureAdjointFormulation != "proxyNormalized"
 && proModeClosureAdjointFormulation != "transportedResidualDraft"
 && proModeClosureAdjointFormulation != "transportedCrossCoupledV3"
)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.proModeClosureAdjointFormulation must "
        << "be proxyNormalized, transportedResidualDraft, or "
        << "transportedCrossCoupledV3. Found "
        << proModeClosureAdjointFormulation
        << exit(FatalError);
}
if
(
    proModeAdjointClosureModel != "sst"
 && proModeAdjointClosureModel != "sa"
)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.proModeAdjointClosureModel must be "
        << "sst or sa. Found "
        << proModeAdjointClosureModel
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
if (proModeFlowAdjointCorrectors < 1 || proModeHeatAdjointCorrectors < 1)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.proModeFlowAdjointCorrectors and "
        << "proModeHeatAdjointCorrectors must both be >= 1. Found flow="
        << proModeFlowAdjointCorrectors
        << ", heat=" << proModeHeatAdjointCorrectors
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
if (proModeClosureAdjointInvTau <= SMALL)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.proModeClosureAdjointInvTau must be > 0. Found "
        << proModeClosureAdjointInvTau
        << exit(FatalError);
}
if (proModeClosureCrossCoupling < 0.0 || proModeClosureCrossCoupling > 1.0)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.proModeClosureCrossCoupling must be in [0, 1]. Found "
        << proModeClosureCrossCoupling
        << exit(FatalError);
}
if (proModeAdjointMaxMagnitude <= SMALL)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.proModeAdjointMaxMagnitude must be > 0. Found "
        << proModeAdjointMaxMagnitude
        << exit(FatalError);
}
if (proModeAdjointNoImproveMinCorrectors < 1)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.proModeAdjointNoImproveMinCorrectors must be >= 1. Found "
        << proModeAdjointNoImproveMinCorrectors
        << exit(FatalError);
}
if
(
    proModeAdjointNoImproveMinCorrectors
 > Foam::max(proModeFlowAdjointCorrectors, proModeHeatAdjointCorrectors)
)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.proModeAdjointNoImproveMinCorrectors must be <= "
        << "max(proModeFlowAdjointCorrectors, proModeHeatAdjointCorrectors). Found minPass="
        << proModeAdjointNoImproveMinCorrectors
        << ", flow=" << proModeFlowAdjointCorrectors
        << ", heat=" << proModeHeatAdjointCorrectors
        << exit(FatalError);
}
if (proModeAdjointNoImprovePatience < 1)
{
    FatalErrorInFunction
        << "customTurbulenceProperties.proModeAdjointNoImprovePatience must be >= 1. Found "
        << proModeAdjointNoImprovePatience
        << exit(FatalError);
}
if
(
    (adjointCycleThermalAdjointOnly && adjointCycleFullAdjointVerification)
 || (adjointCycleThermalAdjointOnly && optimizationCycleTurbulentOptimizationMode)
 || (adjointCycleFullAdjointVerification && optimizationCycleTurbulentOptimizationMode)
 || (adjointCycleThermalAdjointOnly && turbulenceProModeActive)
 || (adjointCycleFullAdjointVerification && turbulenceProModeActive)
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
        << " rollbackGuards=" << (adjointRollbackGuardsEffective ? "true" : "false")
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
        << "via turbulenceProMode=false. The solver will assemble "
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
        << " globalRollbackGuards=" << (adjointRollbackGuardsEffective ? "true" : "false")
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
        << "and closure adjoints ("
        << (proModeSAAdjointClosureSelected ? "ka (SA-family)" : "ka/omegaa")
        << ") are solved/reconstructed for the "
        << "adjoint state, while optimizer updates use laminar-style "
        << "cost/sensitivity/filter/update logic without turbulent sensitivity terms."
        << endl;
    Info<< "Pro mode closure adjoints="
        << (closureAdjointTransportedActive ? "transported-PDE" : "proxy-reconstruction")
        << " (invTau=" << closureAdjointInvTauActive
        << ", crossCoupling=" << closureAdjointCrossCouplingActive
        << ", formulation=" << closureAdjointFormulationActive
        << ", dualBC=" << (closureAdjointDualConsistentBCActive ? "true" : "false")
        << ", model="
        << (proModeSAAdjointClosureSelected ? "sa(ka-only, nuTilda-family)" : "sst(ka+omegaa)")
        << "). "
        << "Pro closure-adjoint solves="
        << (proModeClosureAdjointSolvesActive ? "enabled" : "disabled")
        << ". "
        << "StabilityPack="
        << (proModeAdjointStabilityPackActive ? "on" : "off")
        << ". "
        << "Adjoint correctors per loop: flow=" << proModeFlowAdjointCorrectors
        << ", heat=" << proModeHeatAdjointCorrectors
        << ". Pro adjoint tuning: reset="
        << (proModeResetAdjointState ? "true" : "false")
        << " noImproveGuard=" << (proModeAdjointNoImproveCheck ? "true" : "false")
        << " (minPass=" << proModeAdjointNoImproveMinCorrectors
        << ", patience=" << proModeAdjointNoImprovePatience << ")"
        << " rollbackGuards=" << (adjointRollbackGuardsEffective ? "true" : "false")
        << " |Uadj|max<=" << proModeAdjointMaxMagnitude
        << " viscosity=" << (proModeAdjointUseNuEff ? "nuEff" : "nu")
        << " legacyPressureRelax="
        << (proModeAdjointUseLegacyPressureRelax ? "true" : "false")
        << " warmupIters=" << proModeAdjointWarmupIterations
        << " relax(Ua)=" << proModeUaEqnRelaxMin << "->" << proModeUaEqnRelaxMax
        << " relax(Ub)=" << proModeUbEqnRelaxMin << "->" << proModeUbEqnRelaxMax
        << " relax(closure)=" << proModeClosureEqnRelaxMin << "->" << proModeClosureEqnRelaxMax
        << " pseudoFloor=" << proModePseudoInvTauFloorRatio
        << " robustNS=PBiCGStab/DILU"
        << " robustSym=PCG/DIC"
        << "."
        << endl;
}
else if (customTurbulenceActive)
{
    FatalErrorInFunction
        << "Custom kOmegaSST is currently implemented only through the Stage-B "
        << "primal and staged verification paths, plus the Optimization-Cycle "
        << "mode and turbulence Pro mode. Select exactly one of adjointCycleThermalAdjointOnly, "
        << "adjointCycleFullAdjointVerification, or use turbulenceProMode "
        << "(false=Optimization-Cycle, true=Pro mode)."
        << exit(FatalError);
}
if (customTurbulenceActive)
{
    Info<< "Turbulence transport diffusion scales: primal(k/omega)="
        << turbulencePrimalTransportDiffusionScaleActive
        << ", closureAdjoint(ka/omegaa)="
        << turbulenceAdjointTransportDiffusionScaleActive
        << ", wallTreatmentMode=" << wallTreatmentMode
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

if (customTurbulenceActive)
{
    auto isWallFunctionPatchType =
    [&](const fvPatchScalarField& patchField) -> bool
    {
        const std::string patchType(patchField.type());
        return patchType.find("WallFunction") != std::string::npos;
    };

    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].type() != "wall")
        {
            continue;
        }

        const fvPatchScalarField& kPatch = k.boundaryField()[patchi];
        const fvPatchScalarField& omegaPatch = omega.boundaryField()[patchi];
        const fvPatchScalarField& nutPatch = nut.boundaryField()[patchi];

        const bool kWallFunctionPatch = isWallFunctionPatchType(kPatch);
        const bool omegaWallFunctionPatch = isWallFunctionPatchType(omegaPatch);
        const bool nutWallFunctionPatch = isWallFunctionPatchType(nutPatch);
        const bool kFixedValuePatch = kPatch.type() == word("fixedValue");
        const bool omegaFixedValuePatch = omegaPatch.type() == word("fixedValue");
        const bool nutFixedValuePatch = nutPatch.type() == word("fixedValue");

        if
        (
            kWallFunctionPatch
         || omegaWallFunctionPatch
         || nutWallFunctionPatch
        )
        {
            FatalErrorInFunction
                << "OpenFOAM wall-function BCs are not supported in this solver branch because "
                << "no standard turbulenceModel object is registered in objectRegistry. Patch="
                << mesh.boundary()[patchi].name()
                << " has kType=" << kPatch.type()
                << ", omegaType=" << omegaPatch.type()
                << ", nutType=" << nutPatch.type()
                << ". Use fixedValue BCs and rely on wallTreatmentMode="
                << wallTreatmentMode
                << " internal wall handling."
                << exit(FatalError);
        }

        if
        (
            wallTreatmentModeLegacyFixedKOmega
         && (!kFixedValuePatch || !omegaFixedValuePatch)
        )
        {
            WarningInFunction
                << "wallTreatmentMode=legacyFixedKOmega overrides wall k/omega values every iteration. "
                << "Expected fixedValue BCs on patch="
                << mesh.boundary()[patchi].name()
                << " (kType=" << kPatch.type()
                << ", omegaType=" << omegaPatch.type() << ")."
                << endl;
        }

        if
        (
            wallTreatmentModeWallFunctionCompatible
         && (!kFixedValuePatch || !omegaFixedValuePatch || !nutFixedValuePatch)
        )
        {
            FatalErrorInFunction
                << "wallTreatmentMode=wallFunctionCompatible requires fixedValue BCs for "
                << "k, omega, and nut (solver-internal wall-function-compatible treatment). Patch="
                << mesh.boundary()[patchi].name()
                << " has kType=" << kPatch.type()
                << ", omegaType=" << omegaPatch.type()
                << ", nutType=" << nutPatch.type()
                << exit(FatalError);
        }
    }
}

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
    gammaField*PkField
    /
    (
      static_cast<bool>(optClosureSmoothDenominatorEnable)
      ? Foam::sqrt
        (
          sqr(nut)
          + sqr
            (
              dimensionedScalar
              (
                "nutRegularizationFloor",
                nut.dimensions(),
                Foam::max
                (
                  nutFloor.value(),
                  optClosureNutRegularizationScale*Foam::max(nu.value(), SMALL)
                )
              )
            )
        )
      : Foam::max
        (
          nut,
          dimensionedScalar
          (
            "nutRegularizationFloor",
            nut.dimensions(),
            Foam::max
            (
              nutFloor.value(),
              optClosureNutRegularizationScale*Foam::max(nu.value(), SMALL)
            )
          )
        )
    )
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
