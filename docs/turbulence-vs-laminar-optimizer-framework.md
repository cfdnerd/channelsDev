# Turbulence vs Laminar Optimizer Framework: 8 Major Differences

## Scope
This comparison covers optimizer-side framework behavior only (adjoints, sensitivities, filtering, constraints, continuation, guards/repair, diagnostics). Primal physics-solver setup is intentionally excluded.

## 1) Optimization Pipeline Structure
Laminar uses a single fixed optimization path each outer iteration:
`Tb -> Ub -> Ua -> costfunction -> sensitivity`.

Turbulence uses mode-gated execution paths, including:
- thermal-adjoint-only verification,
- full-adjoint verification,
- turbulent optimization mode,
- baseline fallback.

Key files:
- `laminarOptimizer/src/MTO_ThermalFluid.C`
- `turbulenceOptimizer/src/MTO_ThermalFluid.C`

## 2) Thermal Adjoint Formulation Options
Laminar is effectively fixed to `Tb`.

Turbulence supports thermal-adjoint mode selection:
- `ta`
- `tbOriginal`
- `tbLaminar`

This changes which thermal adjoint drives momentum adjoints and sensitivity paths.

Key files:
- `turbulenceOptimizer/src/createFields.H`
- `turbulenceOptimizer/src/Adjoint_Ta.H`
- `turbulenceOptimizer/src/AdjointHeat_Tb.H`
- `turbulenceOptimizer/src/AdjointHeat_TbLaminarIdentical.H`

## 3) Closure Adjoint Extension (k/omega)
Laminar has no turbulence-closure adjoint equations.

Turbulence adds explicit closure adjoints (`ka`, `omegaa`) and closure formulation controls (proxy vs transported variants), with dual-consistency options.

Key files:
- `turbulenceOptimizer/src/Adjoint_k.H`
- `turbulenceOptimizer/src/Adjoint_omega.H`
- `turbulenceOptimizer/src/readTurbulenceProperties.H`

## 4) Sensitivity Framework Depth
Laminar uses direct objective/power/volume sensitivity assembly.

Turbulence adds:
- turbulent closure contributions,
- gradient L2 normalization,
- sensitivity quality scoring/gating,
- optional reuse of prior gradient snapshots,
- adaptive MMA move-limit handling from sensitivity health.

Key files:
- `laminarOptimizer/src/sensitivity.H`
- `turbulenceOptimizer/src/sensitivity.H`
- `turbulenceOptimizer/src/TurbulentSensitivity.H`

## 5) Filtering and Projection Robustness
Both use PDE filtering + Heaviside projection.

Turbulence adds stronger safeguards:
- finite/NaN-safe normalization,
- denominator protection in projection derivatives,
- controlled ramp of projection sharpness tied to continuation state.

Key files:
- `laminarOptimizer/src/filter_chainrule.H`
- `turbulenceOptimizer/src/filter_chainrule.H`
- `laminarOptimizer/src/filter_x.H`
- `turbulenceOptimizer/src/filter_x.H`

## 6) Continuation, Move-Limit, and Recovery Control
Laminar continuation/update logic is comparatively simple (scheduled `alphaMax`, `qu`, and interpolation updates).

Turbulence introduces adaptive and guard-aware control:
- continuation freeze/advance logic,
- recovery mode,
- adaptive MMA asymptotes and move limits,
- multiphysics-health-coupled ramp controls.

Key files:
- `laminarOptimizer/src/update.H`
- `turbulenceOptimizer/src/update.H`

## 7) Objective/Constraint Framework Breadth
Laminar MMA setup is two constraints (power, volume).

Turbulence can extend to thermal envelope constraints and includes bounded thermal objective machinery (KS-style envelope quantities), in addition to power and volume.

Key files:
- `laminarOptimizer/src/opt_initialization.H`
- `turbulenceOptimizer/src/opt_initialization.H`
- `turbulenceOptimizer/src/costfunction.H`

## 8) Guard, Rollback, Repair, and Diagnostics Infrastructure
Laminar has lighter optimizer diagnostics and stop logic.

Turbulence has broader protective infrastructure:
- adjoint rollback guards (magnitude/no-improve/growth checks),
- optional first-failure snapshot persistence,
- optional bounded omega limiter in closure model,
- richer debug logs and health metrics.

Key files:
- `turbulenceOptimizer/src/AdjointFlow_Ua.H`
- `turbulenceOptimizer/src/AdjointHeat_Ub.H`
- `turbulenceOptimizer/src/Primal_k.H`
- `turbulenceOptimizer/src/Primal_omega.H`
- `turbulenceOptimizer/src/CustomKOmegaSSTModel.H`
- `turbulenceOptimizer/src/debugOptimizer.H`
- `turbulenceOptimizer/src/opt_initialization.H`

## Notes on Current Runtime Toggles (Turbulence)
Recent turbulence `tuneOptParameters` controls include:
- `optAdjointRollbackGuardsEnable`
- `optFirstFailureSnapshotEnable`
- `optClosureBoundedOmegaLimiterEnable`

Key file:
- `turbulenceOptimizer/app/constant/tuneOptParameters`
