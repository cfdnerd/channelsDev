# Optimizer Patch Backlog: Continuation and Model-Agreement

## Goal
Improve optimizer progress and gradient quality for:
- Item 5: continuation-induced non-smoothness
- Item 3: predicted-vs-actual objective mismatch

## Patch 1: Continuation Rate Control (Configurable Schedules)
- Problem: continuation ramps are currently hard-coded and can sharpen the landscape too fast.
- Change:
- Make `alphaMax`, `qu`, and `del` ramps fully configurable from `tuneOptParameters`.
- Replace fixed constants in `update.H` and `filter_x.H` with dict-driven rates/caps.
- Target files:
- `turbulenceOptimizer/src/update.H`
- `turbulenceOptimizer/src/filter_x.H`
- `turbulenceOptimizer/src/createFields.H`
- `turbulenceOptimizer/app/constant/tuneOptParameters`
- Expected benefit:
- Smoother transition from relaxed to sharp design mapping.
- Fewer abrupt drops in trust ratio and fewer consistency spikes after continuation jumps.

## Patch 2: Agreement-Aware Continuation Governor
- Problem: schedule currently advances by iteration count, even when local linear model agreement is poor.
- Change:
- Use previous-iteration `trustRatio` and `gradientConsistencyRatio` to modulate continuation schedule speed.
- Rules should affect schedule progression only (advance/hold/slow), not step acceptance.
- Keep this separate from gate/accept logic.
- Target files:
- `turbulenceOptimizer/src/update.H`
- `turbulenceOptimizer/src/debugOptimizer.H` (read existing metrics/state only)
- `turbulenceOptimizer/src/opt_initialization.H` (state variables if needed)
- `turbulenceOptimizer/src/createFields.H`
- `turbulenceOptimizer/app/constant/tuneOptParameters`
- Expected benefit:
- Prevents sharpening the model while gradients are directionally inconsistent.
- Better alignment between predicted and realized objective changes.

## Patch 3: Objective-Gradient Temporal Filtering (EMA)
- Problem: iteration-to-iteration objective sensitivity can flip direction under noisy/stiff regimes.
- Change:
- Add exponential moving average filtering for objective gradient `dfdx`.
- Blend current `dfdx` with previous filtered `dfdx` before MMA.
- Optionally make smoothing strength adaptive to `gradientConsistencyRatio`.
- Target files:
- `turbulenceOptimizer/src/sensitivity.H`
- `turbulenceOptimizer/src/opt_initialization.H` (persistent filtered gradient storage)
- `turbulenceOptimizer/src/createFields.H`
- `turbulenceOptimizer/app/constant/tuneOptParameters`
- Expected benefit:
- Reduces directional jitter in `dfdx`.
- Improves trust-ratio stability and lowers model-agreement oscillations.

## Patch 4: Secant-Based Objective Gradient Magnitude Correction
- Problem: objective model often has wrong scale even when direction is partly correct.
- Change:
- Compute bounded scalar correction from previous predicted-vs-actual objective reduction.
- Apply correction to objective gradient magnitude for the next iteration.
- Keep correction bounded to avoid instability.
- Target files:
- `turbulenceOptimizer/src/debugOptimizer.H` (compute correction factor)
- `turbulenceOptimizer/src/sensitivity.H` (apply factor to objective gradient)
- `turbulenceOptimizer/src/opt_initialization.H` (store previous correction state)
- `turbulenceOptimizer/src/createFields.H`
- `turbulenceOptimizer/app/constant/tuneOptParameters`
- Expected benefit:
- Better first-order model scaling.
- Lower persistent mismatch in predicted vs realized objective reduction.

## Patch 5: Activate Adaptive MMA Move Limit
- Problem: fixed move limit is less robust when local model quality varies strongly over continuation.
- Change:
- Enable adaptive move limits with conservative defaults.
- Keep move-limit adaptation active as a damping mechanism only.
- Target files:
- `turbulenceOptimizer/app/constant/tuneOptParameters`
- Existing code path already present in `turbulenceOptimizer/src/update.H`.
- Expected benefit:
- Limits damage from temporarily poor local models.
- Stabilizes optimization trajectory while other formulation patches improve fidelity.

## Suggested Rollout Order
- 1) Patch 1
- 2) Patch 5
- 3) Patch 2
- 4) Patch 3
- 5) Patch 4

## Primary Validation Metrics (from `debugOptimizer.H` logs)
- `trustRatio`
- `objectiveDelta`
- `predictedObjectiveDelta`
- `gradientConsistencyRatio`
- `objectiveConsistencyRatio`
- convergence trend of objective and active constraints

