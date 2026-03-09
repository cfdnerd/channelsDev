# High-Re Wall-Function Fixes (1-4) and Literature Basis

This note documents the implemented fixes for high-Re turbulent topology optimization and the references used to justify the formulation choices.

## Fix 1: Runtime wall-treatment switch

- Implemented `customTurbulenceProperties.wallTreatmentMode`:
  - `legacyFixedKOmega`
  - `wallFunctionCompatible`
- Added strict runtime validation and mode reporting.

Rationale:
- High-Re TO with Darcy/porous penalization requires explicit treatment of near-wall modeling and wall-distance effects in turbulence closure.
- Reference: Wu & Zhang (2024), *Aerospace* 11(7):525, doi:10.3390/aerospace11070525.

## Fix 2: Do not override wall k/omega in wall-function mode

- Legacy mode keeps explicit wall overwrite (`k=kMin`, `omega~nu/y^2`).
- Wall-function-compatible mode now preserves BC-owned wall values (`correctBoundaryConditions()` only).
- Added `nut` and `alphat` boundary correction in wall-function mode.

Rationale:
- Frozen/shortcut closure treatments can produce inconsistent gradients; consistent turbulence treatment is required for stable adjoint-driven TO.
- Reference: Dilgen et al. (2018), *Computer Methods in Applied Mechanics and Engineering* 347, doi:10.1016/j.cma.2017.11.029.

## Fix 3: Case-level wall-function BC workflow

- Updated templates:
  - `0/k`: `kqRWallFunction`
  - `0/omega`: `omegaWallFunction`
  - new `0/nut`: `nutkWallFunction`
- Added startup BC consistency checks per wall patch in wall-function mode.

Rationale:
- High-Re wall-function workflows must keep turbulence-wall BCs consistent across `k`, `omega`, and `nut` to avoid nonphysical wall closure behavior.
- References:
  - OpenFOAM wall-function BC docs (`kqRWallFunction`, `nutkWallFunction`).
  - Yoon (2016), *Computer Methods in Applied Mechanics and Engineering* 303, doi:10.1016/j.cma.2016.01.014.

## Fix 4: y+ diagnostics and warnings

- Added tuning controls:
  - `optWallFunctionYPlusDiagnosticsEnable`
  - `optWallFunctionYPlusMin`
  - `optWallFunctionYPlusMax`
  - `optWallFunctionYPlusInBandWarnFraction`
- Added global/patch y+ stats and warning triggers to table + JSON diagnostics.

Rationale:
- Coarse-mesh high-Re TO requires active checks that wall-function assumptions are being met during optimization.
- References:
  - Wu & Zhang (2024), *Aerospace* 11(7):525, doi:10.3390/aerospace11070525.
  - Bayat et al. (2026, preprint), arXiv:2601.02202, doi:10.48550/arXiv.2601.02202.

## Frozen-Turbulence Context

- Frozen turbulence paths were removed from the turbulent TO closure-adjoint flow and replaced by transported closure-adjoint requirements.
- This follows published evidence that frozen turbulence can yield inconsistent and even incorrect-sign sensitivities.
- Reference: Alexandersen et al. (2020), *Fluids* 5(1):29, doi:10.3390/fluids5010029.
