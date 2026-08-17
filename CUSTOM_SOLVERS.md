# Custom Solvers — Bed Load & Bridge Scour

This fork adds custom OpenFOAM/CFDEM solvers on top of upstream CFDEMcoupling-PUBLIC for coupled CFD-DEM simulation of sediment transport: bed load transport and local scour around bridge piers. See `applications/solvers/` for the full CFDEMcoupling-PUBLIC solver set; the ones below are the additions.

These solvers extend `cfdemSolverPiso` with an additional phase-turbulence closure (`DPMTurbulenceModels`, providing a custom `kEpsilon_TD` RAS model) and case-specific bed-elevation / suspension-tracking logic (`deltaY*.H`), used to study particle entrainment and suspension in turbulent open-channel flow.

## Solvers

| Solver | Notes |
|---|---|
| `PhasecfdemSolverPiso` | Baseline custom solver: `cfdemSolverPiso` + `DPMIncompressibleTurbulenceModel` phase-turbulence closure. |
| `PhasecfdemSolverPiso_turb` | Turbulence-closure development iteration. |
| `PhasecfdemSolverPiso_NUMAP` | Variant used for the NUMAP campaign. |
| `PhasecfdemSolverPiso_20260108` | Dated development snapshot. |
| `PhasecfdemSolverPiso_20260502` | Dated development snapshot. |
| `PhasecfdemSolverPiso_20260510` | Dated development snapshot; adds `createFields_UP.H`. |
| `PhasecfdemSolverPiso_20260605` | Dated development snapshot. |
| `PhasecfdemSolverPiso_20260731_BridgeScour` | Solver used for the bridge-pier local-scour case studies. |

Earlier, superseded iterations (`PhasecfdemSolverPiso_old`, `_20260108_backup`, `_20260108 - Copy`, `_turb_trial`) are kept for reference under `applications/solvers/_archive/`.

## Building

With the CFDEM environment sourced (this repo's own `README.md` covers the full CFDEMcoupling-PUBLIC install):

```bash
cd $CFDEM_SOLVER_DIR/applications/solvers/<solverName>
wmake ./DPMTurbulenceModels   # build the custom turbulence library first
wmake .                        # build the solver
```

Developed as part of ongoing CFD-DEM research into bed load transport and bridge-pier local scour, built on OpenFOAM-2506 + CFDEMcoupling-PUBLIC.
