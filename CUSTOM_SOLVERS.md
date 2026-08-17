# Custom Solver — Bed Load & Bridge Scour

This fork adds a custom OpenFOAM/CFDEM solver, `PhasecfdemSolverPiso`, on top of upstream CFDEMcoupling-PUBLIC for coupled CFD-DEM simulation of sediment transport: bed load transport and local scour around bridge piers.

`PhasecfdemSolverPiso` extends `cfdemSolverPiso` with an additional phase-turbulence closure (`DPMTurbulenceModels`, providing a custom `kEpsilon_TD` RAS model) and bed-interface tracking logic (`deltaY_20260601.H`, computing the `neeta` suspension field), used to study particle entrainment and suspension in turbulent open-channel flow.

## OpenFOAM compatibility

Updated and verified to build and run against **OpenFOAM-2506** (this base CFDEMcoupling-PUBLIC repository was last updated by upstream for OpenFOAM-6).

## Building

With the CFDEM environment sourced (this repo's own `README.md` covers the full CFDEMcoupling-PUBLIC install):

```bash
cd $CFDEM_SOLVER_DIR/applications/solvers/PhasecfdemSolverPiso
wmake ./DPMTurbulenceModels   # build the custom turbulence library first
wmake .                        # build the solver
```

Developed as part of ongoing CFD-DEM research into bed load transport and bridge-pier local scour, built on OpenFOAM-2506 + CFDEMcoupling-PUBLIC.
