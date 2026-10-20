# ocmfPEG35kDa

Twelve 3-D DNS cases of concentration polarization of PEG 35 kDa in the
Ortho-Centrifugal Membrane Filter (O-CMF), run with `centriMembSolver` from this
repository. Case `R<rpm>_C<C0>`: rotor speed [rpm] and feed concentration [kg/m³];
the operating conditions of each case are in `constant/simulationConditions`.

## Requirements

* OpenFOAM v2512 (https://www.openfoam.com) with MPI.
* About 6 GB of RAM per case.

## 1. Build the solver

From the root of this repository:

```bash
source <OpenFOAM-v2512>/etc/bashrc
wmake src/membraneBCs
wmake src/transportModels
wmake src/functionObjects/GamaPatchAvgTimeWindow
wmake solver
```

The `Make/files` install the executable and libraries in `$FOAM_SITE_APPBIN` and
`$FOAM_SITE_LIBBIN`. To install them in your user directories instead, replace
`FOAM_SITE_` by `FOAM_USER_` in the four `Make/files` and in `solver/Make/options`.

## 2. Run a case

```bash
cd tutorials/ocmfPEG35kDa/R7070_C10
./Allrun
```

`Allrun` generates the mesh (`makeMesh`: `blockMesh`, four passes of
`refineWallLayer`, `renumberMesh`; 2 467 584 cells), decomposes the case, runs
`centriMembSolver` in parallel and reconstructs the last time. `./Allclean` brings the
case back to its initial state.

## Notes

* Run OpenFOAM commands from inside the case folder: the dictionaries
  `#include "./constant/simulationConditions"`, which is resolved from the working
  directory.
* The number of MPI ranks is `numberOfSubdomains` in `system/decomposeParDict`; it can
  be changed.
* `endTime` is 30 s and fields are written every 0.01 s with `purgeWrite 0`: about
  3000 time directories of about 290 MB each (~870 GB per case). Set `purgeWrite N` in
  `system/controlDict` to keep only the last N.
* Quick test: with `endTime 0.03` the solver runs about 15 time steps.
* Monitoring output (Γ on the membrane, solute balance, intrinsic rejection) is
  written to `postProcessing/`.

## Reference

José Pedro Bernardo da Costa Pereira, *Scaled-Down Screening of Ultrafiltration
Membranes with Predictable Concentration Polarization: The Ortho-Centrifugal Membrane
Filter against Cross-Flow*, MSc thesis, Instituto Superior Técnico, Universidade de
Lisboa, 2026.
