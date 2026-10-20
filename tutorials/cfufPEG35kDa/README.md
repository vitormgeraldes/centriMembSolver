# cfufPEG35kDa

Thirty-six axisymmetric (wedge) simulations of concentration polarization of
PEG 35 kDa in the radial cross-flow ultrafiltration cell (CFUF), run with
`centriMembSolver` from this repository. Case `Q<Q>_C<C0>_P<TMP>`: feed flow rate
[L/min × 100], feed concentration [kg/m³] and transmembrane pressure [bar × 10];
for example `Q282_C10_P12` is 2.82 L/min, 10 kg/m³, 1.2 bar. The operating
conditions of each case are in `constant/simulationConditions`.

## Requirements

* OpenFOAM v2512 (https://www.openfoam.com) with MPI.
* About 2 GB of RAM per case.

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
cd tutorials/cfufPEG35kDa/Q282_C10_P12
./Allrun
```

`Allrun` generates the mesh (`makeMesh`: `blockMesh`, `surfaceFeatureExtract`,
`snappyHexMesh`, `createPatch`, `extrudeMesh`, nine passes of `refineWallLayer`
at the membrane, `renumberMesh`; a 1.5° wedge with 64 940 cells), decomposes the
case, runs `centriMembSolver` in parallel and reconstructs the last time.
`./Allclean` brings the case back to its initial state.

## Notes

* Run OpenFOAM commands from inside the case folder: the dictionaries
  `#include "./constant/simulationConditions"`, which is resolved from the working
  directory.
* The number of MPI ranks is `numberOfSubdomains` in `system/decomposeParDict`; it
  can be changed.
* `makeMesh` takes about 5 minutes (`snappyHexMesh` runs in serial).
* `endTime` is 30 s and fields are written every 0.01 s with `purgeWrite 0`: about
  3000 time directories of about 3.7 MB each (~11 GB per case). Set `purgeWrite N` in
  `system/controlDict` to keep only the last N.
* Quick test: with `endTime 0.03` the solver runs about 1900 (Q282) to 3700
  (Q564) time steps.
* The solver log contains many warnings `Wedge patch 'sym1' is not planar` (and
  `sym2`), with deviations of order 1e-13 from the rounding of the point
  coordinates; they are harmless.
* Monitoring output (Γ on the membrane, solute balance, intrinsic rejection, wall
  shear stress) is written to `postProcessing/`.

## Reference

José Pedro Bernardo da Costa Pereira, *Scaled-Down Screening of Ultrafiltration
Membranes with Predictable Concentration Polarization: The Ortho-Centrifugal Membrane
Filter against Cross-Flow*, MSc thesis, Instituto Superior Técnico, Universidade de
Lisboa, 2026.
