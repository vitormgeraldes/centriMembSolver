# BSA Ultrafiltration Case Download Bundle

This folder is a download-ready bundle of the BSA ultrafiltration tutorial case.

## Included case

- `BSA_nG200_C1_P20_J20_Lx5_Ly1_Lz25/`

## How to run

```bash
cd BSA_nG200_C1_P20_J20_Lx5_Ly1_Lz25
export solute=BSA
chmod +x makeMesh
./makeMesh
centriMembSolver
```

## Notes

- This case expects OpenFOAM with `centriMembSolver` available.
- `transportProperties` selects species data through `${solute}`, so keep `solute=BSA`.
