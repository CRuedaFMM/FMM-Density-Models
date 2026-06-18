# FMM Density Models for Circular and Toroidal Data

Code accompanying the paper:

**FMM Density Models for Circular and Toroidal Data: Applications to Ramachandran Geometry**

The repository is organized in two parts:

- `R/`: reusable functions.
- `scripts/`: reproducible scripts for the analyses in the paper.

The core S-FMM fitting function is assumed to be available from the FMM project codebase:

https://www.eio.uva.es/the-fmm-project/

In particular, the toroidal scripts expect a function compatible with:

```r
S_fmm(M, n_iter, K, L)
```

returning at least:

```r
list(r2 = ..., par_current = ...)
```

## Main scripts

Run from the repository root.

```r
source("scripts/run_circular_benchmarks.R")
source("scripts/run_toroidal_protein_simulations.R")
source("scripts/run_toroidal_synthetic_simulations.R")
source("scripts/make_tables_figures.R")
```

## Contents

### Circular analyses

- `R/circular_density_models.R`: circular densities, likelihoods, mixtures, and Fourier2 density.
- `R/circular_benchmark_models.R`: model fitting and circular benchmark comparison functions.
- `scripts/run_circular_benchmarks.R`: template script for the real circular datasets used in the paper.

### Toroidal analyses

- `R/toroidal_utils.R`: grids, normalization, circular distances, sampling, metrics.
- `R/toroidal_pdb.R`: construction of Ramachandran KDE surfaces from PDB structures.
- `R/toroidal_perturbations.R`: noise-based and structured toroidal perturbations.
- `R/toroidal_fourier33.R`: truncated two-dimensional Fourier representation with `K=L=3`.
- `R/toroidal_vm4.R`: four-component mixture of bivariate product von Mises distributions (VM4).
- `R/toroidal_simulation_pipeline.R`: fitting loop for Fourier33, S-FMM2, and VM4.
- `scripts/run_toroidal_protein_simulations.R`: protein-based simulations for 1CAG and 1V9E.
- `scripts/run_toroidal_synthetic_simulations.R`: synthetic toroidal simulations.

## Notes

The scripts are designed to avoid local user-specific paths. Outputs are written to `outputs/`.
