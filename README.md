# FMM-Density-Models
Software accompanying the paper "FMM Density Models for Circular and Toroidal Data: Applications to Ramachandran Geometry". The repository contains implementations of FMM-based density models for circular and toroidal data, simulation studies, model fitting procedures, and code to reproduce the analyses and figures presented in the manuscript.

## Repository structure

### `R/`

Contains the reusable functions used throughout the analyses. These files
implement the circular and toroidal models, fitting procedures, perturbation
mechanisms, simulation utilities, performance measures, and other supporting
functions required by the main analyses.

### `scripts/`

Contains the main analysis programs used to reproduce the results reported in
the manuscript and Supplementary Material. These scripts call the functions
defined in `R/` and run the circular simulations, toroidal simulations,
sensitivity analyses, Ramachandran applications, and the generation of tables
and figures.
