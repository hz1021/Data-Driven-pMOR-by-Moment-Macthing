# Data-Driven-pMOR-by-Moment-Macthing

# PMOR — Parametric Model Order Reduction Scripts

Compact MATLAB toolbox of scripts and examples for parametric model order
reduction (PMOR) used in research and experiment workflows.

## What this project does

- Provides multiple PMOR workflows and basis-construction methods implemented
  as MATLAB scripts (Loewner, BasisFunction, POD, Series Expansion,
  Matrix Interpolation, etc.).
- Includes example parameterized environments and precomputed matrices used
  for experiments and reproducible comparisons.

## Why this project is useful

- Collection of tested PMOR approaches useful for research, teaching, and
  rapid prototyping of reduced-order models for parameterized systems.
- Scripts automate dataset loading, basis construction, and error metrics
  computation (H2, L2, and related norms are included as `.mat` outputs).

## Key files and folders

- `mainPMORBasisFunction.m`, `mainPMORLoewner.m`, `mainPMORMatrixInterpolation.m` — primary entry scripts.
- `mainPMORPOD.m`, `mainPMORSeriesExpansion.m` — alternative reduction approaches.
- `environments/` — example parameterized model sets and helper scripts.
- `*.mat` — precomputed matrices (A0.mat, Ae.mat, B.mat, C.mat, etc.) used by examples.

## Requirements

- MATLAB (recommended recent release). The code uses standard MATLAB scripts and
  `.mat` data files; no external package manager is required.

## Quick start

1. Open MATLAB and set the project folder on the MATLAB path (point to the
   repository root).
2. Run one of the main entry scripts from the MATLAB command window. For
   example, to run the basis-function workflow:

```matlab
cd('<path-to-repo>')
addpath(genpath(pwd))
mainPMORBasisFunction
```

3. Try other main scripts to compare methods:

```matlab
mainPMORLoewner
mainPMORMatrixInterpolation
mainPMORPOD
mainPMORSeriesExpansion
```

4. Results and intermediate data are saved as `.mat` files in the project
   folder (see files like `H2NormBF.mat`, `l2NormSE_build_6.mat`, etc.).

## Examples and experiments

- The `environments/BuildingModel/` and
  `environments/SyntheticParametricModel/` folders contain example models,
  data preprocessing helpers, and sample experiment scripts (look for
  `main_*_build.m` and `main_modelBasedParametrizedPI.m`).

## Where to get help

- For questions or to report issues, open an issue in this repository.
- If this code was provided as part of a paper or project, check the
  associated publication or contact the project maintainer listed below.

## Maintainers & Contributing

- Maintainer: Please update this file with maintainer name and contact.
- Contributions: please open issues or pull requests. For larger
  contribution guidelines, add a `CONTRIBUTING.md` file and link it here.

## Notes & Next steps

- This repository is research-oriented. Consider adding the following for
  broader adoption:
  - `CONTRIBUTING.md` with coding and testing guidelines
  - Small example script that runs a minimal end-to-end demo
  - Unit or regression tests that validate core scripts

## License

Add a `LICENSE` file at the repository root and update this section with the
chosen license.

---
