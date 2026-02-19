# Data-Driven-pMOR-by-Moment-Macthing

This project provides compact MATLAB scripts and examples for the novel parametric model order reduction (pMOR) framework established in the arXiv paper "[Data-Driven Model Reduction by Moment Matching for Linear and Nonlinear Parametric Systems](https://arxiv.org/abs/2506.10866)".

## What this project does

- Provides the novel pMOR workflows in the arXiv paper, *i.e.*, series expansion and basis function approaches.
- Benchmarks some existing and well-established pMOR workflows in the literature, *i.e.*, multi-parameter moment matching, global Proper Orthognal Decomposition (POD), stability-preserving matrix interpolation, and parametrised Loewner framework.
- Includes precomputed system matrices used for experiments and reproducible comparisons.

## Note for key files

- `main*.m`, *e.g.*, `mainPMORBasisFunction.m` and `mainPMORSeriesExpansion.m` are primary entry scripts implementing aforementioned pMOR methods.
- `H_norm_error_SEBF.m` is the script for computing the relative norm of the error system, *i.e.*, the decrepency between the original full-order parametric model and the resulting reduced-order parametric model.
- The other `*.m` are common or specific funtions called in `main*.m`.
- `*.mat` are precomputed system matrices used by examples.

## Requirements

- MATLAB (recommended recent release). The code uses standard MATLAB scripts and
  `.mat` data files; Control System Toolbox may be required.

## Quick start

1. Open MATLAB and set the project folder on the MATLAB path, *i.e.*, point to the
   repository root.
2. Run one of the main entry scripts from the MATLAB command window. For
   example, one can apply pMOR via the basis function approach using the following commands.

```matlab
current_folder = pwd
addpath(genpath(pwd))
mainPMORBasisFunction
```

3. Apply several pMOR techniques to the given system and evaluate their performance by computing the relative $\mathcal{H}_2$-norm error between each reduced-order model and the full-order model. Users can select different pMOR techniques `*` by uncommenting the corresponding `main*` call at the beginning of `H_norm_error_SEBF.m`.

```matlab
H_norm_error_SEBF
```

4. Results and intermediate data can be saved as `.mat` files in the project
   folder, depending on user requirements.


## Maintainers & Contributing

- Maintainer: Hanqing Zhang and Junyu Mao - [MAC-X Lab](https://giordanoscarciotti.com/mac-x-lab/).
- Contributions: please open issues or pull requests. For larger
  contribution guidelines, add a `CONTRIBUTING.md` file and link it here.


## License

This project is licensed under the MIT License.

---
