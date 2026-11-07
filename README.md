# o-Xylene Packed-Bed Reactor Runaway Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23852766.svg)](https://doi.org/10.5281/zenodo.23852766)

**Author:** Aryan Saxena
**Affiliation:** School of Chemical Engineering, Vellore Institute of Technology, Vellore, India

This repository contains a one-dimensional pseudo-homogeneous packed-bed reactor model for studying thermal behaviour during o-xylene oxidation. It includes the model, frozen base-case parameters, numerical verification routines and extension-study scripts.

## Contents

- `src/pbr/` — reactor model and supporting calculations.
- `params/` — frozen base-case parameters.
- `verification/` — numerical tests, reference calculations and tolerances.
- `phase4/` — pressure-decomposition and oxygen-model studies.
- `extension_study/` — additional pressure-order, coolant, parameter-dataset and reaction-network studies.
- `checks/` and `runs/` — parameter checks and base-case run scripts.

The manuscript, supplementary information, literature-screening records, Scopus exports, publisher PDFs and full derived result datasets are not included in this software release.

## Requirements

- Python 3.12
- Packages listed in `requirements.txt`

Create an environment and install the dependencies:

```bash
python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1
# macOS/Linux
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Run the base verification

From the repository root:

```bash
bash run_all.sh
```

On Windows without Bash, run these commands in PowerShell:

```powershell
python checks/phase2_plausibility.py params/base_case.yaml
python verification/run_verification.py
python verification/make_report.py
python runs/base_case_sanity.py
```

The base verification sequence performs the parameter plausibility check, the A1/V1–V13 verification suite, report generation, and the base-case sanity run. Generated outputs are written locally under `results/` and are not included in this release.

For the additional studies, see [`extension_study/README.md`](extension_study/README.md) for the run order and commands. Scripts write generated files locally; those outputs are not bundled with this release.

## Model scope

The model represents a steady-state, one-dimensional pseudo-homogeneous plug-flow packed-bed reactor with wall cooling. The runaway boundary is based on a steady-state parametric-sensitivity criterion. It is not a dynamic-stability, saddle-node or bifurcation criterion. The model has not been experimentally validated and should not be used as a plant-safety assessment.

## Citation

Please cite version 1.0.0 using its [Zenodo version DOI](https://doi.org/10.5281/zenodo.23852766). The [concept DOI](https://doi.org/10.5281/zenodo.23852765) resolves to the record covering all versions.

## Licence

Copyright © 2026 Aryan Saxena. All rights reserved. No open-source licence is granted. Third-party publications and materials remain subject to their respective rights and licensing terms.
