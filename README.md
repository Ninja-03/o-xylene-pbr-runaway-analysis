# o-Xylene packed-bed reactor runaway analysis

**Author: Aryan Saxena**  
School of Chemical Engineering, Vellore Institute of Technology, Vellore, India

This repository contains the model implementation, frozen parameter file, numerical verification code, and extension-study scripts associated with the research project:

> Separating pressure-level and axial pressure-decline effects on a sensitivity-based parametric runaway boundary in a model o-xylene packed-bed reactor.

## Scope of this public software release

Included:
- the one-dimensional pseudo-homogeneous packed-bed reactor implementation (`src/pbr/`);
- the frozen base parameter file (`params/base_case.yaml`);
- numerical verification code and tolerances (`verification/`);
- pressure-decomposition and oxygen-model extension code (`phase4/`);
- extension-study scripts for generalized pressure order, a second parameter dataset, and multi-step reaction networks (`extension_study/`).

Not included in this pre-publication software release:
- the manuscript, supplementary information, or unpublished manuscript drafts;
- Scopus search exports, forward-citation CSVs, screening worksheets, abstracts, or literature-screening reports;
- downloaded publisher PDFs or other copyrighted source papers;
- internal review notes, intermediate document-generation scripts, failed-attempt logs, or private working files;
- manuscript tables/figures and the full derived result datasets. The included scripts generate numerical outputs locally when run.

These exclusions are deliberate. The complete internal research archive is maintained separately and is not this public release package.

## Model scope

The base implementation models o-xylene oxidation in a one-dimensional pseudo-homogeneous plug-flow packed-bed reactor with wall cooling. The pressure-drop parameterisation and the sensitivity-based parametric runaway boundary are documented in the source code and verification suite. The criterion is a steady-state parametric-sensitivity criterion; it is not a dynamic-stability, saddle-node, or bifurcation criterion.

The `phase4/` directory contains the controlled pressure-decomposition and oxygen-model study implementation. `extension_study/` contains separate implementations and scripts for pressure-order, co-current-coolant, second-dataset, and multi-step-network calculations.

## Requirements

- Python 3.12
- Dependencies pinned in `requirements.txt`

Create an environment and install dependencies:

```bash
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# macOS/Linux:
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Run the base verification workflow

From the repository root:

```bash
bash run_all.sh
```

The script runs the parameter plausibility check, the A1/V1–V13 verification suite, and the base-case sanity run. Generated outputs are written under `results/` and are not included in this pre-publication release.

On Windows, run the commands in `run_all.sh` individually in PowerShell if Bash is unavailable:

```powershell
python checks/phase2_plausibility.py params/base_case.yaml
python verification/run_verification.py
python verification/make_report.py
python runs/base_case_sanity.py
```

## Extension-study scripts

See [`extension_study/README.md`](extension_study/README.md) for the recommended run order and model-specific commands. The extension scripts write JSON, log, and figure outputs locally; those generated outputs are excluded from this pre-publication release.

## Reproducibility

- The base parameter file and verification tolerances are frozen inputs.
- Verification outputs contain code and parameter hashes when generated.
- Do not edit the frozen parameter file without updating the corresponding verification records and explaining the change.
- Results from the extension scripts should be regenerated from the included scripts rather than treated as bundled result data.

## Citation

A `CITATION.cff` file is provided for this software release. Once a versioned GitHub release has been archived by Zenodo, cite the version-specific Zenodo DOI in research outputs.

## Rights and reuse

Copyright © 2026 Aryan Saxena. All rights reserved.

No open-source licence is granted for this software. Users may view, use, display, and fork this public repository through GitHub as permitted by GitHub's Terms of Service. No additional licence is granted for independent reuse, modification, redistribution, or incorporation of the software into other projects. Any permissions beyond those provided through GitHub's platform functionality require separate permission from the copyright holder.

Third-party publications are cited in the parameter file and source code but are not redistributed in this repository. Third-party materials remain subject to their respective rights and licensing terms.
