# Extension-study code

**Author: Aryan Saxena**

This directory contains only model implementation and numerical-analysis scripts. It intentionally excludes generated numerical outputs, manuscript-building scripts, Scopus/literature-screening scripts, and research documents.

## Contents

- `pbr.py`: generalized pressure-order model and sensitivity/root calculations.
- `pbr2.py`: second parameter dataset.
- `pbr3.py`: three- and five-reaction network implementations (imports `pbr2.py`).
- `sweep_n.py`, `coolant_pressure.py`, `kscan_ext.py`, `assemble.py`, `figs.py`: generalized pressure-order and coolant studies.
- `ds2_*.py`: second-dataset calculations and source-comparison checks.
- `nets_*.py`, `cc_*.py`, `delasa_*.py`: multi-step network runs and parameter checks.
- `pi121_*.py`: historical parametric comparison calculations. These values are not the adopted reference-state pressure parameter for the 2 mm source-based case.
- `validate_fast.py`, `verify_*.py`, `re_benchmark_sensitivity.py`: numerical verification and benchmark checks.

## Environment

Use the repository-root Python environment installed from `requirements.txt` (Python 3.12; NumPy, SciPy, PyYAML, and Matplotlib are pinned there).

## Recommended runs

Run commands from this directory:

```bash
python sweep_n.py
python coolant_pressure.py
python kscan_ext.py
python assemble.py
python figs.py
python ds2_run.py
python ds2_rhoL_run.py
python nets_run_fast.py delasa
python nets_run_fast.py cc
python nets_rhoL_run.py
python cc_scale_run.py
```

Additional targeted checks:

```bash
python delasa_dH22_check.py
python delasa_heats_check.py
python ds2_critfeed_check.py
python ds2_jiang_check.py
python ds2_vwf_fig3_check.py
python validate_fast.py
python verify_T0_scan.py
python verify_table7.py
python re_benchmark_sensitivity.py
```

Some network calculations are computationally intensive and may take several minutes. The scripts write result JSON and logs into this directory; generated files are not committed in this pre-publication release. Inspect each script's command-line interface before running a targeted check. For `pi121_run.py`, use one model argument at a time: `ds2`, `delasa`, or `cc`.

## Interpretation cautions

- The five-reaction network's reactor-basis loading factor remains unverified; the corresponding checks are exploratory and do not establish validation of that factor.
- `Pi = 0.011` and `Pi = 0.121` are retained only as historical parametric comparison points. The source-based second-dataset case adopted in the manuscript is 2 mm particles and approximately `Pi = 0.20` at the specified model reference state using the inertia-dominated Ergun term.
- The model has no experimental validation. Do not interpret its extrapolative high-temperature results as physical plant predictions.
