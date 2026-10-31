#!/usr/bin/env bash
# Reproduce Phase 3 from a clean checkout:  pip install -r requirements.txt && bash run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
python3 checks/phase2_plausibility.py params/base_case.yaml
python3 verification/run_verification.py          # A1, V1-V13 -> results/verification/*.json
python3 verification/make_report.py               # -> results/verification_report.{json,md}
python3 runs/base_case_sanity.py                  # -> results/sanity/base_case_sanity.{json,png}
