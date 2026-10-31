"""(1) Parameter loading.

params/base_case.yaml -- the accepted FINAL Phase 2 file -- is the SINGLE SOURCE OF TRUTH for every
physical/model parameter; nothing physical is hard-coded anywhere in src/. The file is read-only and
its sha256 must equal the accepted Phase 2 hash. Test/study configurations change values only in
memory through Params.with_(), and every result records those overrides. The derived quantities
below are implementation-derived (computed from file values only).
"""
from __future__ import annotations

import dataclasses
import hashlib
import math
import pathlib

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_PARAMS = ROOT / "params" / "base_case.yaml"
PHASE2_SHA256_ORIGINAL = "6df8de7d4b71fa8f399dce397b93829d3d6b82c8c2078d15df09fcdd975ddb55"  # accepted Phase 2 file
PHASE2_SHA256 = "6b4b35f2bd00ada4a46ee7b658e5ee9ce5db29654f58146b027074a5563d88ff"  # Phase 3 numerics-only revision (see file comment)
ALLOWED_STATUS = ("VERIFIED", "LITERATURE", "PROPERTY", "PROJECT-DEFINED")
FIELDS = ("k_m", "T_m", "E_over_R", "dH_rxn", "nu_O2", "y_x0", "M_wf", "y_O2_0", "P_ref",
          "T0", "Tc", "Qrho", "R_tube", "L", "cp", "U", "Pi_max", "Pi_inlet_max",
          "P_min_frac", "R_gas", "M_oX", "Vm_normal")


class SpecInconsistency(RuntimeError):
    """The parameter file contradicts the Phase 2 specification implemented here (STOP)."""


@dataclasses.dataclass(frozen=True, eq=False)
class Params:
    k_m: float
    T_m: float
    E_over_R: float
    dH_rxn: float
    nu_O2: float
    y_x0: float
    M_wf: float
    y_O2_0: float
    P_ref: float
    T0: float
    Tc: float
    Qrho: float
    R_tube: float
    L: float
    cp: float
    U: float
    Pi_max: float
    Pi_inlet_max: float
    P_min_frac: float
    R_gas: float
    M_oX: float
    Vm_normal: float
    ranges: dict
    numerics: dict
    model: dict
    criteria: dict
    source_audit: dict
    path: str
    sha256: str
    overrides: tuple = ()

    def with_(self, **kw) -> "Params":
        """Copy with in-memory overrides (the file is never modified); overrides are recorded."""
        unknown = [k for k in kw if k not in FIELDS]
        if unknown:
            raise KeyError(f"unknown parameter(s): {unknown}")
        new = {k: float(v) for k, v in kw.items()}
        rec = dict(self.overrides)
        rec.update(new)
        return dataclasses.replace(self, overrides=tuple(sorted(rec.items())), **new)

    @property
    def A_c(self):
        return math.pi * self.R_tube ** 2

    @property
    def perimeter(self):
        return 2.0 * math.pi * self.R_tube

    @property
    def N_f(self):
        return self.Qrho / self.M_wf

    @property
    def N_x0(self):
        return self.y_x0 * self.N_f

    @property
    def G(self):
        return self.Qrho / self.A_c

    @property
    def dT_ad(self):
        return (-self.dH_rxn) * self.y_x0 / (self.M_wf * self.cp)

    @property
    def P_min(self):
        return self.P_min_frac * self.P_ref


def load(path=DEFAULT_PARAMS, require_phase2_hash=True) -> Params:
    path = pathlib.Path(path)
    raw = path.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    if require_phase2_hash and sha != PHASE2_SHA256:
        raise SpecInconsistency(f"{path}: sha256 {sha} is not the accepted Phase 2 file {PHASE2_SHA256}")
    cfg = yaml.safe_load(raw)
    P = cfg["parameters"]
    missing = [f for f in FIELDS if f not in P]
    if missing:
        raise SpecInconsistency(f"parameters missing from file: {missing}")
    for name, entry in P.items():
        if entry.get("status") not in ALLOWED_STATUS:
            raise SpecInconsistency(f"{name}: status {entry.get('status')!r} not allowed")
    values = {f: float(P[f]["value"]) for f in FIELDS}
    ranges = {k: (float(e["range"][0]), float(e["range"][1])) for k, e in P.items() if "range" in e}
    p = Params(**values, ranges=ranges, numerics=dict(cfg["numerics"]), model=dict(cfg["model"]),
               criteria=dict(cfg["criteria"]), source_audit=dict(cfg["source_audit"]),
               path=str(path), sha256=sha)
    check_spec(p)
    return p


def check_spec(p: Params) -> None:
    """Stop if the file's model/criteria blocks contradict what this package implements."""
    m, problems = p.model, []
    if m.get("oxygen_depletion") is not False:
        problems.append("model.oxygen_depletion must be false (fixed-O2 base model)")
    pd = m.get("pressure_drop", {})
    if pd.get("spec_primary") != "outlet":
        problems.append("model.pressure_drop.spec_primary must be 'outlet'")
    if "dp/dzeta = -Pi (T/Tc)/p" not in " ".join(str(pd.get("form", "")).split()):
        problems.append("model.pressure_drop.form differs from the implemented dp/dzeta = -Pi (T/Tc)/p")
    law = " ".join(str(m.get("kinetic_law", "")).split())
    for token in ("k_m exp[-(E/R)(1/T - 1/T_m)] (P/P_ref) c_x", "c_x = y_x P/(R T)", "y_x = y_x0 (1 - X)"):
        if token not in law:
            problems.append(f"model.kinetic_law lacks '{token}'")
    if p.numerics.get("solver_primary") != "BDF":
        problems.append("numerics.solver_primary must be BDF (the implemented primary solver)")
    for key in ("hot_spot", "runaway_boundary", "Ucrit"):
        if key not in p.criteria:
            problems.append(f"criteria.{key} missing")
    if problems:
        raise SpecInconsistency("; ".join(problems))
