"""Phase 4 common layer: frozen-implementation guard, one machine-readable record per simulation, JSONL storage.

The Phase 3 implementation (src/pbr, params/base_case.yaml, verification/tolerances.yaml) is used READ-ONLY.
meta() stops the study if any of their hashes differs from the frozen Phase 3 values.
"""
import hashlib
import json
import math
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402

from pbr import dimensionless, hotspot, numerics, params, reactor, sensitivity, shooting, storage  # noqa: E402
from pbr.guards import RunStatus  # noqa: E402

FROZEN = {"params_sha256": "6b4b35f2bd00ada4a46ee7b658e5ee9ce5db29654f58146b027074a5563d88ff",
          "code_sha256": "d6b1d5dfb9d392416567eb5c28592da6a17a47796ae21503aa16292b85448c8d",
          "tolerances_sha256": "8d58e4814e12b7f490c34fc0f9256829718c2bc686b8dcf4f2caa0c201055a22"}
OUT = ROOT / "results" / "phase4"
_META = {}


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


STUDY_MODULES = ("__init__.py", "common.py", "o2_variant.py", "studies.py")  # modules that PRODUCE data


def phase4_code_hash():
    """Hash of the data-producing Phase 4 modules (figures.py / tables.py only read stored records)."""
    h = hashlib.sha256()
    for f in sorted((ROOT / "phase4" / n) for n in STUDY_MODULES):
        h.update(f.name.encode())
        h.update(f.read_bytes())
    return h.hexdigest()


def meta():
    if not _META:
        actual = {"params_sha256": sha(ROOT / "params" / "base_case.yaml"), "code_sha256": storage.code_hash(),
                  "tolerances_sha256": sha(ROOT / "verification" / "tolerances.yaml")}
        bad = {k: (actual[k], v) for k, v in FROZEN.items() if actual[k] != v}
        if bad:
            raise RuntimeError(f"frozen Phase 3 implementation changed: {bad}")
        _META.update(actual, frozen_check="PASS", phase4_code_sha256=phase4_code_hash(),
                     environment=storage.environment())
    return _META


def short_meta():
    m = meta()
    return {k: m[k] for k in ("params_sha256", "code_sha256", "tolerances_sha256", "phase4_code_sha256")}


def classify(g):
    if g is None or not isinstance(g, float) or not math.isfinite(g):
        return "unclassified"
    return "amplified (g > 0)" if g > 0 else ("attenuated (g < 0)" if g < 0 else "boundary (g = 0)")


def run_point(p, study, label, mode="none", Pi=0.0, P0=None, config="physical", with_criterion=True,
              solver=None, rtol=None, atol_scaled=None, variant="base", keep_profile=False):
    """One simulation (+ optional criterion) -> one record with inputs, derived groups, solver settings, initial and
    final state, hot spot, conversion, cooling duty, pressure and sensitivity metrics, classification and hashes."""
    t_start = time.perf_counter()
    kw = dict(config=config, solver=solver, rtol=rtol, atol_scaled=atol_scaled)
    shoot = None
    if mode == "outlet":
        sr = shooting.solve_outlet(p, Pi, **kw)
        s = sr.solution
        shoot = {"status": sr.status.value, "detail": sr.detail, "P0": sr.P0, "residual": sr.residual,
                 "bracket": sr.bracket, "monotonic_scan": sr.monotonic, "sign_changes": sr.sign_changes,
                 "brent_iterations": sr.iterations}
        ok = sr.status == RunStatus.PASS and s is not None and s.status == RunStatus.PASS
        status = RunStatus.PASS.value if ok else (sr.status.value if sr.status != RunStatus.PASS else s.status.value)
    else:
        s = reactor.integrate(p, mode=mode, Pi=Pi, P0=P0, **kw)
        ok = s.status == RunStatus.PASS
        status = s.status.value
    rt, a_s = numerics.tolerances(p, rtol, atol_scaled)
    rec = {"study": study, "label": label, "variant": variant, "config": config, "mode": mode, "Pi": Pi,
           "P0_specified_Pa": P0, "inputs": {f: getattr(p, f) for f in params.FIELDS}, "overrides": dict(p.overrides),
           "derived": {**dimensionless.groups(p), "Pi": Pi, "dT_ad_K": p.dT_ad, "N_f_mol_s": p.N_f, "G_kg_m2s": p.G},
           "solver": solver or numerics.PRIMARY, "rtol": rt, "atol_scaled": a_s, "status": status,
           "detail": s.detail if s is not None else shoot["detail"],
           "guard_event": s.guard_event if s is not None else None, "shooting": shoot}
    if ok:
        hs = hotspot.detect(s, p)
        rec.update(initial_state=s.y[:, 0].tolist(), final_state=s.y[:, -1].tolist(), X_L=float(s.X[-1]),
                   T_L=float(s.T[-1]), T_max=hs["T_hs"], z_hs=hs["z_hs"], hs_kind=hs["kind"], dT_hs=hs["T_hs"] - p.Tc,
                   X_at_hs=float(reactor.evaluate(s, [hs["z_hs"]])[0][0]), Q_cool_W=float(s.Q[-1]),
                   P0_Pa=float(s.P[0]), PL_Pa=float(s.P[-1]), dP_Pa=float(s.P[0] - s.P[-1]),
                   P_mean_Pa=reactor.step_quadrature(s, lambda z, y: y[2]) / p.L, nfev=s.nfev)
        if keep_profile:
            idx = np.unique(np.linspace(0, s.z.size - 1, 301).astype(int))
            rec["profile"] = {"z": s.z[idx], "X": s.X[idx], "T": s.T[idx], "P": s.P[idx]}
    rec["boundary_class"] = "unclassified"
    if with_criterion and ok:
        c = sensitivity.criterion(p, mode=mode, Pi=Pi, P0=P0 if mode == "constant" else None, **kw)
        rec["criterion_status"] = c["status"].value
        if c["status"] == RunStatus.PASS:
            dense = c["res"]["dense"]
            dTdT0 = float(dense([rec["z_hs"]])[5][0]) if rec["hs_kind"] == "interior" else None
            rec.update(S_int=c["S_int"], g=c["g"], z_Sint=c["z_Sint"], S_T_max_profile=float(c["res"]["S"][0][1].max()),
                       dS0=c["dS0"], dP0dT0=c["dP0dT0"], criterion_rule=c["rule"], dTmax_dT0_forward=dTdT0,
                       boundary_class=classify(float(c["g"])))
    rec["wall_s"] = time.perf_counter() - t_start
    rec["meta"] = short_meta()
    return rec


ROOT_KEYS = ("detail", "root", "U_crit", "T0_crit", "xtol", "bracket", "coarse_bracket", "iterations", "function_calls",
             "converged", "g_at_root", "sign_changes", "subscan_sign_changes", "local_monotonic_at_root",
             "globally_monotonic", "scan_x", "scan_g", "at_root", "local_g")


def root_summary(r, study, label, kind, **extra):
    """kind = 'U' (g decreasing in U) or 'T0' (g increasing in T0). ROOT_FAILURE with a single-signed scan is
    reported as a boundary outside the searched range (never silently dropped)."""
    out = {"study": study, "label": label, "kind": kind, **extra, **{k: r[k] for k in ROOT_KEYS if k in r}}
    st = r["status"]
    out["status"] = st.value if hasattr(st, "value") else st
    out["outcome"] = "root found" if out["status"] == "PASS" else "root failure"
    g = np.asarray(r.get("scan_g", []), float)
    if out["status"] != "PASS" and g.size and np.all(np.isfinite(g)):
        if np.all(g > 0):
            out["outcome"] = ("g > 0 at every scanned U: U_crit above the searched range" if kind == "U"
                              else "g > 0 at every scanned T0: T0_crit below the searched range")
        elif np.all(g < 0):
            out["outcome"] = ("g < 0 at every scanned U: U_crit below the searched range" if kind == "U"
                              else "g < 0 at every scanned T0: T0_crit above the searched range")
    out["meta"] = short_meta()
    return out


def write_jsonl(name, recs):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / f"{name}.jsonl", "w") as f:
        for r in recs:
            f.write(json.dumps(storage.jsonable(r), allow_nan=False) + "\n")


def write_json(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(storage.jsonable(obj), indent=1, allow_nan=False))


def read_jsonl(name):
    return [json.loads(line) for line in (OUT / f"{name}.jsonl").read_text().splitlines() if line.strip()]


def read_json(name):
    return json.loads((OUT / f"{name}.json").read_text())
