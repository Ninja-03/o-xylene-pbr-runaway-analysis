"""(14) Verification tests V1-V13 and provenance audit A1 for the Phase 3 implementation.

Each test returns {id, purpose, expected, actual, error, tolerance, pass, comments, details}.
Acceptance tolerances come from verification/tolerances.yaml (pre-registered before the first run).
Test configurations change parameters only in memory (Params.with_); the overrides are listed.
"""
import ast
import hashlib
import math
import pathlib
import sys

import numpy as np
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pbr import (dimensionless, energy, guards, hotspot, kinetics, numerics, params,  # noqa: E402
                 pressure, reactor, sensitivity, shooting)
from pbr.guards import GuardError, RunStatus  # noqa: E402
from verification import reference_re_code as RE  # noqa: E402

TOL_PATH = ROOT / "verification" / "tolerances.yaml"
TOL = yaml.safe_load(TOL_PATH.read_text())
PASS = RunStatus.PASS
ALL = ["A1", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"]
REJECTED_V01 = (19.837, 1300.0, 0.0254, 82.7, 96.18, 0.237, 4684.0, 713.15, 900.0, 0.95, 2500.0, 29.48,
                1922.6, -1361.0)


PWL = "PASS WITH DOCUMENTED LIMITATION"


def record(tid, purpose, expected, actual, error, tolerance, ok, comments, details=None, disposition=None):
    return {"id": tid, "purpose": purpose, "expected": expected, "actual": actual, "error": error,
            "tolerance": tolerance, "pass": bool(ok), "disposition": disposition or ("PASS" if ok else "FAIL"),
            "comments": comments, "details": details or {}}


def metric_tols(p):
    return {"X_L": TOL["X_abs"], "T_L": TOL["T_abs_K"], "T_max": TOL["T_abs_K"], "z_hs": TOL["z_rel_L"] * p.L,
            "P0": TOL["P_rel"] * p.P_ref, "P_L": TOL["P_rel"] * p.P_ref, "Q_cool": TOL["Q_rel"] * p.Qrho * p.cp * p.Tc}


def compare(ref, other, p):
    rows, ok = {}, True
    for k, tol in metric_tols(p).items():
        d = abs(other[k] - ref[k])
        rel = d / abs(ref[k]) if ref[k] != 0 else (0.0 if d == 0 else float("inf"))
        rows[k] = {"ref": ref[k], "other": other[k], "abs_diff": d, "rel_diff": rel, "tol_abs": tol, "pass": d <= tol}
        ok = ok and d <= tol
    return rows, ok


def worst(rows):
    k = max(rows, key=lambda k: rows[k]["abs_diff"] / rows[k]["tol_abs"])
    return k, rows[k]["abs_diff"], rows[k]["tol_abs"]


def solve(p, mode="none", Pi=0.0, **kw):
    if mode == "outlet":
        sr = shooting.solve_outlet(p, Pi, **kw)
        return sr.solution, sr
    return reactor.integrate(p, mode=mode, Pi=Pi, **kw), None


def good(s, sr=None):
    return s is not None and s.status == PASS and (sr is None or sr.status == PASS)


def CASES(p):
    return [("baseline Pi = 0", "none", 0.0), ("outlet P(L) = P_ref, Pi = Pi_max", "outlet", p.Pi_max)]


# ------------------------------------------------------------------------------------------------ A1
def A1(p):
    purpose = "Provenance audit: accepted Phase 2 file in use; no parameter-file or rejected value hard-coded in src/"
    sha = hashlib.sha256(pathlib.Path(p.path).read_bytes()).hexdigest()
    cfg = yaml.safe_load(pathlib.Path(p.path).read_text())
    file_vals = set()
    for e in cfg["parameters"].values():
        file_vals.add(float(e["value"]))
        file_vals.update(float(v) for v in e.get("range", []))
    file_vals = {v for v in file_vals if v not in (0.0, 1.0, 2.0)}
    hits, n_lit = [], 0
    for f in sorted((ROOT / "src" / "pbr").glob("*.py")):
        text = f.read_text()
        for node in ast.walk(ast.parse(text)):
            if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
                if isinstance(node.value, int) and abs(node.value) < 10:
                    continue
                n_lit += 1
                v = float(node.value)
                for ref in list(file_vals) + list(REJECTED_V01):
                    if v != 0.0 and abs(v - ref) <= 1e-9 * abs(ref):
                        kind = "parameter-file value" if ref in file_vals else "REJECTED v0.1 value"
                        hits.append(f"{f.name}:{node.lineno} literal {node.value!r} equals {kind} {ref}")
        if "v0.1" in text or "archive" in text:
            hits.append(f"{f.name} references archived v0.1 material")
    rej_in_params = sorted({v for v in REJECTED_V01 for w in file_vals if abs(v - w) <= 1e-9 * abs(v)})
    audit_only = sorted(k for k in cfg.get("source_audit", {}))
    ok = sha == params.PHASE2_SHA256 and not hits and not rej_in_params
    return record("A1", purpose, "sha256 = accepted Phase 2 hash; 0 hard-coded parameter/rejected literals in src/",
                  f"sha256 {sha[:12]}...; {n_lit} numeric literals scanned, {len(hits)} matches; rejected values in "
                  f"parameters block: {rej_in_params or 'none'}", f"{len(hits)} matches", "0 matches", ok,
                  "source_audit entries (e.g. ZB2020 Table 3 literals, R&E Table 6.5 misprints) are audit records "
                  "only: the model never reads them; V9 reads only re_documented_runaway_feed_T_K.",
                  {"sha256": sha, "matches": hits, "source_audit_keys": audit_only,
                   "tolerances_sha256": hashlib.sha256(TOL_PATH.read_bytes()).hexdigest()})


# ------------------------------------------------------------------------------------------------ V1
def V1(p):
    purpose = ("Isothermal, no pressure drop (heat effects disabled: dH = 0, U = 0): X(z) vs analytic "
               "1 - exp(-kappa z), kappa = A_c k(T0) P_ref/(N_f R T0); P^2 dependence at constant P = 1.3 P_ref")
    rows, ok = [], True
    for T0, fac in ((600.0, 1.0), (625.0, 1.0), (650.0, 1.0), (625.0, 1.3)):
        q = p.with_(dH_rxn=0.0, U=0.0, T0=T0)
        s = reactor.integrate(q, mode="none") if fac == 1.0 else reactor.integrate(q, mode="constant", P0=fac * q.P_ref)
        k = q.k_m * math.exp(-q.E_over_R * (1.0 / T0 - 1.0 / q.T_m))
        kappa = q.A_c * k * q.P_ref / (q.N_f * q.R_gas * T0) * fac ** 2
        Xa = 1.0 - np.exp(-kappa * s.z)
        e = np.abs(s.X - Xa)
        row = {"case": f"T0 = {T0:g} K, P = {fac:g} P_ref", "status": s.status.value,
               "X_L_analytic": float(Xa[-1]), "X_L_numerical": float(s.X[-1]),
               "abs_err_L": float(abs(s.X[-1] - Xa[-1])), "rel_err_L": float(abs(s.X[-1] - Xa[-1]) / Xa[-1]),
               "max_abs_err_profile": float(e.max()), "T_drift_K": float(np.max(np.abs(s.T - T0)))}
        if fac != 1.0:
            ratio = kinetics.rate(0.0, T0, fac * q.P_ref, q) / kinetics.rate(0.0, T0, q.P_ref, q)
            row.update(rate_ratio=float(ratio), rate_ratio_expected=fac ** 2)
            ok = ok and abs(ratio - fac ** 2) <= 1e-12 * fac ** 2
        ok = ok and s.ok and row["max_abs_err_profile"] <= TOL["X_abs"] and row["T_drift_K"] <= TOL["T_abs_K"]
        rows.append(row)
    w = max(rows, key=lambda r: r["max_abs_err_profile"])
    return record("V1", purpose, f"|X_num - X_ana| <= {TOL['X_abs']:g} on every profile; rate ratio = (P/P_ref)^2",
                  "X(L) num/ana: " + "; ".join(f"{r['X_L_numerical']:.9f}/{r['X_L_analytic']:.9f}" for r in rows),
                  f"max abs {w['max_abs_err_profile']:.2e}; max rel(L) {max(r['rel_err_L'] for r in rows):.2e}",
                  f"abs {TOL['X_abs']:g}", ok,
                  f"rate ratio at 1.3 P_ref = {rows[-1]['rate_ratio']:.15f} vs 1.69: implemented rate scales as P^2",
                  {"rows": rows, "overrides": "dH_rxn = 0, U = 0, T0 as listed (in memory)"})


# ------------------------------------------------------------------------------------------------ V2
def V2(p):
    purpose = "Adiabatic invariant (U = 0, no pressure drop): T - T0 = dT_ad X; bounds and monotonicity"
    rows, ok = [], True
    for T0 in (625.0, 600.0):
        q = p.with_(U=0.0, T0=T0)
        s = reactor.integrate(q, mode="none")
        d = q.dT_ad
        inv = np.concatenate([s.T - T0 - d * s.X, s.y_steps[1] - T0 - d * s.y_steps[0]])
        xt = numerics.x_bound_tol(s.rtol)
        i50 = np.nonzero(s.X >= 0.5)[0]
        z50 = float(np.interp(0.5, s.X[i50[0] - 1:i50[0] + 1], s.z[i50[0] - 1:i50[0] + 1])) if i50.size and i50[0] > 0 else float("nan")
        row = {"case": f"T0 = {T0:g} K, U = 0", "status": s.status.value, "max_invariant_err_K": float(np.abs(inv).max()),
               "X_min": float(s.X.min()), "X_max": float(s.X.max()),
               "X_monotonic": bool(np.all(np.diff(s.X) >= -xt)),
               "T_monotonic": bool(np.all(np.diff(s.T) >= -TOL["invariant_rel"] * d)),
               "X_L": float(s.X[-1]), "T_L": float(s.T[-1]), "T0_plus_dTad_XL": float(T0 + d * s.X[-1]),
               "z50_light_off_descriptor_m": z50}
        ok = (ok and s.ok and row["max_invariant_err_K"] <= TOL["invariant_rel"] * d and row["X_min"] >= -xt
              and row["X_max"] <= 1.0 + xt and row["X_monotonic"] and row["T_monotonic"])
        rows.append(row)
    tol = TOL["invariant_rel"] * p.dT_ad
    return record("V2", purpose, f"max |T - T0 - dT_ad X| <= {tol:.2e} K (invariant_rel x dT_ad); 0<=X<=1; X,T monotonic",
                  "; ".join(f"{r['case']}: max err {r['max_invariant_err_K']:.2e} K, X(L) {r['X_L']:.9f}, "
                            f"T(L) {r['T_L']:.3f} K" for r in rows),
                  f"{max(r['max_invariant_err_K'] for r in rows):.2e} K", f"{tol:.2e} K", ok,
                  "Adiabatic outlet temperatures are MATHEMATICAL excursions (T0 + dT_ad X), far beyond any documented "
                  "operating temperature; they verify the invariant only. z50 is the Phase 2 light-off descriptor "
                  "(descriptive, not a criterion).", {"rows": rows, "dT_ad_K": p.dT_ad})


# ------------------------------------------------------------------------------------------------ V3
def V3(p):
    purpose = ("Pressure-drop MATHEMATICAL validation of dp/dzeta = -Pi (T/T_c)/p (no catalyst geometry used): "
               "isothermal analytic profiles, Pi relation, outlet and inlet specifications, heated identity")
    q = p.with_(dH_rxn=0.0, U=0.0)          # T = T0 = T_c exactly -> T/T_c = 1
    rows, ok = [], True
    for Pi in (0.20, 0.45, p.Pi_max):
        sr = shooting.solve_outlet(q, Pi)
        s = sr.solution
        if not good(s, sr):
            rows.append({"case": f"outlet Pi = {Pi:g}", "status": sr.status.value, "detail": sr.detail})
            ok = False
            continue
        P0a = q.P_ref * math.sqrt(1.0 + 2.0 * Pi)
        err = np.abs(s.P - np.sqrt(P0a ** 2 - 2.0 * Pi * q.P_ref ** 2 * s.z / q.L))
        Pi_rec = ((s.P[0] / q.P_ref) ** 2 - (s.P[-1] / q.P_ref) ** 2) / 2.0
        row = {"case": f"outlet P(L) = P_ref, Pi = {Pi:g}", "status": s.status.value, "P0_numerical_Pa": sr.P0,
               "P0_analytic_Pa": P0a, "P0_rel_err": abs(sr.P0 - P0a) / P0a, "max_abs_err_Pa": float(err.max()),
               "max_rel_err": float(err.max() / q.P_ref), "outlet_residual": sr.residual,
               "Pi_recovered": float(Pi_rec), "Pi_abs_err": float(abs(Pi_rec - Pi)),
               "strictly_decreasing": bool(np.all(np.diff(s.P) < 0.0)), "P_min_over_P_ref": float(s.P.min() / q.P_ref)}
        ok = ok and row["max_rel_err"] <= TOL["P_rel"] and row["Pi_abs_err"] <= TOL["solution_rel"] and row["strictly_decreasing"]
        rows.append(row)
    for Pi in (0.10, p.Pi_inlet_max):
        s = reactor.integrate(q, mode="inlet", Pi=Pi)
        err = np.abs(s.P - q.P_ref * np.sqrt(1.0 - 2.0 * Pi * s.z / q.L))
        Pi_rec = (1.0 - (s.P[-1] / q.P_ref) ** 2) / 2.0
        row = {"case": f"inlet P(0) = P_ref, Pi = {Pi:g}", "status": s.status.value,
               "PL_numerical_Pa": float(s.P[-1]), "PL_analytic_Pa": q.P_ref * math.sqrt(1.0 - 2.0 * Pi),
               "max_abs_err_Pa": float(err.max()), "max_rel_err": float(err.max() / q.P_ref),
               "Pi_recovered": float(Pi_rec), "Pi_abs_err": float(abs(Pi_rec - Pi)),
               "strictly_decreasing": bool(np.all(np.diff(s.P) < 0.0))}
        ok = (ok and s.ok and row["max_rel_err"] <= TOL["P_rel"] and row["Pi_abs_err"] <= TOL["solution_rel"]
              and row["strictly_decreasing"])
        rows.append(row)
    Pi = 0.45
    sr = shooting.solve_outlet(p, Pi)
    s = sr.solution
    I_T = reactor.step_quadrature(s, lambda z, y: y[1] / p.Tc) / p.L
    Pi_id = ((s.P[0] / p.P_ref) ** 2 - (s.P[-1] / p.P_ref) ** 2) / (2.0 * I_T)
    heated = {"case": "heated (base parameters), outlet Pi = 0.45: P(L)^2 = P0^2 - 2 Pi P_ref^2 int(T/T_c)",
              "status": sr.status.value, "int_T_over_Tc": I_T, "Pi_from_identity": Pi_id,
              "rel_err": abs(Pi_id - Pi) / Pi}
    ok = ok and good(s, sr) and heated["rel_err"] <= TOL["quadrature_rel"]
    rows.append(heated)
    m = max(r.get("max_rel_err", 0.0) for r in rows)
    return record("V3", purpose, f"max |P_num - P_ana|/P_ref <= {TOL['P_rel']:g}; Pi recovered to "
                  f"{TOL['solution_rel']:g}; heated identity to {TOL['quadrature_rel']:g}",
                  f"max profile rel err {m:.2e} over 5 isothermal cases; heated identity rel err {heated['rel_err']:.2e}",
                  f"{m:.2e}", f"{TOL['P_rel']:g}", ok,
                  "Mathematical/implementation validation only: Pi is a study parameter; no particle diameter, voidage "
                  "or bulk density is used, so the absolute industrial pressure drop is NOT validated here.",
                  {"rows": rows, "overrides": "isothermal cases: dH_rxn = 0, U = 0 (in memory)"})


# ------------------------------------------------------------------------------------------------ V4
def V4(p):
    purpose = "Global energy balance: Qrho c_p [T(L) - T0] = (-dH) N_x0 X(L) - Q_cool (Phase 2 energy equation)"
    rows, ok = [], True
    for name, mode, Pi in CASES(p):
        s, sr = solve(p, mode, Pi)
        if not good(s, sr):
            rows.append({"case": name, "status": "FAILED"})
            ok = False
            continue
        heat = (-p.dH_rxn) * p.N_x0 * s.X[-1]
        sens = energy.heat_capacity_flow(p) * (s.T[-1] - p.T0)
        Qa = float(s.Q[-1])
        Qq = reactor.step_quadrature(s, lambda z, y: energy.q_cool(y[1], p))
        Hq = reactor.step_quadrature(s, lambda z, y: energy.q_rx(kinetics.rate(y[0], y[1], y[2], p), p))
        row = {"case": name, "heat_released_W": heat, "sensible_heat_W": sens, "Q_cool_aux_W": Qa,
               "Q_cool_quadrature_W": Qq, "residual_aux_W": sens - (heat - Qa),
               "rel_residual_aux": abs(sens - (heat - Qa)) / heat, "residual_quad_W": sens - (heat - Qq),
               "rel_residual_quad": abs(sens - (heat - Qq)) / heat, "heat_release_quadrature_W": Hq,
               "rel_err_heat_release": abs(Hq - heat) / heat}
        ok = (ok and row["rel_residual_aux"] <= TOL["invariant_rel"] and row["rel_residual_quad"] <= TOL["quadrature_rel"]
              and row["rel_err_heat_release"] <= TOL["quadrature_rel"])
        rows.append(row)
    return record("V4", purpose, f"rel. residual <= {TOL['invariant_rel']:g} (auxiliary Q) and <= "
                  f"{TOL['quadrature_rel']:g} (independent quadrature)",
                  "; ".join(f"{r['case']}: Q_cool {r.get('Q_cool_aux_W', float('nan')):.4f} W, residual "
                            f"{r.get('residual_aux_W', float('nan')):.2e} W (aux) / {r.get('residual_quad_W', float('nan')):.2e} W (quad)"
                            for r in rows),
                  f"max rel {max(max(r.get('rel_residual_aux', 1), r.get('rel_residual_quad', 1)) for r in rows):.2e}",
                  f"{TOL['invariant_rel']:g} / {TOL['quadrature_rel']:g}", ok,
                  "Sign convention: q_cool = 2 pi R U (T - T_c) > 0 is heat removed from the gas to the salt; "
                  "heat release = (-dH) N_x0 X(L) > 0.", {"rows": rows})


# ------------------------------------------------------------------------------------------------ V5
def V5(p):
    purpose = "Solver-tolerance independence (BDF): production tolerance vs 100x tighter; loose run shows convergence"
    r2, a2 = numerics.tolerances(p)
    settings = (("loose", 1e-6, 1e-8), ("production", r2, a2), ("tight", TOL["tight_rtol"], TOL["tight_rtol"] * 1e-2))
    rows, ok = [], True
    for name, mode, Pi in CASES(p):
        m = {}
        for lab, rt, at in settings:
            s, sr = solve(p, mode, Pi, rtol=rt, atol_scaled=at)
            m[lab] = hotspot.metrics(s, p) if good(s, sr) else None
        if m["production"] is None or m["tight"] is None:
            rows.append({"case": name, "status": "FAILED"})
            ok = False
            continue
        d, g = compare(m["tight"], m["production"], p)
        dl = compare(m["tight"], m["loose"], p)[0] if m["loose"] else {}
        ok = ok and g
        rows.append({"case": name, "phase2_vs_tight": d, "loose_vs_tight_abs": {k: v["abs_diff"] for k, v in dl.items()},
                     "metrics_tight": m["tight"]})
    w = [worst(r["phase2_vs_tight"]) for r in rows if "phase2_vs_tight" in r]
    return record("V5", purpose, "every metric |production - tight| within its tolerance",
                  "; ".join(f"{r['case']}: worst {k} diff {d:.2e} (tol {t:.2e})" for r, (k, d, t) in zip(rows, w)),
                  f"max diff/tol = {max(d / t for _, d, t in w):.2e}" if w else "n/a", "diff/tol <= 1", ok,
                  "Metrics: X(L), T(L), T_max, z_hs, P0, P(L), Q_cool. Loose (1e-6) vs tight differences are larger, "
                  "showing convergence toward the tight result.", {"rows": rows, "settings": settings})


# ------------------------------------------------------------------------------------------------ V6

def V6(p):
    purpose = ("Independent stiff-solver comparison: BDF (primary) vs Radau (REQUIRED cross-check); LSODA run and "
               "recorded as an additional solver")
    rows, radau_ok, bdf_ok, lsoda_limits, lsoda_disagree = [], True, True, [], []
    for name, mode, Pi in CASES(p):
        m, runs = {}, {}
        for solver in ("BDF", "Radau", "LSODA"):
            s, sr = solve(p, mode, Pi, solver=solver)
            runs[solver] = (s, sr)
            m[solver] = hotspot.metrics(s, p) if good(s, sr) else None
        if m["BDF"] is None:
            bdf_ok = False
            rows.append({"case": name, "solver": "BDF", "status": "NO VALID SOLUTION"})
            continue
        for solver in ("Radau", "LSODA"):
            if m[solver] is None:
                s, sr = runs[solver]
                trials = []
                if sr is not None:
                    for P0, st in zip(sr.scan_P0, sr.scan_status):
                        if st != PASS:
                            t_ = reactor.integrate(p, mode="outlet", Pi=Pi, P0=float(P0), solver=solver)
                            tb = reactor.integrate(p, mode="outlet", Pi=Pi, P0=float(P0), solver="BDF")
                            trials.append({"P0_over_Pref": float(P0 / p.P_ref), "status": t_.status.value,
                                           "detail": t_.detail, "BDF_same_trial_status": tb.status.value,
                                           "BDF_same_trial_Xmax_minus_1": float(tb.X.max() - 1.0)})
                row = {"case": name, "solver": solver, "status": "NO VALID SOLUTION",
                       "detail": sr.detail if sr is not None else s.detail, "failed_trials": trials}
                rows.append(row)
                if solver == "Radau":
                    radau_ok = False
                else:
                    lsoda_limits.append(row)
                continue
            d, g = compare(m["BDF"], m[solver], p)
            rows.append({"case": name, "solver": solver, "vs_BDF": d, "pass": g})
            if solver == "Radau":
                radau_ok = radau_ok and g
            elif not g:
                lsoda_disagree.append(name)
    guard_only = all(r["failed_trials"] and all(t_["status"] == "INVALID_PHYSICAL_STATE" for t_ in r["failed_trials"])
                     for r in lsoda_limits)
    ok = bdf_ok and radau_ok and not lsoda_disagree and guard_only
    disp = "FAIL" if not ok else (PWL if lsoda_limits else "PASS")
    w = [worst(r["vs_BDF"]) for r in rows if "vs_BDF" in r and r["solver"] == "Radau"]
    lim = "; ".join(f"{r['case']}: " + ", ".join(f"P0 = {t_['P0_over_Pref']:.4f} P_ref -> {t_['status']} ({t_['detail']}); "
                                                   f"BDF on the same trial: {t_['BDF_same_trial_status']}, "
                                                   f"max X - 1 = {t_['BDF_same_trial_Xmax_minus_1']:.1e}"
                                                   for t_ in r["failed_trials"]) for r in lsoda_limits)
    return record("V6", purpose, "BDF and Radau agree within every metric tolerance in all cases (Radau = required "
                  "independent solver); LSODA results recorded",
                  "Radau: " + "; ".join(f"{r['case']} worst {k} {d:.2e} (tol {t:.2e})" for r, (k, d, t) in
                                        zip([r for r in rows if r.get("solver") == "Radau" and "vs_BDF" in r], w)) +
                  f" | LSODA: valid and in agreement in {sum(1 for r in rows if r.get('solver') == 'LSODA' and r.get('pass'))}"
                  f"/{len(CASES(p))} cases",
                  f"Radau max diff/tol = {max(d / t for _, d, t in w):.2e}" if w else "n/a", "diff/tol <= 1", ok,
                  "BDF + Radau is the primary independent-solver verification pair; the Phase 3 requirement of 'at least "
                  "one additional suitable solver' is satisfied by Radau. LSODA is retained as a DOCUMENTED UNSUITABLE "
                  "cross-check for the deep-runaway outlet-pressure regime: " + (lim or "none") + ". This is a "
                  "solver-regime limitation (LSODA's interpolated conversion overshoots the equilibrium X = 1 by its "
                  "error allowance and the guard correctly rejects the state); it is NOT evidence that the physical "
                  "model is invalid.",
                  {"rows": rows, "radau_pass": radau_ok, "lsoda_limitations": lsoda_limits,
                   "lsoda_disagreements_on_valid_solutions": lsoda_disagree}, disposition=disp)


# ------------------------------------------------------------------------------------------------ V7
def V7(p):
    purpose = "Dimensional model vs independently coded dimensionless model (own groups, equations and shooting)"
    cases = [("base T0 = T_c = 625 K, Pi = 0", p, "none", 0.0),
             ("base, outlet P(L) = P_ref, Pi = Pi_max", p, "outlet", p.Pi_max),
             ("base, inlet P(0) = P_ref, Pi = Pi_inlet_max", p, "inlet", p.Pi_inlet_max),
             ("physical T0 = 640 K, T_c = 625 K, Pi = 0", p.with_(T0=640.0), "none", 0.0)]
    rows, ok = [], True
    for name, q, mode, Pi in cases:
        s, sr = solve(q, mode, Pi)
        d = dimensionless.shoot_outlet(q, Pi)["result"] if mode == "outlet" else dimensionless.integrate(q, mode, Pi)
        if not good(s, sr) or d is None or d["status"] != PASS:
            rows.append({"case": name, "status": "FAILED"})
            ok = False
            continue
        m_nd = {"X_L": float(d["X"][-1]), "T_L": float(d["T"][-1]), "T_max": float(d["T_hs"]), "z_hs": float(d["z_hs"]),
                "P0": float(d["P"][0]), "P_L": float(d["P"][-1]), "Q_cool": float(d["Q_cool"])}
        rm, g = compare(hotspot.metrics(s, q), m_nd, q)
        Pd = reactor.evaluate(s, d["z"])[2]
        dP, dp = float(np.max(np.abs(Pd - d["P"]))), float(np.max(np.abs(Pd / q.P_ref - d["p"])))
        prof = dP <= TOL["P_rel"] * q.P_ref and dp <= TOL["P_rel"]
        ok = ok and g and prof
        rows.append({"case": name, "metrics": rm, "P_profile_max_abs_diff_Pa": dP, "p_profile_max_abs_diff": dp,
                     "groups": d["groups"]})
    w = [worst(r["metrics"]) for r in rows if "metrics" in r]
    return record("V7", purpose, "every metric and the (normalized) pressure profile within tolerance",
                  "; ".join(f"{r['case']}: worst {k} {dd:.2e} (tol {t:.2e}), max|dp| {r['p_profile_max_abs_diff']:.1e}"
                            for r, (k, dd, t) in zip(rows, w)),
                  f"max diff/tol = {max(dd / t for _, dd, t in w):.2e}" if w else "n/a", "diff/tol <= 1", ok,
                  "Groups at base: " + ", ".join(f"{k} = {v:.5g}" for k, v in rows[0].get("groups", {}).items()),
                  {"rows": rows})


# ------------------------------------------------------------------------------------------------ V8
def V8(p):
    purpose = "Guard system: valid and intentionally invalid cases (P > 0.05 P_ref, 0 <= X <= 1, dP/dz < 0, T > 0, failures)"
    cases = []

    def add(name, expected, actual, detail=""):
        cases.append({"case": name, "expected": expected.value, "actual": actual.value, "pass": expected == actual,
                      "detail": detail})

    s = reactor.integrate(p)
    add("valid baseline", PASS, s.status, s.detail)
    sr = shooting.solve_outlet(p, p.Pi_max)
    s2 = sr.solution
    dP = pressure.dPdz(s2.T, s2.P, p, p.Pi_max)
    add("valid outlet Pi = Pi_max", PASS, s2.status,
        f"min P/P_ref {s2.P.min() / p.P_ref:.4f} (> 0.05); max dP/dz {dP.max():.4g} Pa/m (< 0); X in [{s2.X.min():.2e}, {s2.X.max():.9f}]")
    add("inlet spec with Pi = 0.40 > Pi_inlet_max (Phase 2 feasible range)", RunStatus.OUT_OF_MODEL_RANGE,
        reactor.integrate(p, mode="inlet", Pi=0.40).status)
    q = p.with_(dH_rxn=0.0, U=0.0)
    s3 = reactor.integrate(q, mode="inlet", Pi=0.60, config="verification", check_ranges=False)
    z_ana = q.L * (1.0 - q.P_min_frac ** 2) / (2.0 * 0.60)
    zg = s3.guard_event["z"] if s3.guard_event else float("nan")
    loc_ok = abs(zg - z_ana) <= TOL["guard_z_rel_L"] * q.L
    add("FORCED inlet Pi = 0.60 isothermal (guard test): P reaches 0.05 P_ref", RunStatus.INVALID_PHYSICAL_STATE,
        s3.status, f"guard at z = {zg:.9f} m, analytic {z_ana:.9f} m, |err| {abs(zg - z_ana):.1e} m "
        f"({'OK' if loc_ok else 'LOCATION ERROR'}); classified as invalid state, not runaway")
    add("Pi = -0.10 (range check)", RunStatus.OUT_OF_MODEL_RANGE, reactor.integrate(p, mode="inlet", Pi=-0.10).status)
    s4 = reactor.integrate(q, mode="inlet", Pi=-0.10, config="verification", check_ranges=False)
    add("FORCED Pi = -0.10 (dP/dz > 0 must be detected)", RunStatus.INVALID_PHYSICAL_STATE, s4.status, s4.detail)
    add("initial X0 = 1.2", RunStatus.INVALID_PHYSICAL_STATE, reactor.integrate(p, X0=1.2).status)
    add("initial X0 = -0.1", RunStatus.INVALID_PHYSICAL_STATE, reactor.integrate(p, X0=-0.1).status)
    xt = numerics.x_bound_tol(numerics.tolerances(p)[0])
    zz, Tok, Pok = np.zeros(5), np.full(5, p.T0), np.full(5, p.P_ref)
    # excursions defined RELATIVE to the allowance (10 x rtol): just inside -> PASS, clearly outside -> INVALID
    for xv, exp in ((1.0 + 10.0 * xt, RunStatus.INVALID_PHYSICAL_STATE), (1.0 + 0.1 * xt, PASS),
                    (-10.0 * xt, RunStatus.INVALID_PHYSICAL_STATE)):
        v = guards.check_solution(np.array([0.0, 0.2, 0.5, 0.9, xv]), Tok, Pok, zz, p, False, xt)
        add(f"synthetic profile X_end = {xv!r} (allowance {xt:.0e})", exp, v[0][0] if v else PASS)
    v = guards.check_solution(np.linspace(0, 0.5, 5), np.array([p.T0, p.T0, 0.0, p.T0, p.T0]), Pok, zz, p, False, xt)
    add("synthetic profile with T = 0", RunStatus.INVALID_PHYSICAL_STATE, v[0][0] if v else PASS)
    v = guards.check_solution(np.linspace(0, 0.5, 5), Tok, np.linspace(p.P_ref, 1.1 * p.P_ref, 5), np.full(5, 1.0), p, True, xt)
    add("synthetic pressure-drop profile with rising P", RunStatus.INVALID_PHYSICAL_STATE, v[0][0] if v else PASS)
    add("T0 = -10 K", RunStatus.INVALID_INPUT, reactor.integrate(p.with_(T0=-10.0)).status)
    add("k_m = NaN", RunStatus.INVALID_INPUT, reactor.integrate(p.with_(k_m=float("nan"))).status)
    add("T0 = 700 K (outside Phase 2 range 600-680 K)", RunStatus.OUT_OF_MODEL_RANGE, reactor.integrate(p.with_(T0=700.0)).status)
    add("normalized configuration with T_c != T0", RunStatus.OUT_OF_MODEL_RANGE,
        reactor.integrate(p.with_(T0=640.0), config="normalized").status)
    s5 = reactor.integrate(p, fault_after_z=0.7)
    add("NaN injected into the RHS for z >= 0.7 m", RunStatus.SOLVER_FAILURE, s5.status, s5.detail)
    try:
        reactor.evaluate(s, [1.1 * p.L])
        st = PASS
    except GuardError as e:
        st = e.status
    add("dense-output interpolation outside [0, L]", RunStatus.SOLVER_FAILURE, st)
    sr6 = shooting.solve_outlet(p, p.Pi_max, bracket=(p.P_ref, 1.01 * p.P_ref), n_scan=3)
    add("shooting with an invalid bracket", RunStatus.ROOT_FAILURE, sr6.status, sr6.detail)
    r7 = sensitivity.find_Ucrit(p, U_bracket=(600.0, 750.0), n_scan=2)
    add("U_crit bracket without a sign change", RunStatus.ROOT_FAILURE, r7["status"], r7["detail"])
    no_runaway_status = all("RUNAWAY" not in s.value for s in RunStatus)
    ok = all(c["pass"] for c in cases) and loc_ok and no_runaway_status
    n_ok = sum(c["pass"] for c in cases)
    return record("V8", purpose, "every case returns its expected status; guard location = analytic; no status means runaway",
                  f"{n_ok}/{len(cases)} cases as expected; guard location error {abs(zg - z_ana):.1e} m",
                  f"{len(cases) - n_ok} mismatches", "0 mismatches", ok,
                  "Statuses: " + ", ".join(s.value for s in RunStatus) + ". Pressure-guard violations are "
                  "INVALID_PHYSICAL_STATE and numerical failures SOLVER_FAILURE/ROOT_FAILURE, never runaway.",
                  {"cases": cases})


# ------------------------------------------------------------------------------------------------ V9

def V9(p):
    purpose = ("Benchmark reproduction: published Rawlings & Ekerdt companion code (Figs 6.35/6.36); diagnostic "
               "comparison of published-code behaviour, the documented '631 K' statement and the Phase 2 boundary")
    rt = TOL["fd_solver_rtol"]
    Nxf = RE.yxfeed * RE.Nf
    q_src = p.with_(R_gas=RE.Rg)          # the source's own gas constant, for like-for-like comparison
    rows, okA, claim_ok = [], True, True
    for Tf in RE.Tfeed:                   # (a)-(c) the four published feed temperatures
        _, Nx_ref, T_ref, sref = RE.run(Tf, rtol=rt, atol=[rt * Nxf, rt * RE.Tm])
        Y = reactor.evaluate(reactor.integrate(q_src.with_(T0=Tf)), RE.z)
        eN = float(np.max(np.abs(q_src.N_x0 * (1.0 - Y[0]) * 1e-3 - Nx_ref)) / Nxf)
        eT = float(np.max(np.abs(Y[1] - T_ref) / T_ref))
        Y2 = reactor.evaluate(reactor.integrate(p.with_(T0=Tf)), RE.z)
        eN2 = float(np.max(np.abs(p.N_x0 * (1.0 - Y2[0]) * 1e-3 - Nx_ref)) / Nxf)
        eT2 = float(np.max(np.abs(Y2[1] - T_ref) / T_ref))
        _, Nx_pub, T_pub, _ = RE.run(Tf)
        rows.append({"T_feed_K": Tf, "a_rel_err_Nx": eN, "a_rel_err_T": eT, "b_CODATA_R_rel_err_Nx": eN2,
                     "b_CODATA_R_rel_err_T": eT2, "b_CODATA_R_max_abs_dT_K": float(np.max(np.abs(Y2[1] - T_ref))),
                     "c_published_settings_max_abs_dT_K": float(np.max(np.abs(T_pub - T_ref)))})
        okA = okA and sref.status == 0 and eN <= TOL["V9_code_rel"] and eT <= TOL["V9_code_rel"]
        claim_ok = claim_ok and eN2 <= TOL["V9_code_rel"] and eT2 <= TOL["V9_code_rel"]
    # (e) diagnostic at the requested feed temperatures
    diag = []
    for Tf in (615.0, 620.0, 625.0, 628.0, 630.0, 631.0):
        _, Npub, Tpub, spub = RE.run(Tf)                          # published script, published settings
        s_ex = reactor.integrate(q_src.with_(T0=Tf))              # identical equations/constants (V9a), exact hot spot
        hs = hotspot.detect(s_ex, q_src)
        c_src = sensitivity.criterion(q_src.with_(T0=Tf))
        c_prod = sensitivity.criterion(p.with_(T0=Tf))
        diag.append({"T_feed_K": Tf, "published_script_status": int(spub.status),
                     "published_T_max_on_200pt_grid_K": float(Tpub.max()), "published_X_L": float(1.0 - Npub[-1] / Nxf),
                     "exact_T_max_K": hs["T_hs"], "exact_z_hs_m": hs["z_hs"], "hs_kind": hs["kind"],
                     "criterion_S_int_sourceR": c_src["S_int"], "criterion_g_sourceR": c_src["g"],
                     "criterion_S_int_CODATA_R": c_prod["S_int"], "criterion_g_CODATA_R": c_prod["g"],
                     "criterion_side": "amplified (runaway side)" if c_prod["g"] > 0 else "attenuated (safe side)"})
    for a, b in zip(diag[:-1], diag[1:]):
        b["dTmax_dT0_exact_from_previous_K_per_K"] = (b["exact_T_max_K"] - a["exact_T_max_K"]) / (b["T_feed_K"] - a["T_feed_K"])
    T_doc = float(p.source_audit["re_documented_runaway_feed_T_K"])
    rT = sensitivity.find_T0crit(p, U=p.U, config="physical", n_scan=9, xtol=TOL["T0_xtol"])
    rTs = sensitivity.find_T0crit(q_src, U=q_src.U, config="physical", n_scan=9, xtol=TOL["T0_xtol"])
    T0c, T0cs = rT.get("T0_crit", float("nan")), rTs.get("T0_crit", float("nan"))
    d = {r["T_feed_K"]: r for r in diag}
    exc = d[630.0]["published_T_max_on_200pt_grid_K"] - d[625.0]["published_T_max_on_200pt_grid_K"]
    steep = max(diag[1:], key=lambda r: r["dTmax_dT0_exact_from_previous_K_per_K"])
    classification = {
        "visible thermal excursion/onset": (
            f"NOT SUPPORTED by the accessible code: the published script's own hot spot rises from "
            f"{d[625.0]['published_T_max_on_200pt_grid_K']:.1f} K (625 K) to {d[628.0]['published_T_max_on_200pt_grid_K']:.1f} K "
            f"(628 K) and is {d[630.0]['published_T_max_on_200pt_grid_K']:.1f} K with X(L) = {d[630.0]['published_X_L']:.5f} at "
            f"630 K; the steepest listed interval is {steep['T_feed_K'] - 3 if steep['T_feed_K'] == 628.0 else steep['T_feed_K']:g}"
            f"-{steep['T_feed_K']:g} K ({steep['dTmax_dT0_exact_from_previous_K_per_K']:.1f} K per K). The steep rise is "
            f"therefore already present below {T_doc:g} K in the accessible code."),
        "plotting-based threshold": (
            "NOT SUPPORTED by the accessible code: the published script plots only T_feed = 615, 620, 625, 630 K, and "
            f"its 630 K curve is already a full excursion (+{exc:.0f} K over the 625 K curve). The printed figure itself "
            "was not accessible, so the figure cannot be checked."),
        "rounded engineering statement": "CANNOT BE DETERMINED from accessible evidence (no stated rounding or tolerance).",
        "different definition of runaway": ("CANNOT BE DETERMINED: the accessible record gives only the phrase 'runaway "
                                            "above about 631 K' with no definition; it is not the Phase 2 sensitivity-based "
                                            f"boundary ({T0c:.3f} K)."),
        "conclusion": (f"The documented {T_doc:g} K statement is INCONSISTENT with the accessible published code (whose steep "
                       f"hot-spot rise lies between 625 and 628 K and whose own 630 K case gives T_max "
                       f"{d[630.0]['published_T_max_on_200pt_grid_K']:.1f} K on its grid with X(L) = {d[630.0]['published_X_L']:.5f}) "
                       "and is NOT reproduced by the adopted sensitivity-based boundary. Its "
                       "definition cannot be established from the accessible sources; it is recorded as a benchmark-"
                       "documentation contradiction, and no parameter, equation or criterion was changed in response."),
    }
    disp = PWL if okA else "FAIL"
    return record("V9", purpose, f"(a) code-level reproduction <= {TOL['V9_code_rel']:g} relative (pass criterion); "
                  f"(e) documented '{T_doc:g} K' statement compared and classified (no pass criterion - documented contradiction)",
                  f"(a) max rel err N_x {max(r['a_rel_err_Nx'] for r in rows):.1e}, T {max(r['a_rel_err_T'] for r in rows):.1e}; "
                  f"(e) sensitivity-based boundary {T0c:.3f} K (CODATA R) / {T0cs:.3f} K (source R) vs documented {T_doc:g} K",
                  f"(a) {max(max(r['a_rel_err_Nx'], r['a_rel_err_T']) for r in rows):.1e}; (e) {T_doc - T0c:.3f} K (not a failure criterion)",
                  f"{TOL['V9_code_rel']:g}", okA,
                  "REPRODUCED: published code equations + constants at all four published feed temperatures. "
                  f"NOT REPRODUCED: the documented '{T_doc:g} K' statement (contradicted by the published code itself). "
                  "NOT ACCESSIBLE: published figure data. NOT COMPARABLE: ZB2020 critical values (different dataset). "
                  "Pre-registered (e) acceptance (|T0_crit - 631 K| <= 2.5 K) was NOT met; by the Phase 3 decision it is "
                  "recorded as a benchmark-documentation contradiction, not an implementation failure. Morbidelli-Varma "
                  "generalized criterion not evaluated because the primary source was not accessed.",
                  {"rows": rows, "diagnostic": diag, "classification_of_631K": classification,
                   "boundary_CODATA_R": rT, "boundary_source_R": rTs,
                   "pre_registered_e_check": {"criterion": "|T0_crit - 631 K| <= 2.5 K", "met": bool(abs(T0c - T_doc) <= TOL["V9_onset_halfwidth_K"])},
                   "documentation_findings": {
                       "R_note": (f"Phase 2 note 'R difference within the V9 tolerance' is {'CONFIRMED' if claim_ok else 'NOT CONFIRMED'}: "
                                  f"CODATA R vs 8.314 changes the published-case profiles by up to "
                                  f"{max(r['b_CODATA_R_max_abs_dT_K'] for r in rows):.2f} K pointwise"),
                       "published_settings_accuracy_K": max(r["c_published_settings_max_abs_dT_K"] for r in rows)},
                   "comparison_configuration": "R_gas = 8.314 (source value) for (a) and the source-R columns"},
                  disposition=disp)


# ------------------------------------------------------------------------------------------------ V10
def _fd_rows(fwd, fd, labels, steps, h2):
    """Per quantity: the largest step h* whose FD agrees with FD(h*/10) to fd_stability_rel is the numerically
    stable FD estimate; forward vs FD(h*) must agree to fd_rel. The Phase 2-step comparison is also recorded."""
    rows, ok = [], True
    for j, lab in enumerate(labels):
        errs = {f"{h:g}": abs(fwd[j] - fd[h][j]) / max(abs(fd[h][j]), 1e-300) for h in steps}
        stab = {}
        for h in steps:
            h10 = [k for k in steps if math.isclose(k, h / 10.0, rel_tol=1e-9)]
            if h10:
                stab[f"{h:g}"] = abs(fd[h][j] - fd[h10[0]][j]) / max(abs(fd[h10[0]][j]), 1e-300)
        stable = [h for h in steps if stab.get(f"{h:g}", float("inf")) <= TOL["fd_stability_rel"]]
        h_star = max(stable) if stable else None
        e_star = errs[f"{h_star:g}"] if h_star is not None else float("inf")
        good_row = e_star <= TOL["fd_rel"]
        ok = ok and good_row
        rows.append({"quantity": lab, "forward": float(fwd[j]), "fd": {f"{h:g}": float(fd[h][j]) for h in steps},
                     "rel_err": errs, "fd_stability_h_vs_h/10": stab, "stable_step": h_star,
                     "rel_err_at_stable_step": e_star, "rel_err_at_phase2_step": errs.get(f"{h2:g}"),
                     "phase2_step_is_stable": stab.get(f"{h2:g}", float("inf")) <= TOL["fd_stability_rel"],
                     "pass": good_row})
    return rows, ok


def V10(p):
    purpose = "Forward sensitivity S = d(X,T,P)/dT0 (T_c fixed) vs central finite differences at several steps"
    steps, h2 = list(TOL["fd_steps_rel"]), TOL["fd_step_rel_phase2"]
    fdk = dict(solver="Radau", rtol=TOL["fd_solver_rtol"], atol_scaled=TOL["fd_solver_rtol"] * 1e-2)
    out, ok = {}, True
    # Case A: baseline
    c = sensitivity.criterion(p)
    zp = np.array([0.25 * p.L, 0.5 * p.L, 0.75 * p.L, p.L, c["z_Sint"]])
    Sf = c["res"]["dense"](zp)[4:7]
    fwd = list(Sf[1]) + [Sf[0][3]]
    labels = [f"S_T(z = {z:.4f} m)" for z in zp] + ["S_X(L)"]
    fd = {}
    for hr in steps:
        h = hr * p.T0
        yp = reactor.evaluate(reactor.integrate(p.with_(T0=p.T0 + h), **fdk), zp)
        ym = reactor.evaluate(reactor.integrate(p.with_(T0=p.T0 - h), **fdk), zp)
        d = (yp - ym) / (2.0 * h)
        fd[hr] = list(d[1]) + [d[0][3]]
    rowsA, okA = _fd_rows(fwd, fd, labels, steps, h2)
    out["A_baseline"] = {"rows": rowsA, "z_Sint": c["z_Sint"], "S_int": c["S_int"]}
    # Case B: outlet P(L) = P_ref (total derivative incl. dP0/dT0)
    Pi = 0.20
    cB = sensitivity.criterion(p, mode="outlet", Pi=Pi)
    zpB = np.array([0.25 * p.L, 0.5 * p.L, 0.75 * p.L, p.L, cB["z_Sint"]])
    SB = cB["res"]["dense"](zpB)[4:7]
    fwdB = list(SB[1]) + [cB["dP0dT0"]]
    labelsB = [f"S_T(z = {z:.4f} m)" for z in zpB] + ["dP0/dT0 [Pa/K]"]
    stepsB = [h for h in steps if h <= 1e-4]
    fdB = {}
    for hr in stepsB:
        h = hr * p.T0
        res = []
        for sgn in (1.0, -1.0):
            brk = (cB["P0"] * (1 - 2e-3), cB["P0"] * (1 + 2e-3))
            sr = shooting.solve_outlet(p.with_(T0=p.T0 + sgn * h), Pi, n_scan=2, bracket=brk,
                                       xtol_rel=1e-14, res_tol=1e-10, **fdk)
            if sr.status != PASS:
                raise RuntimeError(f"FD shooting failed: {sr.detail}")
            res.append((reactor.evaluate(sr.solution, zpB), sr.P0))
        d = (res[0][0] - res[1][0]) / (2.0 * h)
        fdB[hr] = list(d[1]) + [(res[0][1] - res[1][1]) / (2.0 * h)]
    rowsB, okB = _fd_rows(fwdB, fdB, labelsB, stepsB, h2)
    out["B_outlet_Pi_0.20"] = {"rows": rowsB, "z_Sint": cB["z_Sint"], "S_int": cB["S_int"], "P0": cB["P0"]}
    ok = okA and okB
    allrows = rowsA + rowsB
    e = max(r["rel_err_at_stable_step"] for r in allrows)
    e2 = max(r["rel_err_at_phase2_step"] for r in allrows)
    n_unstable2 = sum(not r["phase2_step_is_stable"] for r in allrows)
    return record("V10", purpose, f"per quantity: rel err <= {TOL['fd_rel']:g} at the largest numerically stable FD step "
                  f"(FD(h) vs FD(h/10) <= {TOL['fd_stability_rel']:g})",
                  f"{len(allrows)} quantities; max rel err at stable step {e:.1e}; at the Phase 2 step {e2:.1e} "
                  f"({n_unstable2} quantities for which the Phase 2 step is NOT a stable FD step)",
                  f"{e:.1e}", f"{TOL['fd_rel']:g}", ok,
                  "Parameter: T0 (T_c held fixed). Case B holds P(L) = P_ref, so the FD re-solves the shooting and the "
                  "forward result includes dP0/dT0. FD reference solves: Radau, rtol 1e-12.", out)


# ------------------------------------------------------------------------------------------------ V11
def V11(p):
    purpose = "U_crit root finding (Phase 2 criterion, physical configuration, T0 = T_c = 625 K, Pi = 0)"
    xt = TOL["Ucrit_xtol"]
    r = sensitivity.find_Ucrit(p, n_scan=16, xtol=xt)
    if r["status"] != PASS:
        return record("V11", purpose, "unique bracketed root", r["status"].value, "n/a", xt, False, r["detail"], r)
    U = r["U_crit"]
    r_full = sensitivity.find_Ucrit(p, n_scan=2, xtol=xt)
    r_shift = sensitivity.find_Ucrit(p, U_bracket=(U - 37.0, U + 53.0), n_scan=4, xtol=xt)
    # offset bracket kept in the root neighbourhood: BDF at rtol 1e-12 can fail (step-size underflow) on the
    # augmented system in DEEP runaway (e.g. U = 362.9, g ~ 700), which is irrelevant to the root's tolerance shift
    r_tight = sensitivity.find_Ucrit(p, U_bracket=(U - 4.7, U + 5.3), n_scan=3, xtol=xt,
                                     rtol=TOL["tight_rtol"], atol_scaled=TOL["tight_rtol"] * 1e-2)
    r_rad = sensitivity.find_Ucrit(p, U_bracket=(U - 20.0, U + 20.0), n_scan=3, xtol=xt, solver="Radau")
    r_rep = sensitivity.find_Ucrit(p, n_scan=16, xtol=xt)
    alt = {"full bracket [0, 750]": r_full, "shifted bracket": r_shift, f"rtol {TOL['tight_rtol']:g}": r_tight,
           "Radau": r_rad, "identical rerun": r_rep}
    diffs = {k: (abs(v["U_crit"] - U) if v["status"] == PASS else float("inf")) for k, v in alt.items()}
    ok = (r["converged"] and r["sign_changes"] == 1 and r["subscan_sign_changes"] == 1 and r["local_monotonic_at_root"]
          and all(d <= xt for d in diffs.values()) and r_rep["U_crit"] == U)
    return record("V11", purpose, f"valid bracket, one sign change (scan and 9-point sub-scan), strictly monotonic through the root, "
                  f"converged; all variants within {xt} W/(m2 K)",
                  f"U_crit = {U:.4f} W/(m2 K); g at root {r['g_at_root']:.1e}; "
                  f"{r['iterations']} Brent iterations; bracket {r['bracket']}",
                  "max |dU_crit| = " + f"{max(diffs.values()):.2e}", f"{xt} W/(m2 K)", ok,
                  f"Scan g(U): globally monotonic = {r['globally_monotonic']} (g jumps on the SAFE side where the interior "
                  f"maximum of S_T disappears; no root there); local g through root = {[round(x, 5) for x in r['local_g']]}; "
                  f"at root: {r['at_root'].get('rule')}, "
                  f"dS_T/dz(0) = {r['at_root'].get('dS0', float('nan')):.3g} 1/m. U_crit is the value of the SENSITIVITY-BASED PARAMETRIC RUNAWAY BOUNDARY, not a dynamic-stability "
                  "boundary, saddle-node or bifurcation point.",
                  {"main": r, "variants": {k: {"U_crit": v.get("U_crit"), "status": v["status"], "abs_diff": diffs[k]}
                                           for k, v in alt.items()}})


# ------------------------------------------------------------------------------------------------ V12
def V12(p):
    purpose = "Consistency of the boundary reached along U (T0 fixed) and along T0 (U fixed), both configurations"
    rU = sensitivity.find_Ucrit(p, n_scan=16, xtol=TOL["Ucrit_xtol_V12"])
    if rU["status"] != PASS:
        return record("V12", purpose, "U root", rU["status"].value, "n/a", "n/a", False, rU["detail"], rU)
    U = rU["U_crit"]
    rows, ok = [], True
    for cfg in ("physical", "normalized"):
        rT = sensitivity.find_T0crit(p, U=U, config=cfg, T0_bracket=(615.0, 635.0), n_scan=5, xtol=TOL["T0_xtol"])
        d = abs(rT.get("T0_crit", float("nan")) - p.T0)
        good_row = rT["status"] == PASS and d <= TOL["T0_boundary_resolution"]
        ok = ok and good_row
        rows.append({"configuration": cfg, "T0_crit": rT.get("T0_crit"), "abs_diff_K": d, "status": rT["status"],
                     "scan_T0": rT.get("scan_x"), "scan_g": rT.get("scan_g"), "pass": good_row})
    return record("V12", purpose, f"T0_crit at U = U_crit(625 K) returns 625 K within {TOL['T0_boundary_resolution']} K",
                  f"U_crit(625 K) = {U:.4f}; " + "; ".join(f"{r['configuration']}: T0_crit = {r['T0_crit']:.4f} K" for r in rows),
                  "; ".join(f"{r['abs_diff_K']:.1e} K" for r in rows), f"{TOL['T0_boundary_resolution']} K", ok,
                  "Both configurations coincide at T0 = T_c = 625 K; they differ along the T0 direction (T_c fixed vs "
                  "T_c = T0), so each returning 625 K checks its own implementation. Single boundary point only (no map).",
                  {"U_crit": U, "rows": rows})


# ------------------------------------------------------------------------------------------------ V13
def V13(p):
    purpose = "Outlet-pressure shooting: P(L) vs guessed P0 scan, monotonicity, bracket, root, residual, tighter tolerance"
    rows, ok = [], True
    for Pi in TOL["V13_Pi_values"]:
        lo, hi = shooting.analytic_bracket(p, Pi)   # includes the rtol margin
        P0s = np.linspace(lo, hi, TOL["V13_scan_points"])
        tab = []
        for P0 in P0s:
            F, s = shooting.residual(P0, p, Pi)
            tab.append({"P0_Pa": float(P0), "status": s.status.value,
                        "P_L_Pa": float(s.P[-1]) if s.ok else None,
                        "guard_z_m": s.guard_event["z"] if s.guard_event else None, "F": float(F)})
        F = np.array([t["F"] for t in tab])
        mono = bool(np.all(np.diff(F) > 0.0))
        sgn = np.sign(F)
        nsc = int(np.sum(sgn[:-1] * sgn[1:] < 0))
        sr = shooting.solve_outlet(p, Pi)
        srt = shooting.solve_outlet(p, Pi, rtol=TOL["tight_rtol"], atol_scaled=TOL["tight_rtol"] * 1e-2)
        s = sr.solution
        prof_ok = good(s, sr) and bool(np.all(np.diff(s.P) < 0.0)) and s.P.min() > p.P_min
        dP0 = abs(sr.P0 - srt.P0) / p.P_ref if srt.status == PASS else float("inf")
        row = {"Pi": Pi, "bracket_Pa": (lo, hi), "bracket_valid": bool(F[0] <= 0.0 <= F[-1]), "monotonic": mono,
               "sign_changes": nsc, "P0_Pa": sr.P0, "P0_over_Pref": sr.P0 / p.P_ref,
               "isothermal_P0_over_Pref": math.sqrt(1.0 + 2.0 * Pi), "outlet_residual": sr.residual,
               "brent_iterations": sr.iterations, "profile_ok": prof_ok, "P0_tight_Pa": srt.P0, "P0_rel_change_tight": dP0,
               "T_max_K": hotspot.detect(s, p)["T_hs"] if good(s, sr) else None, "scan": tab}
        ok = (ok and row["bracket_valid"] and mono and nsc == 1 and sr.status == PASS
              and abs(sr.residual) <= TOL["shoot_residual"] and prof_ok and dP0 <= TOL["P_rel"])
        rows.append(row)
    return record("V13", purpose, f"monotonic P(L)(P0), one sign change, |P(L)/P_ref - 1| <= {TOL['shoot_residual']:g}, "
                  f"P0 stable to {TOL['P_rel']:g} under rtol 1e-10",
                  "; ".join(f"Pi = {r['Pi']}: P0 = {r['P0_Pa'] / 1e3:.4f} kPa ({r['P0_over_Pref']:.4f} P_ref), residual "
                            f"{r['outlet_residual']:.1e}, monotonic {r['monotonic']}" for r in rows),
                  "max residual " + f"{max(abs(r['outlet_residual']) for r in rows):.1e}", f"{TOL['shoot_residual']:g}", ok,
                  "Scan: 41 guessed inlet pressures over the rigorous bracket; trials stopped by the pressure guard are "
                  "invalid states (continued residual used for bracketing only).", {"rows": rows})
