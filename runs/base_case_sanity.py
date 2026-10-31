"""Phase 3 BASE-CASE SANITY RUN (section 8) and dimensional vs dimensionless check (section 9).

NOT a scientific result. One base case (T0 = T_c = 625 K, U = 373 W/(m2 K), y_x0 = 0.019, physical
configuration) in the two Phase 2 pressure configurations: (A) benchmark baseline, no pressure drop;
(B) primary pressure-drop configuration, P(L) = P_ref, at the parameter-file value Pi = Pi_max.
"""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
for _p in (str(ROOT), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402
import yaml  # noqa: E402

from pbr import (composition, dimensionless, energy, hotspot, numerics, params, plotting,  # noqa: E402
                 pressure, reactor, shooting, storage)

TOL = yaml.safe_load((ROOT / "verification" / "tolerances.yaml").read_text())


def run(p, mode, Pi, solver=None):
    if mode == "outlet":
        sr = shooting.solve_outlet(p, Pi, solver=solver)
        return sr.solution, sr
    return reactor.integrate(p, mode=mode, solver=solver), None


def one(p, label, mode, Pi):
    s, sr = run(p, mode, Pi)
    m, hs = hotspot.metrics(s, p), hotspot.detect(s, p)
    xt = numerics.x_bound_tol(s.rtol)
    pd_on = mode in pressure.PD_MODES
    dPdz = pressure.dPdz(s.T, s.P, p, Pi) if pd_on else np.zeros_like(s.z)
    heat = (-p.dH_rxn) * p.N_x0 * m["X_L"]
    resid = energy.heat_capacity_flow(p) * (m["T_L"] - p.T0) - (heat - m["Q_cool"])
    T_lo, T_hi = min(p.T0, p.Tc), max(p.T0, p.Tc) + p.dT_ad
    checks = {
        "0 <= X <= 1": bool(s.X.min() >= -xt and s.X.max() <= 1 + xt),
        "X non-decreasing": bool(np.all(np.diff(s.X) >= -xt)),
        "min(T0,Tc) <= T <= max(T0,Tc) + dT_ad (rigorous bounds)": bool(s.T.min() >= T_lo - TOL["T_abs_K"] and s.T.max() <= T_hi + TOL["T_abs_K"]),
        "P > 0.05 P_ref everywhere": bool(s.P.min() > p.P_min),
        "dP/dz < 0 everywhere" if pd_on else "P = P_ref everywhere (baseline)":
            bool(np.all(dPdz < 0)) if pd_on else bool(np.all(s.P == p.P_ref)),
        "P(L) = P_ref": bool(abs(m["P_L"] / p.P_ref - 1) <= TOL["shoot_residual"]),
        "global energy balance closes": bool(abs(resid) / heat <= TOL["invariant_rel"]),
    }
    s2, sr2 = run(p, mode, Pi, solver="Radau")
    m2 = hotspot.metrics(s2, p)
    cross = {k: abs(m2[k] - m[k]) for k in ("X_L", "T_L", "T_max", "z_hs", "P0", "P_L", "Q_cool")}
    # section 9: independent dimensionless formulation
    d = dimensionless.shoot_outlet(p, Pi)["result"] if mode == "outlet" else dimensionless.integrate(p, mode, Pi)
    Pd = reactor.evaluate(s, d["z"])[2]
    items = {"outlet conversion X(L) [-]": (m["X_L"], float(d["X"][-1]), TOL["X_abs"]),
             "outlet temperature T(L) [K]": (m["T_L"], float(d["T"][-1]), TOL["T_abs_K"]),
             "maximum temperature T_max [K]": (m["T_max"], float(d["T_hs"]), TOL["T_abs_K"]),
             "hot-spot location z_hs [m]": (m["z_hs"], float(d["z_hs"]), TOL["z_rel_L"] * p.L),
             "pressure profile max|P_dim - P_nd| [Pa]": (None, float(np.max(np.abs(Pd - d["P"]))), TOL["P_rel"] * p.P_ref),
             "normalized pressure profile max|p_dim - p_nd| [-]": (None, float(np.max(np.abs(Pd / p.P_ref - d["p"]))), TOL["P_rel"])}
    comp = {}
    for k, (a, b, tol) in items.items():
        diff = abs(b - a) if a is not None else b
        comp[k] = {"dimensional": a, "dimensionless": b if a is not None else None, "abs_diff": diff,
                   "rel_diff": (diff / abs(a)) if a else diff / (p.P_ref if "[Pa]" in k else 1.0),
                   "tolerance": tol, "pass": diff <= tol}
    worst = max(comp.values(), key=lambda r: r["abs_diff"] / r["tolerance"])
    return s, {
        "label": label, "mode": mode, "Pi": Pi, "solver": s.solver, "rtol": s.rtol, "atol": s.atol,
        "solver_status": s.status.value, "solver_detail": s.detail, "nfev": s.nfev,
        "shooting": None if sr is None else {"status": sr.status.value, "P0": sr.P0, "residual": sr.residual,
                                              "bracket": sr.bracket, "monotonic_scan": sr.monotonic,
                                              "sign_changes": sr.sign_changes, "brent_iterations": sr.iterations},
        "guard_event": s.guard_event, "guard_status": "no guard triggered" if s.guard_event is None else s.guard_event,
        "inlet_pressure_Pa": m["P0"], "outlet_pressure_Pa": m["P_L"], "pressure_drop_Pa": m["dP"],
        "outlet_conversion": m["X_L"], "outlet_temperature_K": m["T_L"], "maximum_temperature_K": m["T_max"],
        "hot_spot_location_m": m["z_hs"], "hot_spot_kind": hs["kind"], "cooling_duty_W": m["Q_cool"],
        "energy_residual_W": resid, "O2_that_stoichiometry_would_consume_at_X_L": composition.o2_depletion_fraction(m["X_L"], p),
        "physical_checks": checks, "radau_crosscheck_abs_diff": cross,
        "dimensionless_comparison": comp,
        "dimensionless_max_discrepancy_over_tolerance": worst["abs_diff"] / worst["tolerance"],
        "dimensionless_pass": all(r["pass"] for r in comp.values()),
    }


def main():
    p = params.load()
    out, sols = {}, []
    for key, label, mode, Pi in (("A_baseline", "baseline, no pressure drop (benchmark)", "none", 0.0),
                                 ("B_outlet_Pi_max", f"P(L) = P_ref, Pi = Pi_max = {p.Pi_max:g}", "outlet", p.Pi_max)):
        s, r = one(p, label, mode, Pi)
        out[key] = r
        sols.append((s, label))
    out["B_outlet_Pi_max"]["interpretation"] = (
        "MATHEMATICAL/DIAGNOSTIC result, not a validated industrial operating prediction. With the outlet "
        "specification, Pi = 0.70 requires P0 = {:.4f} P_ref, so the difference from run A combines a PRESSURE-LEVEL "
        "effect (higher inlet pressure; rate ~ P^2) and a PRESSURE-DECLINE effect (pressure falling along the bed). "
        "It is NOT a pure pressure-drop effect; the two are to be separated in Phase 4 (constant P = P0 comparison). "
        "T_max = {:.1f} K is a runaway-branch mathematical temperature excursion.".format(
            out["B_outlet_Pi_max"]["inlet_pressure_Pa"] / p.P_ref, out["B_outlet_Pi_max"]["maximum_temperature_K"]))
    out["A_baseline"]["interpretation"] = "Diagnostic benchmark-baseline result (no pressure drop); not a scientific result."
    out["note"] = ("Sanity run only - not a scientific result. Temperatures above the documented operating "
                   "temperatures (salt 618-658 K; observed industrial hot spot ~713 K, context only) are mathematical "
                   "model outputs: the kinetic validity range is undocumented and no engineering limit is sourced.")
    res = ROOT / "results" / "sanity"
    storage.save(res / "base_case_sanity.json", out, p)
    plotting.profiles(sols, p, res / "base_case_sanity.png", "Phase 3 SANITY RUN - base case T0 = Tc = 625 K, U = 373 W/(m2 K)",
                      "Diagnostic only, not a scientific result; kinetic validity range undocumented (runaway-branch temperatures are mathematical outputs)")
    for key in ("A_baseline", "B_outlet_Pi_max"):
        r = out[key]
        print("\n   INTERPRETATION: " + r["interpretation"])
        print(f"\n== {r['label']} ==  solver {r['solver']} rtol {r['rtol']:g}  status {r['solver_status']}  guard: {r['guard_status']}")
        if r["shooting"]:
            sh = r["shooting"]
            print(f"   shooting {sh['status']}: P0 = {sh['P0']:.4f} Pa, residual {sh['residual']:.2e}, monotonic {sh['monotonic_scan']}, sign changes {sh['sign_changes']}")
        print(f"   P_in {r['inlet_pressure_Pa']/1e3:.4f} kPa | P_out {r['outlet_pressure_Pa']/1e3:.4f} kPa | dP {r['pressure_drop_Pa']/1e3:.4f} kPa")
        print(f"   X(L) {r['outlet_conversion']:.6f} | T(L) {r['outlet_temperature_K']:.4f} K | T_max {r['maximum_temperature_K']:.4f} K at z = {r['hot_spot_location_m']:.5f} m ({r['hot_spot_kind']})")
        print(f"   cooling duty {r['cooling_duty_W']:.4f} W | energy residual {r['energy_residual_W']:.1e} W | O2 fraction stoichiometry would consume {r['O2_that_stoichiometry_would_consume_at_X_L']:.3f}")
        print("   physical checks: " + ", ".join(f"{k}: {'OK' if v else 'FAIL'}" for k, v in r["physical_checks"].items()))
        print("   Radau cross-check |diff|: " + ", ".join(f"{k} {v:.1e}" for k, v in r["radau_crosscheck_abs_diff"].items()))
        print(f"   dimensional vs dimensionless: {'PASS' if r['dimensionless_pass'] else 'FAIL'}, max discrepancy/tolerance {r['dimensionless_max_discrepancy_over_tolerance']:.2e}")
        for k, c in r["dimensionless_comparison"].items():
            print(f"      {k:52s} abs {c['abs_diff']:.2e}  rel {c['rel_diff']:.2e}  tol {c['tolerance']:.2e}  {'PASS' if c['pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
