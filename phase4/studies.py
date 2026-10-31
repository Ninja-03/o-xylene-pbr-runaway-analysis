"""Phase 4 study drivers. Usage: python phase4/studies.py [study ...] (see STUDIES at the end), or 'all' 

Study design (fixed here, before results are interpreted):
  base operating point: T0 = T_c = 625 K, U = 373 W/(m2 K), y_x0 = 0.019, physical configuration (T_c fixed).
  T0 study: 600-680 K step 5 K (Phase 2 range), refined to 0.5 K in every coarse interval where T_max changes by
            > 15 K, X(L) by > 0.05, or g changes sign (refinement rule fixed a priori).
  U study:  0-750 W/(m2 K) step 25 (Phase 2 range) + base U and a 2.5-step refinement of the g sign-change interval.
  Pressure: configurations A constant P = P0 (level only), B inlet P(0) = P_ref (decline only, Pi <= 0.30),
            C outlet P(L) = P_ref (level + decline, P0 solved). Decomposition at matched thermal conditions:
            level = A(P0 = P0_C) - baseline, decline = C - A, combined = C - baseline, decline@P_ref = B - baseline.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np  # noqa: E402

from phase4 import common  # noqa: E402
from phase4.common import root_summary, run_point, write_json, write_jsonl  # noqa: E402
from phase4 import o2_variant  # noqa: E402
from pbr import params, sensitivity, shooting  # noqa: E402

P = params.load()
REP_T0 = {605.0, 615.0, 620.0, 625.0, 627.0, 628.0, 630.0, 640.0}
REP_U = {250.0, 350.0, 373.0, 380.0, 400.0, 500.0}
PIS = [0.1, 0.2, 0.35, 0.5, 0.7]
# Pressure-drop U_crit roots use U in [100, 750] instead of [0, 750]: at U = 0 (adiabatic) with the pressure equation
# active, SciPy's event location fails on the extremely steep front ("f(a) and f(b) must have different signs" ->
# SOLVER_FAILURE, correctly rejected by the validated root finder). g >> 0 there (g ~ 4e3-6e3 at U = 125), so the
# root is unaffected. Constant-pressure (no pressure equation) roots keep the full Phase 2 range [0, 750].
PD_BRACKET = (100.0, 750.0)


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------------------------------------ T0 study
def t0():
    recs = []

    def sweep(tag, cfg, mode, Pi, grid, refine, profiles):
        pts = {}
        for T0 in grid:
            q = P.with_(T0=T0, Tc=T0) if cfg == "normalized" else P.with_(T0=T0)
            pts[T0] = run_point(q, "T0", f"{tag} T0={T0:g}", mode=mode, Pi=Pi, config=cfg,
                                keep_profile=profiles and T0 in REP_T0)
        if refine:
            ks = sorted(pts)
            add = []
            for a, b in zip(ks[:-1], ks[1:]):
                ra, rb = pts[a], pts[b]
                if (abs(rb["T_max"] - ra["T_max"]) > 15 or abs(rb["X_L"] - ra["X_L"]) > 0.05
                        or np.sign(ra["g"]) != np.sign(rb["g"])):
                    add += list(np.arange(a + 0.5, b, 0.5))
            for T0 in add:
                pts[T0] = run_point(P.with_(T0=T0), "T0", f"{tag} T0={T0:g}", mode=mode, Pi=Pi, config=cfg,
                                    keep_profile=profiles and T0 in REP_T0)
        for T0 in sorted(pts):
            pts[T0]["series"] = tag
            recs.append(pts[T0])
        log(f"T0 sweep {tag}: {len(pts)} points")

    sweep("physical Pi=0", "physical", "none", 0.0, np.arange(600.0, 680.1, 5.0), True, True)
    sweep("normalized Pi=0", "normalized", "none", 0.0, np.arange(600.0, 680.1, 5.0), False, False)
    sweep("physical outlet Pi=0.2", "physical", "outlet", 0.2, np.arange(600.0, 680.1, 10.0), False, False)
    write_jsonl("t0_study", recs)


# ------------------------------------------------------------------------------------------------ U study
def u():
    recs, roots = [], []
    grid = sorted(set(np.arange(0.0, 750.1, 25.0)) | {373.0, 380.0})
    pts = {U: run_point(P.with_(U=U), "U", f"physical Pi=0 U={U:g}", keep_profile=U in REP_U) for U in grid}
    ks = sorted(pts)
    for a, b in zip(ks[:-1], ks[1:]):
        if np.sign(pts[a]["g"]) != np.sign(pts[b]["g"]):
            for U in np.arange(a + 2.5, b, 2.5):
                pts[U] = run_point(P.with_(U=U), "U", f"physical Pi=0 U={U:g}")
    for U in sorted(pts):
        pts[U]["series"] = "physical Pi=0"
        recs.append(pts[U])
    roots.append(root_summary(sensitivity.find_Ucrit(P), "U", "U_crit T0=625 physical Pi=0", "U", T0=625.0,
                              config="physical", mode="none", Pi=0.0))
    for U in np.arange(0.0, 750.1, 75.0):
        r = run_point(P.with_(U=U), "U", f"physical outlet Pi=0.2 U={U:g}", mode="outlet", Pi=0.2)
        r["series"] = "physical outlet Pi=0.2"
        recs.append(r)
    write_jsonl("u_study", recs)
    write_json("u_roots", roots)
    log(f"U study: {len(recs)} points, U_crit = {roots[0].get('U_crit')}")


# ------------------------------------------------------------------------------------------------ boundary maps
def maps_main():
    out = []
    for T0 in np.arange(600.0, 680.1, 5.0):
        out.append(root_summary(sensitivity.find_Ucrit(P.with_(T0=T0)), "map", f"A physical Ucrit(T0={T0:g})", "U",
                                series="A: Ucrit(T0), physical, Pi=0", T0=T0, config="physical", mode="none", Pi=0.0))
    log("map A physical done")
    for U in np.arange(200.0, 750.1, 50.0):
        out.append(root_summary(sensitivity.find_T0crit(P, U=U), "map", f"B physical T0crit(U={U:g})", "T0",
                                series="B: T0crit(U), physical, Pi=0", U=U, config="physical", mode="none", Pi=0.0))
    log("map B physical done")
    for T0 in np.arange(600.0, 680.1, 10.0):
        out.append(root_summary(sensitivity.find_Ucrit(P.with_(T0=T0, Tc=T0), config="normalized"), "map",
                                f"A normalized Ucrit(T0={T0:g})", "U", series="A: Ucrit(T0), normalized Tc=T0, Pi=0",
                                T0=T0, config="normalized", mode="none", Pi=0.0))
    log("map A normalized done")
    write_json("boundary_maps", out)


def maps_outlet():
    out = []
    for T0 in (610.0, 625.0, 640.0):
        out.append(root_summary(sensitivity.find_Ucrit(P.with_(T0=T0), mode="outlet", Pi=0.2, n_scan=7,
                                                        U_bracket=PD_BRACKET), "map",
                                f"A outlet Pi=0.2 Ucrit(T0={T0:g})", "U", series="A: Ucrit(T0), physical, outlet Pi=0.2",
                                T0=T0, config="physical", mode="outlet", Pi=0.2, U_bracket=PD_BRACKET))
        log(f"map A outlet Pi=0.2 T0={T0:g}: {out[-1].get('U_crit')}")
    write_json("boundary_maps_outlet", out)


# ------------------------------------------------------------------------------------------------ pressure study
def pressure_points():
    recs, decomp = [], []
    for plabel, U in (("U=373 (base)", 373.0), ("U=420", 420.0)):
        q = P.with_(U=U)
        ref = run_point(q, "pressure", f"{plabel} baseline P=P_ref", keep_profile=True)
        ref["config_name"] = "baseline"
        recs.append(ref)
        for Pi in PIS:
            C = run_point(q, "pressure", f"{plabel} C outlet Pi={Pi:g}", mode="outlet", Pi=Pi, keep_profile=True)
            C["config_name"] = "C outlet"
            A = run_point(q, "pressure", f"{plabel} A constant P=P0_C Pi={Pi:g}", mode="constant", P0=C["P0_Pa"],
                          keep_profile=True)
            A["config_name"] = "A constant P0"
            B = None
            if Pi <= P.Pi_inlet_max:
                B = run_point(q, "pressure", f"{plabel} B inlet Pi={Pi:g}", mode="inlet", Pi=Pi, keep_profile=True)
                B["config_name"] = "B inlet"
            recs += [C, A] + ([B] if B else [])
            row = {"point": plabel, "U": U, "Pi": Pi, "P0_C_over_Pref": C["P0_Pa"] / P.P_ref}
            for k in ("X_L", "T_max", "dT_hs", "z_hs", "X_at_hs", "S_int", "g", "Q_cool_W"):
                row[k] = {"baseline": ref[k], "A": A[k], "C": C[k], "B": B[k] if B else None,
                          "level_effect(A-baseline)": A[k] - ref[k], "decline_effect(C-A)": C[k] - A[k],
                          "combined(C-baseline)": C[k] - ref[k],
                          "decline_at_Pref(B-baseline)": (B[k] - ref[k]) if B else None}
            decomp.append(row)
        log(f"pressure decomposition {plabel} done")
    write_jsonl("pressure_study", recs)
    write_json("pressure_decomposition", decomp)


def pressure_ci():
    roots = []
    for lev in (1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6):
        roots.append(root_summary(sensitivity.find_Ucrit(P, mode="constant", P0=lev * P.P_ref), "pressure_boundary",
                                  f"A constant P={lev:g} P_ref", "U", series="A constant level", level=lev, Pi=0.0))
    for Pi in (0.1, 0.2, 0.3):
        roots.append(root_summary(sensitivity.find_Ucrit(P, mode="inlet", Pi=Pi, U_bracket=PD_BRACKET),
                                  "pressure_boundary", f"B inlet Pi={Pi:g}", "U", series="B inlet decline", Pi=Pi,
                                  level=1.0, U_bracket=PD_BRACKET))
        log(f"inlet root Pi={Pi:g}: {roots[-1].get('U_crit')}")
    write_json("pressure_boundary_ci", roots)


def pressure_out(pis, name):
    roots = []
    for Pi in pis:
        rc = sensitivity.find_Ucrit(P, mode="outlet", Pi=Pi, n_scan=7, U_bracket=PD_BRACKET)
        rs = root_summary(rc, "pressure_boundary", f"C outlet Pi={Pi:g}", "U", series="C outlet combined", Pi=Pi,
                          U_bracket=PD_BRACKET)
        if rs["status"] == "PASS":
            P0c = shooting.solve_outlet(P.with_(U=rc["U_crit"]), Pi).P0
            rs["level"] = P0c / P.P_ref
            roots.append(rs)
            roots.append(root_summary(sensitivity.find_Ucrit(P, mode="constant", P0=P0c), "pressure_boundary",
                                      f"A' constant P=P0_C(Ucrit) Pi={Pi:g}", "U", series="A' constant at matched P0",
                                      Pi=Pi, level=P0c / P.P_ref))
        else:
            roots.append(rs)
        log(f"outlet root Pi={Pi:g}: {rs.get('U_crit')} ({rs['outcome']})")
    write_json(name, roots)


# ------------------------------------------------------------------------------------------------ O2 model-form check
def o2():
    recs, roots = [], []
    ver = o2_variant.verify(P)
    for T0 in np.arange(600.0, 680.1, 5.0):
        b = run_point(P.with_(T0=T0), "o2", f"base T0={T0:g}", keep_profile=T0 == 625.0)
        b["series"] = "fixed-O2 (base)"
        with o2_variant.o2_depletion():
            v = run_point(P.with_(T0=T0), "o2", f"O2-depletion T0={T0:g}", variant="O2-depletion",
                          keep_profile=T0 == 625.0)
        v["series"] = "O2-depletion variant"
        recs += [b, v]
    for T0 in (610.0, 625.0, 640.0):
        roots.append(root_summary(sensitivity.find_Ucrit(P.with_(T0=T0)), "o2", f"base Ucrit(T0={T0:g})", "U",
                                  series="fixed-O2 (base)", T0=T0))
        with o2_variant.o2_depletion():
            roots.append(root_summary(sensitivity.find_Ucrit(P.with_(T0=T0)), "o2", f"O2 Ucrit(T0={T0:g})", "U",
                                      series="O2-depletion variant", T0=T0))
    b = run_point(P, "o2", "base outlet Pi=0.2", mode="outlet", Pi=0.2)
    with o2_variant.o2_depletion():
        v = run_point(P, "o2", "O2-depletion outlet Pi=0.2", mode="outlet", Pi=0.2, variant="O2-depletion")
    recs += [b, v]
    write_jsonl("o2_study", recs)
    write_json("o2_roots", {"verification": ver, "roots": roots})
    log("O2 study done; verification pass =", ver["pass"])


# ------------------------------------------------------------------------------------------------ robustness
def _fd(f, x0, one_sided, steps=(1e-2, 1e-3, 1e-4), agree=1e-2):
    """Central (or one-sided 2nd-order) differences at several relative steps; the reported value is the largest
    step whose estimate agrees with the next smaller step to `agree` (Phase 3 lesson: no blind step choice)."""
    f0 = f(x0) if one_sided else None
    est = {}
    for h in steps:
        d = h * x0
        if one_sided:
            f1, f2 = f(x0 + d), f(x0 + 2 * d)
            est[h] = {k: (-3 * f0[k] + 4 * f1[k] - f2[k]) / (2 * d) for k in f1}
        else:
            fp, fm = f(x0 + d), f(x0 - d)
            est[h] = {k: (fp[k] - fm[k]) / (2 * d) for k in fp}
    out = {}
    for k in est[steps[0]]:
        chosen = None
        for h, h2 in zip(steps[:-1], steps[1:]):
            a, b = est[h][k], est[h2][k]
            if abs(a - b) <= agree * max(abs(b), 1e-300):
                chosen = h
                break
        out[k] = {"derivative": est[chosen][k] if chosen else est[steps[-1]][k], "stable_step": chosen,
                  "estimates": {f"{h:g}": est[h][k] for h in steps}}
    return out


def _q(r):
    return {"T_max": r["T_max"], "dT_hs": r["dT_hs"], "X_L": r["X_L"], "S_int": r["S_int"]}


def robust_points():
    res = []
    pnames = ["T0", "U", "cp", "dH_rxn", "E_over_R", "k_m", "y_x0", "Tc"]
    for plabel, U in (("P1 base T0=625 U=373 (g>0)", 373.0), ("P2 T0=625 U=420 (g<0)", 420.0)):
        base = P.with_(U=U)
        r0 = run_point(base, "robust", f"{plabel} reference")
        y0 = _q(r0)
        for name in pnames:
            x0 = getattr(base, name)
            fd = _fd(lambda x: _q(run_point(base.with_(**{name: x}), "robust", "fd")), x0, one_sided=(name == "cp"))
            row = {"point": plabel, "parameter": name, "value": x0,
                   "method": "FD one-sided 2nd order (value at range bound)" if name == "cp" else "FD central"}
            for k, v in fd.items():
                row[k] = {**v, "normalized": v["derivative"] * x0 / y0[k] if y0[k] else None}
            if name == "T0":
                row["dTmax_dT0_forward_sensitivity"] = r0["dTmax_dT0_forward"]
            res.append(row)
        rPi = run_point(base, "robust", f"{plabel} outlet Pi=0.2", mode="outlet", Pi=0.2)
        fdPi = _fd(lambda x: _q(run_point(base, "robust", "fd", mode="outlet", Pi=x)), 0.2, one_sided=False)
        row = {"point": plabel, "parameter": "Pi (outlet, at 0.2)", "value": 0.2, "method": "FD central"}
        for k, v in fdPi.items():
            row[k] = {**v, "normalized": v["derivative"] * 0.2 / _q(rPi)[k]}
        res.append(row)
        log(f"robust point sensitivities {plabel} done")
    write_json("robustness_points", res)


def robust_boundary():
    Uc = sensitivity.find_Ucrit(P, xtol=1e-6)["U_crit"]
    res = {"boundary_reference_Ucrit": Uc, "boundary_sensitivities": {}}
    for name in ["T0", "cp", "dH_rxn", "E_over_R", "k_m", "y_x0", "Tc"]:
        x0 = getattr(P, name)

        def ucrit(x, name=name):
            r = sensitivity.find_Ucrit(P.with_(**{name: x}), U_bracket=(Uc - 80.0, Uc + 80.0), n_scan=9, xtol=1e-6)
            return {"U_crit": r["U_crit"]}
        v = _fd(ucrit, x0, one_sided=(name == "cp"), steps=(1e-2, 1e-3))["U_crit"]
        res["boundary_sensitivities"][name] = {**v, "value": x0, "normalized": v["derivative"] * x0 / Uc}
        log(f"boundary sensitivity {name}: {v['derivative']:.5g} (stable step {v['stable_step']})")
    write_json("robustness_boundary", res)


# ------------------------------------------------------------------------------------------------ numerical robustness
def numerics_():
    out = []
    for T0 in (610.0, 625.0, 640.0):
        q = P.with_(T0=T0)
        ref = sensitivity.find_Ucrit(q)
        rad = sensitivity.find_Ucrit(q, solver="Radau")
        tight = sensitivity.find_Ucrit(q, rtol=1e-11, atol_scaled=1e-13)
        out.append({"case": f"physical Pi=0 T0={T0:g}", "BDF": ref.get("U_crit"), "Radau": rad.get("U_crit"),
                    "BDF_rtol_1e-11": tight.get("U_crit"), "statuses": [x["status"].value for x in (ref, rad, tight)],
                    "abs_diff_Radau": abs(rad.get("U_crit") - ref.get("U_crit")),
                    "abs_diff_tight": abs(tight.get("U_crit") - ref.get("U_crit"))})
        log(f"numerics T0={T0:g} done")
    ref = sensitivity.find_Ucrit(P, mode="outlet", Pi=0.2, n_scan=7, U_bracket=PD_BRACKET)
    rad = sensitivity.find_Ucrit(P, mode="outlet", Pi=0.2, n_scan=7, U_bracket=PD_BRACKET, solver="Radau")
    out.append({"case": "physical outlet Pi=0.2 T0=625", "BDF": ref.get("U_crit"), "Radau": rad.get("U_crit"),
                "statuses": [ref["status"].value, rad["status"].value],
                "abs_diff_Radau": abs(rad.get("U_crit") - ref.get("U_crit"))})
    write_json("numerical_robustness", out)
    log("numerical robustness done")


STUDIES = {"t0": t0, "u": u, "maps_main": maps_main, "maps_outlet": maps_outlet, "pressure_points": pressure_points,
           "pressure_ci": pressure_ci, "pressure_out1": lambda: pressure_out([0.1, 0.2, 0.35], "pressure_boundary_out1"),
           "pressure_out2": lambda: pressure_out([0.5, 0.7], "pressure_boundary_out2"), "o2": o2,
           "robust_points": robust_points, "robust_boundary": robust_boundary, "numerics": numerics_}

if __name__ == "__main__":
    common.write_json("meta", common.meta())
    for name in (sys.argv[1:] if sys.argv[1:] and sys.argv[1] != "all" else list(STUDIES)):
        log(f"=== {name} ===")
        STUDIES[name]()
    log("ALL REQUESTED STUDIES DONE")
