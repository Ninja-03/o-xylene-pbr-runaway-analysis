"""Phase 4 tables (CSV + Markdown), generated only from stored records in results/phase4."""
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from phase4.common import OUT, read_json, read_jsonl  # noqa: E402
from phase4.interpret import interpretation  # noqa: E402


def merged():
    maps = read_json("boundary_maps") + read_json("boundary_maps_outlet")
    if (OUT / "boundary_maps_rebracket.json").exists():
        maps += read_json("boundary_maps_rebracket")
    pb = read_json("pressure_boundary_ci") + read_json("pressure_boundary_out1") + read_json("pressure_boundary_out2")
    rb = read_json("robustness_boundary")
    rob = {"point_sensitivities": read_json("robustness_points"), **rb}
    return maps, pb, rob

COLS = ["label", "series", "config", "mode", "Pi", "T0", "Tc", "U", "status", "P0_Pa", "PL_Pa", "dP_Pa", "X_L", "T_L",
        "T_max", "dT_hs", "z_hs", "X_at_hs", "Q_cool_W", "S_int", "S_T_max_profile", "g", "boundary_class",
        "dTmax_dT0_forward"]


def flat(r):
    d = {k: r.get(k) for k in COLS}
    d.update(T0=r["inputs"]["T0"], Tc=r["inputs"]["Tc"], U=r["inputs"]["U"])
    return d


def fmt(v, nd=4):
    if isinstance(v, float):
        return f"{v:.{nd}g}" if abs(v) < 1e4 or abs(v) == 0 else f"{v:.6g}"
    return "" if v is None else str(v)


def md(rows, cols, nd=5):
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        out.append("| " + " | ".join(fmt(r.get(c), nd) for c in cols) + " |")
    return "\n".join(out)


def main():
    lines = ["# Phase 4 tables (generated from results/phase4 records)", ""]
    meta = read_json("meta")
    lines += [f"Frozen check {meta['frozen_check']}: params `{meta['params_sha256'][:16]}`, code `{meta['code_sha256'][:16]}`, "
              f"tolerances `{meta['tolerances_sha256'][:16]}`, Phase 4 code `{meta['phase4_code_sha256'][:16]}`.", ""]
    allrecs = []
    for name in ("t0_study", "u_study", "pressure_study", "o2_study"):
        recs = read_jsonl(name)
        allrecs += recs
        with open(OUT / f"{name}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS)
            w.writeheader()
            for r in recs:
                w.writerow(flat(r))
    status = {}
    for r in allrecs:
        status[r["status"]] = status.get(r["status"], 0) + 1
    lines += [f"## Status counts over {len(allrecs)} stored simulations", "", str(status), ""]
    t0 = [flat(r) for r in read_jsonl("t0_study")]
    c = ["T0", "X_L", "T_max", "dT_hs", "z_hs", "X_at_hs", "Q_cool_W", "S_int", "g", "dTmax_dT0_forward"]
    for s in ("physical Pi=0", "normalized Pi=0"):
        lines += [f"## T0 study - {s} (U = 373)", "", md([r for r in t0 if r["series"] == s], c), ""]
    lines += ["## T0 study - physical outlet Pi=0.2", "",
              md([r for r in t0 if r["series"] == "physical outlet Pi=0.2"], ["T0", "P0_Pa", "dP_Pa"] + c[1:]), ""]
    us = [flat(r) for r in read_jsonl("u_study")]
    cu = ["U", "X_L", "T_max", "z_hs", "Q_cool_W", "S_int", "g"]
    for s in ("physical Pi=0", "physical outlet Pi=0.2"):
        lines += [f"## U study - {s} (T0 = Tc = 625 K)", "", md([r for r in us if r["series"] == s],
                                                                ["U", "P0_Pa"] + cu[1:] if "outlet" in s else cu), ""]
    lines += ["## U_crit root (T0 = 625 K, physical, Pi = 0)", "", md(read_json("u_roots"), ["label", "status", "U_crit", "g_at_root", "bracket", "iterations"]), ""]
    maps, pb, rob = merged()
    for m in maps + pb:
        m["interpretation"] = interpretation(m)
    lines += ["## Boundary map roots (primary interpreted domain: physical Pi=0, T0 <= 650 K, U <= 450)", "",
              md(maps, ["label", "T0", "U", "status", "U_crit", "T0_crit", "sign_change_location", "interpretation"]), ""]
    dec = read_json("pressure_decomposition")
    rows = []
    for d in dec:
        for k in ("X_L", "T_max", "z_hs", "S_int", "g", "Q_cool_W"):
            v = d[k]
            rows.append({"point": d["point"], "Pi": d["Pi"], "P0_C/P_ref": d["P0_C_over_Pref"], "metric": k,
                         "baseline": v["baseline"], "A const P0": v["A"], "C outlet": v["C"], "B inlet": v["B"],
                         "level (A-base)": v["level_effect(A-baseline)"], "decline (C-A)": v["decline_effect(C-A)"],
                         "combined (C-base)": v["combined(C-baseline)"], "decline@Pref (B-base)": v["decline_at_Pref(B-baseline)"]})
    lines += ["## Pressure-level vs pressure-decline decomposition (T0 = 625 K)", "",
              md(rows, list(rows[0].keys())), ""]
    lines += ["## Boundary under pressure configurations (T0 = 625 K)", "",
              md(pb, ["series", "Pi", "level", "status", "U_crit", "interpretation"]), ""]
    o2 = [flat(r) for r in read_jsonl("o2_study")]
    base = {r["T0"]: r for r in o2 if r["series"] == "fixed-O2 (base)"}
    var = {r["T0"]: r for r in o2 if r["series"] == "O2-depletion variant"}
    rows = [{"T0": T, "X_L base": base[T]["X_L"], "X_L O2": var[T]["X_L"], "T_max base": base[T]["T_max"],
             "T_max O2": var[T]["T_max"], "dT_max": var[T]["T_max"] - base[T]["T_max"], "g base": base[T]["g"],
             "g O2": var[T]["g"], "same side": (base[T]["g"] > 0) == (var[T]["g"] > 0)} for T in sorted(base) if T in var]
    o2r = read_json("o2_roots")
    lines += ["## O2-depletion model-form check", "", "Variant verification: " + str(o2r["verification"]), "", md(rows, list(rows[0].keys())), "",
              md(o2r["roots"] + (read_json("o2_roots_rebracket") if (OUT / "o2_roots_rebracket.json").exists() else []),
                 ["label", "series", "T0", "status", "U_crit", "outcome"]), ""]
    rows = []
    for r in rob["point_sensitivities"]:
        rows.append({"point": r["point"], "parameter": r["parameter"], "value": r["value"],
                     **{f"{k} norm": r[k]["normalized"] for k in ("T_max", "dT_hs", "X_L", "S_int")},
                     **{f"{k} stable step": r[k]["stable_step"] for k in ("dT_hs", "X_L")}})
    lines += ["## Normalized parameter sensitivities (phi/y) dy/dphi", "", md(rows, list(rows[0].keys())), "",
              f"Reference U_crit = {rob['boundary_reference_Ucrit']:.6f} W/(m2 K)", "",
              md([{"parameter": k, **v} for k, v in rob["boundary_sensitivities"].items()],
                 ["parameter", "value", "derivative", "normalized", "stable_step"]), ""]
    lines += ["## Numerical robustness near the boundary", "", md(sum((read_json(f"numerical_robustness_{k}") for k in ("phys", "tight", "outlet")), []),
                                                                   ["case", "BDF", "Radau", "BDF_rtol_1e-11", "abs_diff_Radau", "abs_diff_tight", "statuses"], 8), ""]
    (OUT / "tables.md").write_text("\n".join(lines) + "\n")
    print("tables written; status counts:", status)


if __name__ == "__main__":
    main()
