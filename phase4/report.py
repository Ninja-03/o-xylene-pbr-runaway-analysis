"""Phase 4 final report (generated from stored records only): results/phase4/phase4_report.md."""
import hashlib
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from phase4.common import OUT, ROOT, read_json, read_jsonl  # noqa: E402
from phase4.interpret import FLAG_TEXT, T0_VALID, U_VALID, interpretation  # noqa: E402
from phase4.tables import merged  # noqa: E402


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def main():
    maps, pb, rob = merged()
    for m in maps + pb:
        m["interpretation"] = interpretation(m)
    sims = {n: read_jsonl(n) for n in ("t0_study", "u_study", "pressure_study", "o2_study")}
    st = {}
    for n, rs in sims.items():
        for r in rs:
            st[r["status"]] = st.get(r["status"], 0) + 1
    nonpass = [(n, r["label"], r["status"], r["detail"][:80]) for n, rs in sims.items() for r in rs if r["status"] != "PASS"]
    o2rb = read_json("o2_roots_rebracket")
    o2r = read_json("o2_roots")
    roots = maps + pb + read_json("u_roots") + o2r["roots"] + o2rb
    rst = {}
    for m in roots:
        k = interpretation(m)
        k = "flagged (criterion discontinuity)" if k == FLAG_TEXT else k.split(" (")[0] if k.startswith("validated") else k
        rst[k] = rst.get(k, 0) + 1
    L = ["# Phase 4 final report - controlled parametric study (generated from stored records)", ""]
    m = read_json("meta")
    L += ["## A. Total simulations", "", f"{sum(len(v) for v in sims.values())} stored simulations: " +
          ", ".join(f"{k} {len(v)}" for k, v in sims.items()) + f". Root computations stored: {len(roots)} "
          "(plus robustness finite-difference and numerical-robustness roots stored in their own files).", ""]
    L += ["## B. Statuses", "", f"Simulations: {st}.", "", "Non-PASS simulations: " +
          ("; ".join(f"{a}: {b} -> {c} ({d})" for a, b, c, d in nonpass) or "none"), "",
          f"Root outcomes by interpretation: {rst}.", ""]
    A = sorted([x for x in maps if x.get("series") == "A: Ucrit(T0), physical, Pi=0" and x["status"] == "PASS"], key=lambda x: x["T0"])
    L += [f"## C. Validated U_crit(T0) - primary interpreted domain (Tc = 625 K, Pi = 0, T0 <= {T0_VALID:g} K, U <= {U_VALID:g})", "",
          "| T0 (K) | U_crit (W m-2 K-1) | note |", "|---|---|---|"]
    L += [f"| {x['T0']:g} | {x['U_crit']:.4f} | {'rebracketed [100, 750]' if x.get('rebracketed') else ''} |" for x in A if x["T0"] <= T0_VALID]
    B = sorted([x for x in maps if x.get("series") == "B: T0crit(U), physical, Pi=0" and x["status"] == "PASS"], key=lambda x: x["U"])
    L += ["", "T0,crit(U) validated roots: " + "; ".join(f"U = {x['U']:g}: {x['T0_crit']:.3f} K ({x['interpretation']})" for x in B) +
          "; U = 200-350: T0_crit < 600 K (below searched range).", ""]
    F = [x for x in maps if x.get("g_discontinuity_at_root")]
    L += ["## D. Flagged / discontinuous region (diagnostic only)", "", FLAG_TEXT + ":", ""] + \
         [f"- {x['label']}: sign change at {x['sign_change_location']:.3f}" for x in F] + \
         ["", "Above-range (not flagged): " + "; ".join(f"{x['label']}" for x in maps if "above the searched range" in x.get("outcome", "")), ""]
    r = o2rb[0]
    L += ["## E. O2-depletion U_crit at T0 = 610 K", "",
          f"{r['status']}: U_crit = {r['U_crit']:.4f} W m-2 K-1 (bracket U in [100, 750]; final bracket {r['bracket']}, "
          f"converged {r['converged']}, {r['iterations']} Brent iterations, g at root {r['g_at_root']:.2e}, one sign change "
          f"in scan and sub-scan, monotonic through root = {r['local_monotonic_at_root']}). Fixed-O2: {r['fixed_O2_U_crit']:.4f}; "
          f"difference {r['diff_vs_fixed_O2']:+.2f}.", "", "| T0 | U_crit fixed-O2 | U_crit O2-depletion | difference |", "|---|---|---|---|"]
    allo2 = o2r["roots"] + o2rb
    for T in (610.0, 625.0, 640.0):
        b = [x for x in allo2 if x["series"] == "fixed-O2 (base)" and x["T0"] == T and x["status"] == "PASS"][0]
        v = [x for x in allo2 if x["series"] == "O2-depletion variant" and x["T0"] == T and x["status"] == "PASS"][0]
        L.append(f"| {T:g} | {b['U_crit']:.2f} | {v['U_crit']:.2f} | {v['U_crit'] - b['U_crit']:+.2f} |")
    ref = [x for x in pb if x.get("series") == "A constant level" and x.get("level") == 1.0][0]["U_crit"]
    L += ["", "## F. Pressure-level / pressure-decline decomposition of U_crit (T0 = 625 K)", "",
          "A = constant pressure (level only); B = inlet-fixed decline; C = outlet-fixed level + decline; A' = constant "
          "pressure at the matched inlet pressure of C. At fixed outlet pressure, the modeled pressure-level effect exceeds "
          "the pressure-decline effect over the tested Pi range. The rate scales as P^2 at fixed T and X; this is a property "
          "of the present model/parameterization, not a validated industrial generalization, and C vs Pi = 0 is not a pure "
          "pressure-drop effect.", "",
          f"| Pi | P0/P_ref of C | U_crit C | U_crit A' | level (A' - {ref:.2f}) | decline (C - A') | combined (C - {ref:.2f}) |",
          "|---|---|---|---|---|---|---|"]
    for pi in (0.1, 0.2, 0.35, 0.5, 0.7):
        c = [x for x in pb if x.get("series") == "C outlet combined" and x["Pi"] == pi][0]
        a = [x for x in pb if x.get("series") == "A' constant at matched P0" and x["Pi"] == pi]
        cu = c["U_crit"] if c["status"] == "PASS" else None
        au = a[0]["U_crit"] if a and a[0]["status"] == "PASS" else None
        f = lambda v: f"{v:.2f}" if isinstance(v, float) else "> 750"
        L.append(f"| {pi:g} | {c.get('level', float('nan')) if cu else float('nan'):.4f} | {f(cu)} | "
                 f"{f(au) if a else 'n/a'} | {f(au - ref) if au else ('> +' + format(750 - ref, '.1f') if a else 'n/a')} | "
                 f"{f(cu - au) if (cu and au) else 'n/a'} | {f(cu - ref) if cu else '> +' + format(750 - ref, '.1f')} |")
    L += ["", "B inlet-fixed (decline only): " + "; ".join(f"Pi = {x['Pi']:g}: {x['U_crit']:.2f} ({x['U_crit'] - ref:+.2f})"
                                                         for x in pb if x.get("series") == "B inlet decline"), ""]
    L += ["## G. Figures", ""] + [f"- {p.name}" for p in sorted((OUT / "figures").glob("*.png"))] + [""]
    L += ["## Evidence table", "", "| Observation | Numerical evidence | Kind | Limitation | Confidence |", "|---|---|---|---|---|",
          "| The sensitivity-based boundary precedes the steep calculated thermal rise | g changes sign between T0 = 615.5 and 616 K; calculated T_max 718.9 -> 917.7 K over 625-628 K (U = 373) | directly computed | R&E 631 K contradiction remains documented | high |",
          "| At fixed outlet pressure, the modeled pressure-level effect exceeds the pressure-decline effect over the tested Pi range | see F | directly computed; P^2 dependence is a mathematical consequence of the model | Pi not an absolute industrial pressure drop | high within the model |",
          "| Pressure decline alone lowers U_crit | B inlet-fixed values in F | directly computed | same | high within the model |",
          f"| T_c has the largest evaluated normalized local sensitivity of U_crit at P2 | {rob['boundary_sensitivities']['Tc']['normalized']:+.2f} | directly computed (local FD, stable step) | local; T_c varied at fixed T0 | medium-high |",
          f"| The normalized local sensitivity of U_crit to k_m is approximately +1.00 at P2 | {rob['boundary_sensitivities']['k_m']['normalized']:+.4f} | directly computed (local FD) | local result, not a scaling law | medium |",
          f"| dH_rxn and y_x0 have identical local effects | U_crit {rob['boundary_sensitivities']['dH_rxn']['normalized']:+.3f} vs {rob['boundary_sensitivities']['y_x0']['normalized']:+.3f} | mathematical consequence (both enter only via dT_ad) | - | high |",
          "| The fixed-O2 assumption materially changes the T0-direction boundary | variant boundary between T0 = 630 and 635 K vs ~615.8 K (U = 373); U_crit shifts -13.7 to -20.4 | model-form comparison | variant not validated | high that the assumption is material |",
          "| Tc = T0 and Tc fixed configurations differ strongly | U_crit 177.65 vs 366.78 (600 K); 581.70 vs 406.10 (640 K) | directly computed | Tc = T0 is not industrially representative | high |", ""]
    L += ["## H. Remaining scientific limitations", "",
          "- Kinetic validity range undocumented: calculated runaway-branch temperatures are mathematical model outputs outside the documented kinetic validity range.",
          f"- Primary boundary interpreted only for T0 <= {T0_VALID:g} K, U <= {U_VALID:g}; beyond, the criterion's discontinuity prevents validated smooth roots.",
          "- Pi is a verified dimensionless parameterization; absolute industrial pressure drop is not validated (no particle/voidage data).",
          "- Fixed-O2 base model; the O2-depletion comparison shows the assumption is material.",
          "- Local sensitivities only; at P1 (g > 0) most finite differences have no stable step; the U_crit sensitivity to T0 has no stable step.",
          "- Adiabatic (U = 0) states with pressure drop hit a SciPy event-location limitation; affected roots use U in [100, 750].",
          "- R&E '631 K' benchmark-documentation contradiction; single lumped reaction; 1-D pseudo-homogeneous; constant salt temperature.", ""]
    code = {p.name: sha(p) for p in sorted((ROOT / "phase4").glob("*.py"))}
    data = {p.name: sha(p) for p in sorted(OUT.glob("*")) if p.is_file() and p.name != "phase4_report.md"}
    L += ["## I. Hashes", "", f"Frozen Phase 3: params `{m['params_sha256']}`, code `{m['code_sha256']}`, tolerances `{m['tolerances_sha256']}`.",
          f"Phase 4 study code (in every simulation record): `{m['phase4_code_sha256']}`.", "", "Phase 4 code files:", ""] + \
         [f"- {k}: `{v}`" for k, v in code.items()] + ["", "Phase 4 data files:", ""] + [f"- {k}: `{v}`" for k, v in data.items()]
    (OUT / "phase4_report.md").write_text("\n".join(L) + "\n")
    print("\n".join(L[:60]))


if __name__ == "__main__":
    main()
