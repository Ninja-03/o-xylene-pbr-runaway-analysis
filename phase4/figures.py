"""Phase 4 figures, generated only from stored records in results/phase4 (no manual transcription).
Terminology: 'sensitivity-based parametric runaway boundary' only; no stability language."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from phase4.common import OUT, read_json, read_jsonl  # noqa: E402
from phase4.interpret import FLAG_TEXT, T0_VALID, U_VALID  # noqa: E402


def merged():
    maps = read_json("boundary_maps") + read_json("boundary_maps_outlet")
    if (OUT / "boundary_maps_rebracket.json").exists():
        maps += read_json("boundary_maps_rebracket")
    pb = read_json("pressure_boundary_ci") + read_json("pressure_boundary_out1") + read_json("pressure_boundary_out2")
    rb = read_json("robustness_boundary")
    rob = {"point_sensitivities": read_json("robustness_points"), **rb}
    return maps, pb, rob

FIG = OUT / "figures"
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 150, "savefig.bbox": "tight",
                     "legend.fontsize": 7.5})
FOOT = ("Base: Tc = 625 K, U = 373 W m$^{-2}$ K$^{-1}$, y$_{x0}$ = 0.019, L = 1.5 m (R&E dataset). Mathematical model "
        "outputs; calculated runaway-branch temperatures lie outside the documented kinetic validity range.")


def save(fig, name, foot=FOOT):
    fig.text(0.5, -0.02, foot, ha="center", fontsize=6.5, style="italic")
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / name)
    plt.close(fig)


def series(recs, s):
    r = sorted([x for x in recs if x.get("series") == s and x.get("status") == "PASS"],
               key=lambda x: (x["inputs"]["T0"], x["inputs"]["U"]))
    return r


def mark_side(ax, x, y, g):
    g = np.asarray(g, float)
    ax.plot(np.asarray(x)[g <= 0], np.asarray(y)[g <= 0], "o", ms=3.2, mfc="white", mec="k", lw=0, zorder=3)
    ax.plot(np.asarray(x)[g > 0], np.asarray(y)[g > 0], "o", ms=3.2, mfc="k", mec="k", lw=0, zorder=3)


def main():
    t0 = read_jsonl("t0_study")
    us = read_jsonl("u_study")
    uroot = read_json("u_roots")[0]
    maps, pb, rob = merged()
    pr = read_jsonl("pressure_study")
    o2 = read_jsonl("o2_study")
    o2r = read_json("o2_roots")["roots"]
    Tc = 625.0

    # 1-2 profiles
    rep = sorted([r for r in t0 if r.get("series") == "physical Pi=0" and "profile" in r], key=lambda r: r["inputs"]["T0"])
    for key, name, lab in (("T", "fig01_T_profiles.png", "T (K)"), ("X", "fig02_X_profiles.png", "X (-)")):
        fig, ax = plt.subplots(figsize=(5.2, 3.6))
        for r in rep:
            ax.plot(r["profile"]["z"], r["profile"][key], lw=1.2,
                    label=f"T0 = {r['inputs']['T0']:g} K (g {'>' if r['g'] > 0 else '<'} 0)")
        if key == "T":
            ax.axhline(Tc, color="0.5", ls="--", lw=0.8, label="Tc = 625 K")
        ax.set(xlabel="axial position z (m)", ylabel=lab,
               title=("Temperature" if key == "T" else "Conversion") + " profiles vs feed temperature (P = P_ref)")
        ax.legend(ncol=1, fontsize=6.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
        save(fig, name)

    # 3 pressure profiles
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    sel = [r for r in pr if r["label"].startswith("U=373") and "profile" in r]
    for r in sel:
        pi = float(r["label"].split("Pi=")[1]) if "Pi=" in r["label"] else 0.0
        if r["config_name"] == "baseline" or (pi in (0.2, 0.7) and r["config_name"] in ("C outlet", "A constant P0")) \
                or (r["config_name"] == "B inlet" and pi == 0.2):
            ls = {"baseline": "-", "C outlet": "-", "A constant P0": ":", "B inlet": "--"}[r["config_name"]]
            name = {"baseline": "baseline P = P_ref", "C outlet": "C outlet-fixed", "B inlet": "B inlet-fixed",
                    "A constant P0": "A constant P = P0 of C"}[r["config_name"]]
            ax.plot(r["profile"]["z"], np.asarray(r["profile"]["P"]) / 101000.0, ls=ls, lw=1.3,
                    label=name + ("" if r["config_name"] == "baseline" else f", Pi = {pi:g}"))
    ax.set(xlabel="axial position z (m)", ylabel="P / P_ref (-)  (P_ref = 101 kPa)",
           title="Pressure profiles: A constant P0, B inlet-fixed, C outlet-fixed (U = 373)")
    ax.legend()
    save(fig, "fig03_P_profiles.png")

    # 4-5 Tmax and Xout vs T0
    for key, name, lab in (("T_max", "fig04_Tmax_vs_T0.png", "T_max (K)"), ("X_L", "fig05_Xout_vs_T0.png", "X(L) (-)")):
        fig, ax = plt.subplots(figsize=(5.4, 3.7))
        for s, st in (("physical Pi=0", "-"), ("normalized Pi=0", "--"), ("physical outlet Pi=0.2", ":")):
            r = series(t0, s)
            x, y, g = [q["inputs"]["T0"] for q in r], [q[key] for q in r], [q["g"] for q in r]
            ax.plot(x, y, st, lw=1.1, label=s.replace("physical", "Tc = 625 K").replace("normalized", "Tc = T0"))
            mark_side(ax, x, y, g)
        ax.set(xlabel="feed temperature T0 (K)", ylabel=lab, title=f"{lab.split(' ')[0]} vs T0 (U = 373 W m$^{{-2}}$ K$^{{-1}}$)")
        ax.plot([], [], "o", mfc="k", mec="k", label="g > 0 (inlet perturbations amplified)")
        ax.plot([], [], "o", mfc="white", mec="k", label="g < 0 (attenuated)")
        ax.legend()
        save(fig, name)

    # 6 Tmax vs U
    fig, ax = plt.subplots(figsize=(5.4, 3.7))
    for s, st in (("physical Pi=0", "-"), ("physical outlet Pi=0.2", ":")):
        r = sorted([x for x in us if x.get("series") == s and x["status"] == "PASS"], key=lambda x: x["inputs"]["U"])
        x, y, g = [q["inputs"]["U"] for q in r], [q["T_max"] for q in r], [q["g"] for q in r]
        ax.plot(x, y, st, lw=1.1, label=s.replace("physical", "T0 = Tc = 625 K,"))
        mark_side(ax, x, y, g)
    ax.axvline(uroot["U_crit"], color="C3", lw=0.9, label=f"U_crit = {uroot['U_crit']:.2f} (Pi = 0)")
    ax.set(xlabel="overall heat-transfer coefficient U (W m$^{-2}$ K$^{-1}$)", ylabel="T_max (K)",
           title="Hot-spot temperature vs cooling (T0 = Tc = 625 K)")
    ax.legend()
    save(fig, "fig06_Tmax_vs_U.png")

    # 7 Ucrit vs T0 -- primary: validated smooth roots only, T0 <= 650 K
    o2all = o2r + (read_json("o2_roots_rebracket") if (OUT / "o2_roots_rebracket.json").exists() else [])
    fig, ax = plt.subplots(figsize=(6.2, 4.3))
    ax.axvspan(T0_VALID, 682, color="0.88", zorder=0)
    ax.text(T0_VALID + 1.5, 60, "criterion discontinuity /\nboundary not interpreted", fontsize=7, color="0.3")
    for s, st, mk, lab in (("A: Ucrit(T0), physical, Pi=0", "-", "o", "Tc = 625 K, Pi = 0"),
                           ("A: Ucrit(T0), normalized Tc=T0, Pi=0", "--", "s", "Tc = T0, Pi = 0"),
                           ("A: Ucrit(T0), physical, outlet Pi=0.2", ":", "^", "Tc = 625 K, outlet-fixed Pi = 0.2")):
        r = sorted([m for m in maps if m.get("series") == s and m["status"] == "PASS" and m["T0"] <= T0_VALID],
                   key=lambda m: m["T0"])
        ln, = ax.plot([m["T0"] for m in r], [m["U_crit"] for m in r], st, marker=mk, ms=3.5, lw=1.1, label=lab)
        up = [m for m in maps if m.get("series") == s and "above the searched range" in m.get("outcome", "")]
        if up:
            ax.plot([m["T0"] for m in up], [770] * len(up), marker=mk, ls="", mfc="none", color=ln.get_color(),
                    clip_on=False, label=f"{lab}: U_crit > 750 (above searched range)")
    fl = [m for m in maps if m.get("g_discontinuity_at_root") and m.get("kind") == "U"]
    ax.plot([m["T0"] for m in fl], [m["sign_change_location"] for m in fl], "x", color="0.4", ms=6,
            label="Flagged sign changes associated with discontinuity of the\nsensitivity criterion; not interpreted as validated boundary roots")
    ov = sorted([m for m in o2all if m["series"] == "O2-depletion variant" and m["status"] == "PASS"], key=lambda m: m["T0"])
    ax.plot([m["T0"] for m in ov], [m["U_crit"] for m in ov], "P", color="C4", ms=6, ls="",
            label="O2-depletion variant (model-form check), Tc = 625 K, Pi = 0")
    ax.set(xlim=(598, 682), ylim=(0, 790), xlabel="feed temperature T0 (K)", ylabel="U_crit (W m$^{-2}$ K$^{-1}$)",
           title="U_crit(T0): sensitivity-based parametric runaway boundary")
    ax.legend(fontsize=6.3, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    save(fig, "fig07_Ucrit_vs_T0.png")

    # 8 boundary map -- primary interpreted domain T0 <= 650 K, U <= 450
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    A = sorted([m for m in maps if m.get("series") == "A: Ucrit(T0), physical, Pi=0" and m["status"] == "PASS"
                and m["T0"] <= T0_VALID], key=lambda m: m["T0"])
    B = sorted([m for m in maps if m.get("series") == "B: T0crit(U), physical, Pi=0" and m["status"] == "PASS"
                and m["U"] <= U_VALID and m["T0_crit"] <= T0_VALID], key=lambda m: m["U"])
    FA = [m for m in maps if m.get("g_discontinuity_at_root") and m.get("series") == "A: Ucrit(T0), physical, Pi=0"]
    FB = [m for m in maps if m.get("g_discontinuity_at_root") and m.get("series") == "B: T0crit(U), physical, Pi=0"]
    xa, ya = np.array([m["T0"] for m in A]), np.array([m["U_crit"] for m in A])
    ax.fill_between(xa, 0, ya, color="C3", alpha=0.12, lw=0)
    ax.fill_between(xa, ya, U_VALID, color="C0", alpha=0.12, lw=0)
    ax.axvspan(T0_VALID, 680, color="0.82", zorder=0)
    ax.axhspan(U_VALID, 750, xmax=(T0_VALID - 600) / 80, color="0.82", zorder=0)
    ax.plot(xa, ya, "-k", lw=1.4, label="A: U_crit(T0) roots at fixed T0 (validated)")
    ax.plot([m["T0_crit"] for m in B], [m["U"] for m in B], "o", mfc="none", mec="C3", ms=5,
            label="B: T0,crit(U) roots at fixed U (validated)")
    ax.plot([m["T0"] for m in FA] + [m["sign_change_location"] for m in FB],
            [m["sign_change_location"] for m in FA] + [m["U"] for m in FB], "x", color="0.35", ms=5,
            label="Flagged sign changes (criterion discontinuity):\ndiagnostic only, not boundary roots")
    ax.plot([625], [373], "*", color="C2", ms=9, label="base point (T0 = 625 K, U = 373)")
    ax.text(603, 120, "g > 0: inlet perturbations amplified", fontsize=7)
    ax.text(603, 420, "g < 0: attenuated", fontsize=7)
    ax.text(603, 600, "U > 450: outside validated domain (not interpreted)", fontsize=7, color="0.3")
    ax.text(652, 300, "criterion discontinuity /\nboundary not interpreted", fontsize=7, color="0.3", rotation=90)
    ax.set(xlim=(600, 680), ylim=(0, 750), xlabel="feed temperature T0 (K)", ylabel="U (W m$^{-2}$ K$^{-1}$)",
           title="Sensitivity-based parametric runaway boundary map\n(Tc = 625 K, Pi = 0; validated domain T0 <= 650 K, U <= 450)")
    ax.legend(loc="upper left", fontsize=6.3, bbox_to_anchor=(1.01, 1.0))
    save(fig, "fig08_boundary_map.png")

    # 9 hot-spot location
    fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.4))
    r = series(t0, "physical Pi=0")
    ax[0].plot([q["inputs"]["T0"] for q in r], [q["z_hs"] for q in r], "-", lw=1.1)
    mark_side(ax[0], [q["inputs"]["T0"] for q in r], [q["z_hs"] for q in r], [q["g"] for q in r])
    ax[0].set(xlabel="T0 (K)", ylabel="hot-spot location z_hs (m)", title="z_hs vs T0 (U = 373, Pi = 0)")
    r = sorted([x for x in us if x.get("series") == "physical Pi=0" and x["status"] == "PASS"], key=lambda x: x["inputs"]["U"])
    ax[1].plot([q["inputs"]["U"] for q in r], [q["z_hs"] for q in r], "-", lw=1.1)
    mark_side(ax[1], [q["inputs"]["U"] for q in r], [q["z_hs"] for q in r], [q["g"] for q in r])
    ax[1].set(xlabel="U (W m$^{-2}$ K$^{-1}$)", ylabel="z_hs (m)", title="z_hs vs U (T0 = 625 K, Pi = 0)")
    save(fig, "fig09_hotspot_location.png", FOOT + " Filled: g > 0; open: g < 0.")

    # 10 pressure decomposition
    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.0))
    ref = [m for m in pb if m.get("series") == "A constant level" and m.get("level") == 1.0][0]["U_crit"]
    for s, mk, lab in (("C outlet combined", "^", "C: outlet-fixed P(L) = P_ref (level + decline)"),
                       ("A' constant at matched P0", "s", "A': constant P at matched inlet P0 of C (level only)"),
                       ("B inlet decline", "v", "B: inlet-fixed P(0) = P_ref (decline only)")):
        r = sorted([m for m in pb if m.get("series") == s], key=lambda m: m["Pi"])
        ok = [m for m in r if m["status"] == "PASS"]
        ln, = ax[0].plot([m["Pi"] for m in ok], [m["U_crit"] for m in ok], marker=mk, ms=4.5, label=lab)
        for m in [m for m in r if "above the searched range" in m.get("outcome", "")]:
            ax[0].plot(m["Pi"], 775, marker=mk, mfc="none", color=ln.get_color(), clip_on=False)
            ax[0].annotate("> 750", (m["Pi"], 775), textcoords="offset points", xytext=(5, -3), fontsize=7)
    ax[0].axhline(ref, ls="--", color="0.5", lw=0.8, label=f"Pi = 0, P = P_ref ({ref:.2f})")
    ax[0].set(xlabel="Pi (-)", ylabel="U_crit (W m$^{-2}$ K$^{-1}$)", ylim=(300, 800), title="(a) U_crit vs Pi (T0 = 625 K)")
    ax[0].legend(fontsize=6.3)
    r = sorted([m for m in pb if m.get("series") == "A constant level"], key=lambda m: m["level"])
    ok = [m for m in r if m["status"] == "PASS"]
    ax[1].plot([m["level"] for m in ok], [m["U_crit"] for m in ok], "o-", ms=4, label="A: constant P = P0 (level only)")
    for m in [m for m in r if "above the searched range" in m.get("outcome", "")]:
        ax[1].plot(m["level"], 775, "o", mfc="none", color="C0", clip_on=False)
        ax[1].annotate("> 750", (m["level"], 775), textcoords="offset points", xytext=(5, -3), fontsize=7)
    ax[1].set(xlabel="constant pressure level P0 / P_ref (-)", ylabel="U_crit (W m$^{-2}$ K$^{-1}$)", ylim=(300, 800),
              title="(b) A: pressure-level effect only (T0 = 625 K)")
    ax[1].legend(fontsize=6.3)
    dec = [d for d in read_json("pressure_decomposition") if d["point"] == "U=420"]
    xx = np.arange(len(dec))
    w = 0.27
    for i, (k, lab) in enumerate((("level_effect(A-baseline)", "level: A - baseline"),
                                  ("decline_effect(C-A)", "decline: C - A"), ("combined(C-baseline)", "combined: C - baseline"))):
        ax[2].bar(xx + (i - 1) * w, [d["T_max"][k] for d in dec], w, label=lab)
    ax[2].set_xticks(xx, [f"Pi = {d['Pi']:g}" for d in dec])
    ax[2].set(ylabel="change in calculated T_max (K)", title="(c) T_max decomposition at U = 420, T0 = 625 K")
    ax[2].legend(fontsize=6.3)
    save(fig, "fig10_pressure_level_vs_decline.png",
         "A = constant pressure (pressure-level effect only); B = inlet-fixed pressure decline; C = outlet-fixed combined "
         "pressure-level + pressure-decline effect; A' = constant pressure at the matched inlet pressure of C.\n"
         "At fixed outlet pressure, the modeled pressure-level effect exceeds the pressure-decline effect over the tested "
         "Pi range (property of the present model/parameterization; absolute pressure drop not validated). "
         "Calculated runaway-branch T_max values lie outside the documented kinetic validity range.")

    # 11 O2 comparison
    fig, ax = plt.subplots(1, 3, figsize=(12.6, 3.7))
    for s, st in (("fixed-O2 (base)", "-"), ("O2-depletion variant", "--")):
        r = series(o2, s)
        ax[0].plot([q["inputs"]["T0"] for q in r], [q["T_max"] for q in r], st, lw=1.1, label=s)
        ax[1].plot([q["inputs"]["T0"] for q in r], [q["X_L"] for q in r], st, lw=1.1, label=s)
        rr = sorted([m for m in o2all if m["series"] == s and m["status"] == "PASS"], key=lambda m: m["T0"])
        ax[2].plot([m["T0"] for m in rr], [m["U_crit"] for m in rr], st, marker="o", ms=4, label=s)
        for m in rr:
            ax[2].annotate(f"{m['U_crit']:.1f}", (m["T0"], m["U_crit"]), textcoords="offset points",
                           xytext=(3, 4 if s.startswith("fixed") else -10), fontsize=6.5)
    ax[0].set(xlabel="T0 (K)", ylabel="calculated T_max (K)", title="T_max (U = 373, Pi = 0)")
    ax[1].set(xlabel="T0 (K)", ylabel="X(L) (-)", title="X(L) (U = 373, Pi = 0)")
    ax[2].set(xlabel="T0 (K)", ylabel="U_crit (W m$^{-2}$ K$^{-1}$)", title="U_crit (Tc = 625 K, Pi = 0)")
    ax[0].legend()
    save(fig, "fig11_O2_model_form.png", FOOT + " O2-depletion is a model-form check, not the base model.")

    # 12 sensitivity ranking
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.6))
    ps = [r for r in rob["point_sensitivities"] if r["point"].startswith("P2")]
    def lab(r, key):
        n = "Pi (outlet 0.2; runaway branch)" if r["parameter"].startswith("Pi") else r["parameter"]
        return n + (" *" if r[key]["stable_step"] is None else "")
    for a, key, title in ((ax[0], "dT_hs", "local: dT_hs = T_max - Tc at P2 (U = 420, g < 0)"), (ax[1], "X_L", "local: X(L) at P2")):
        vals = sorted([(lab(r, key), r[key]["normalized"]) for r in ps], key=lambda t: abs(t[1]))
        a.barh([v[0] for v in vals], [v[1] for v in vals], color=["C3" if v[1] > 0 else "C0" for v in vals])
        a.set(xlabel="normalized sensitivity (phi/y) dy/dphi", title=title)
    bs = sorted(rob["boundary_sensitivities"].items(), key=lambda t: abs(t[1]["normalized"]))
    ax[2].barh([k + (" *" if v["stable_step"] is None else "") for k, v in bs], [v["normalized"] for _, v in bs],
               color=["C3" if v["normalized"] > 0 else "C0" for _, v in bs])
    ax[2].set(xlabel="normalized sensitivity (phi/U_crit) dU_crit/dphi", title="local: U_crit at T0 = 625 K")
    save(fig, "fig12_sensitivity_ranking.png", FOOT + " Red: positive; blue: negative. * = no stable finite-difference "
         "step (smallest-step estimate, not validated).")
    print("figures written:", sorted(p.name for p in FIG.glob("*.png")))


if __name__ == "__main__":
    main()
