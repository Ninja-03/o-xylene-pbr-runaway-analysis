#!/usr/bin/env python3
"""
Phase 2 checker v0.2 - provenance, dimensional consistency and plausibility.

Usage:  python checks/phase2_plausibility.py [path/to/params.yaml]

Reads ONLY the parameter file; integrates NO reactor profile (Phase 3).
Check status:
  PASS   - satisfied
  INFO   - evidence recorded, no action needed
  FLAG   - documented limitation or expected effect; resolution printed
  DECIDE - open modelling decision            -> exit code 1
  FAIL   - hard error (provenance/units/dims) -> exit code 1
"""
import math
import pathlib
import sys

import numpy as np
import yaml

HERE = pathlib.Path(__file__).resolve().parent
PARAMS = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "params" / "base_case.yaml"

# --------------------------------------------------------------- dimensions
BASE = ("M", "L", "T", "Th", "N")


class D:
    def __init__(self, **e):
        self.e = tuple(float(e.get(k, 0.0)) for k in BASE)

    @classmethod
    def of(cls, exps):
        d = cls()
        d.e = tuple(exps)
        return d

    def __mul__(self, o):
        return D.of(a + b for a, b in zip(self.e, o.e))

    def __truediv__(self, o):
        return D.of(a - b for a, b in zip(self.e, o.e))

    def __pow__(self, n):
        return D.of(a * n for a in self.e)

    def __eq__(self, o):
        return all(abs(a - b) < 1e-12 for a, b in zip(self.e, o.e))

    def dimless(self):
        return all(abs(a) < 1e-12 for a in self.e)

    def __repr__(self):
        s = " ".join(f"{k}^{v:g}" for k, v in zip(BASE, self.e) if abs(v) > 1e-12)
        return s or "1"


class Q:
    """Value with dimension; adding unlike dimensions raises TypeError."""

    def __init__(self, v, d=None):
        self.v = float(v)
        self.d = d if d is not None else D()

    @staticmethod
    def c(o):
        return o if isinstance(o, Q) else Q(o)

    def __mul__(self, o):
        o = Q.c(o)
        return Q(self.v * o.v, self.d * o.d)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = Q.c(o)
        return Q(self.v / o.v, self.d / o.d)

    def __rtruediv__(self, o):
        return Q.c(o) / self

    def __pow__(self, n):
        return Q(self.v ** n, self.d ** n)

    def __add__(self, o):
        o = Q.c(o)
        if not self.d == o.d:
            raise TypeError(f"inhomogeneous sum [{self.d}] + [{o.d}]")
        return Q(self.v + o.v, self.d)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.v, self.d)

    def __sub__(self, o):
        return self + (-Q.c(o))

    def __rsub__(self, o):
        return Q.c(o) - self


def qexp(x):
    if not x.d.dimless():
        raise TypeError(f"exp() of dimensional quantity [{x.d}]")
    return Q(math.exp(x.v))


JOULE, WATT, PASCAL = D(M=1, L=2, T=-2), D(M=1, L=2, T=-3), D(M=1, L=-1, T=-2)
UNITS = {
    "-": D(), "K": D(Th=1), "m": D(L=1), "Pa": PASCAL, "1/s": D(T=-1), "kg/s": D(M=1, T=-1),
    "kg/mol": D(M=1, N=-1), "m^3/mol": D(L=3, N=-1), "J/mol": JOULE / D(N=1),
    "J/(kg K)": JOULE / D(M=1, Th=1), "J/(mol K)": JOULE / D(N=1, Th=1),
    "W/(m^2 K)": WATT / D(L=2, Th=1),
}
CONV = {  # unit definitions only (no physical data)
    ("kJ/kmol", "J/mol"): 1.0, ("kcal/kmol", "J/mol"): 4.1868, ("kJ/(m^2 s K)", "W/(m^2 K)"): 1000.0,
    ("kJ/(kg K)", "J/(kg K)"): 1000.0, ("kN/m^2", "Pa"): 1000.0, ("1/sec", "1/s"): 1.0, ("m", "m"): 1.0,
}
ALLOWED = ("VERIFIED", "LITERATURE", "PROPERTY", "PROJECT-DEFINED")
PEER = ("primary", "peer_reviewed")
FULL = ("full_text", "code")
RESULTS = []


def record(cid, status, text, resolution=""):
    RESULTS.append((cid, status, text, resolution))


def finish():
    print("\n" + "=" * 100 + "\nCHECK RESULTS  (" + str(PARAMS) + ")\n" + "=" * 100)
    for cid, st, text, res in RESULTS:
        print(f"[{st:6s}] {cid:5s} {text}")
        if res:
            print(f"{'':15s}-> resolution: {res}")
    counts = {s: sum(1 for r in RESULTS if r[1] == s) for s in ("PASS", "INFO", "FLAG", "DECIDE", "FAIL")}
    print("-" * 100)
    print(f"{len(RESULTS)} checks: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    blocked = counts["FAIL"] + counts["DECIDE"]
    print("GATE: " + ("OPEN - no FAIL, no DECIDE" if not blocked else f"CLOSED - {blocked} FAIL/DECIDE item(s)"))
    sys.exit(1 if blocked else 0)


cfg = yaml.safe_load(PARAMS.read_text())
P, REFS, META = cfg.get("parameters", {}), cfg.get("references", {}), cfg.get("meta", {})


def ref_info(key):
    r = REFS.get(key)
    if isinstance(r, dict):
        return r.get("kind", "unknown"), r.get("accessed", "unknown")
    return "unknown", "unknown"  # free-text references carry no access record


# ------------------------------------------------------ C1 provenance (strict)
problems, bad_status = [], {}
critical = META.get("critical")
if not critical:
    problems.append("meta.critical (list of critical physical parameters) is missing")
for name, p in P.items():
    v, st = p.get("value"), p.get("status")
    if not isinstance(v, (int, float)) or isinstance(v, bool):
        problems.append(f"{name}: non-numeric value {v!r}")
    if st not in ALLOWED:
        bad_status[st] = bad_status.get(st, 0) + 1
        continue
    srcs = list(p.get("sources", []))
    for s in srcs + list(p.get("range_sources", [])):
        if s not in REFS:
            problems.append(f"{name}: unknown source key {s}")
    infos = [ref_info(s) for s in srcs if s in REFS]
    if "recall" in str(p.get("note", "")).lower():
        problems.append(f"{name}: note marks the value as recalled")
    if st == "VERIFIED" and not any(k in PEER and a == "full_text" for k, a in infos):
        problems.append(f"{name}: VERIFIED but no accessed primary/peer-reviewed full text")
    if st == "LITERATURE" and not any(a in FULL for _, a in infos):
        problems.append(f"{name}: LITERATURE but no accessed full-text/code source")
    if st == "PROPERTY" and not srcs:
        problems.append(f"{name}: PROPERTY without source")
    if st == "PROJECT-DEFINED" and len(str(p.get("note", ""))) < 40:
        problems.append(f"{name}: PROJECT-DEFINED without written rationale")
    if srcs and all(a == "no" for _, a in infos):
        problems.append(f"{name}: supported only by sources that were not accessed")
    if p.get("unit") not in UNITS:
        problems.append(f"{name}: unit {p.get('unit')!r} not in dimension table")
    if "range" in p:
        lo, hi = p["range"]
        if not (lo < hi and lo <= v <= hi):
            problems.append(f"{name}: malformed range or base value outside range")
if bad_status:
    problems.insert(0, "statuses outside {" + ", ".join(ALLOWED) + "}: " +
                    ", ".join(f"{k} x{n}" for k, n in bad_status.items()))
for c in critical or []:
    if c not in P:
        problems.append(f"critical parameter {c} missing")
    elif P[c].get("status") not in ("VERIFIED", "LITERATURE"):
        problems.append(f"critical parameter {c} has status {P[c].get('status')}; must be VERIFIED or LITERATURE")
n_by = {s: sum(1 for p in P.values() if p.get("status") == s) for s in ALLOWED}
record("C1", "FAIL" if problems else "PASS", "provenance: " + (
    "; ".join(problems) if problems else
    f"{len(P)} parameters ({', '.join(f'{k} {v}' for k, v in n_by.items())}); "
    f"all {len(critical)} critical parameters VERIFIED/LITERATURE from accessed sources; no recalled values"))
if problems:
    finish()

# ----------------------------------------- C2 unit conversions and re-derivations
AUD = cfg["source_audit"]
lines, bad = [], []
for name, p in P.items():
    for key, tol in (("original", 1e-4), ("cross_check", 5e-3)):
        o = p.get(key)
        if not o:
            continue
        f = CONV.get((o["unit"], p["unit"]))
        if f is None:
            bad.append(f"{name}: no conversion {o['unit']} -> {p['unit']}")
            continue
        si = o["value"] * f
        rel = abs(si - p["value"]) / abs(si)
        lines.append(f"{name}[{key}, {o.get('source', '')}] {rel:.1e}")
        if rel > tol:
            bad.append(f"{name}: {key} {o['value']} {o['unit']} -> {si:.6g}, stored {p['value']:.6g}")
m = AUD["re_code_Mwf_formula"]
Mwf_re = (m["x_air"] * (m["y_N2"] * m["M_N2"] + m["y_O2"] * m["M_O2"]) + m["x_oX"] * m["M_oX"]) * 1e-3
if abs(Mwf_re - P["M_wf"]["value"]) > 1e-9:
    bad.append(f"M_wf re-derived {Mwf_re:.7f} != stored")
if abs(m["x_air"] * m["y_O2"] - P["y_O2_0"]["value"]) > 1e-9:
    bad.append("y_O2_0 != x_air*y_O2 of RE_CODE feed")
nu = P["nu_O2"]["value"]  # C8H10 + nu O2 -> C8H4O3 + 3 H2O
if not (8 == 8 and 10 == 4 + 2 * 3 and abs(2 * nu - (3 + 3)) < 1e-12):
    bad.append("atom balance C8H10 + nu O2 -> C8H4O3 + 3 H2O does not close")
Vm, MoX = P["Vm_normal"]["value"], P["M_oX"]["value"]
y_of = lambda g: (g / 1000 / MoX) / (1 / Vm + g / 1000 / MoX)
y_lo, y_hi = (y_of(g) for g in AUD["loading_range_g_per_Nm3"])
rlo, rhi = P["y_x0"]["range"]
if abs(y_lo - rlo) / rlo > 0.01 or abs(y_hi - rhi) / rhi > 0.01:
    bad.append(f"y_x0 range {rlo}-{rhi} != loadings -> {y_lo:.4f}-{y_hi:.4f}")
record("C2", "FAIL" if bad else "PASS", "conversions & re-derivations: " + (
    "; ".join(bad) if bad else "; ".join(lines) + f"; M_wf = {Mwf_re * 1e3:.4f} kg/kmol re-derived; "
    f"30-100 g/Nm3 -> y = {y_lo:.4f}-{y_hi:.4f}; C/H/O atom balance closes with nu_O2 = {nu:g}"))

# ------------------------------------------------------------- C3 dimensions
q = {n: Q(p["value"], UNITS[p["unit"]]) for n, p in P.items()}
R, K = q["R_gas"], D(Th=1)
A_c = math.pi * q["R_tube"] ** 2
N_f = q["Qrho"] / q["M_wf"]


def k_of(T):
    return q["k_m"] * qexp(q["E_over_R"] * (1 / q["T_m"] - 1 / T))


def groups(T0, Tc, U, y, cp):
    gamma = q["E_over_R"] / Tc
    dTad = (-q["dH_rxn"]) * y / (q["M_wf"] * cp)
    B = gamma * dTad / Tc
    Da = A_c * k_of(Tc) * q["P_ref"] * q["L"] / (N_f * R * Tc)
    St = 2 * math.pi * q["R_tube"] * U * q["L"] / (q["Qrho"] * cp)
    return dict(gamma=gamma, dTad=dTad, B=B, Da=Da, St=St, NS=St / (B * Da),
                theta0=gamma * (T0 - Tc) / Tc, lam=q["nu_O2"] * y / q["y_O2_0"])


BASE_ARGS = dict(T0=q["T0"], Tc=q["Tc"], U=q["U"], y=q["y_x0"], cp=q["cp"])
try:
    g = groups(**BASE_ARGS)
    bad = [n for n in ("gamma", "B", "Da", "St", "NS", "theta0", "lam") if not g[n].d.dimless()]
    if not g["dTad"].d == K:
        bad.append("dT_ad not a temperature")
    T, Pp, y = q["T0"], q["P_ref"], q["y_x0"]
    r = k_of(T) * (Pp / q["P_ref"]) * (y * Pp / (R * T))
    dXdz = A_c * r / (y * N_f)
    dTdz = ((-q["dH_rxn"]) * A_c * r - 2 * math.pi * q["R_tube"] * q["U"] * (T - q["Tc"])) / (q["Qrho"] * q["cp"])
    d_p, eps = Q(1.0, D(L=1)), Q(0.5)  # placeholders: only dimensions are tested here
    G = q["Qrho"] / A_c
    dPdz = -(1.75 * G ** 2 * (1 - eps) / ((Pp * q["M_wf"] / (R * T)) * d_p * eps ** 3))
    Pi_f = 1.75 * G ** 2 * (1 - eps) * R * q["Tc"] * q["L"] / (q["M_wf"] * d_p * eps ** 3 * q["P_ref"] ** 2)
    for n_, (v_, d_) in {"dX/dz": (dXdz, D(L=-1)), "dT/dz": (dTdz, K / D(L=1)),
                         "dP/dz": (dPdz, PASCAL / D(L=1)), "Pi": (Pi_f, D())}.items():
        if not v_.d == d_:
            bad.append(f"{n_} is [{v_.d}], expected [{d_}]")
    record("C3", "FAIL" if bad else "PASS", "dimensions: " + ("; ".join(bad) if bad else
           "mole, energy and Ergun balances homogeneous term by term; gamma, B, Da, St, N/S, "
           "theta0, lambda, Pi dimensionless; dT_ad in K"))
except TypeError as exc:
    record("C3", "FAIL", f"dimensional inconsistency: {exc}")
    finish()

# --------------------------------------------------- C4-C7 regime plausibility
Da, B, St, NS, gam, dTad = (g[k].v for k in ("Da", "B", "St", "NS", "gamma", "dTad"))
record("C4", "PASS" if 0.05 < Da < 20 else "FLAG",
       f"kinetic regime: Da = {Da:.4f}; isothermal estimate X(L) = {1 - math.exp(-Da):.2f}")
record("C5", "PASS" if gam > 10 and B > 5 else "FLAG",
       f"thermal regime: gamma = {gam:.3f}, dT_ad = {dTad:.1f} K, B = {B:.2f}")
record("C6", "PASS" if 1 < NS < 10 else "FLAG",
       f"proximity to runaway: St = {St:.3f}, N/S = St/(B Da) = {NS:.3f} (order unity)")
T_run = AUD["re_documented_runaway_feed_T_K"]
lo, hi = P["T0"]["range"]
rows = []
for T_ in np.linspace(lo, hi, 5):
    gg = groups(**{**BASE_ARGS, "T0": Q(T_, K), "Tc": Q(T_, K)})
    rows.append(f"{gg['NS'].v:.2f}@{T_:.0f}K")
ns_vals = [float(s.split("@")[0]) for s in rows]
ok = lo < T_run < hi and min(ns_vals) < 1 < max(ns_vals)
record("C7", "PASS" if ok else "FLAG",
       f"T0 range {lo:.0f}-{hi:.0f} K: physical case (Tc = {q['Tc'].v:.0f} K) contains the documented runaway "
       f"feed temperature ~{T_run:.0f} K; normalized case (Tc = T0) N/S = " + ", ".join(rows))
ylo, yhi = P["y_x0"]["range"]
yrows = []
for yy in (ylo, P["y_x0"]["value"], yhi):
    gg = groups(**{**BASE_ARGS, "y": Q(yy)})
    yrows.append(f"y={yy:.4f}: dT_ad {gg['dTad'].v:.0f} K, B {gg['B'].v:.1f}, N/S {gg['NS'].v:.2f}")
Ulo, Uhi = P["U"]["range"]
ns_Uhi = groups(**{**BASE_ARGS, "U": Q(Uhi, UNITS['W/(m^2 K)'])})["NS"].v
record("C7b", "PASS", "study ranges: " + "; ".join(yrows) + f"; U {Ulo:.0f}-{Uhi:.0f} W/(m2 K) -> N/S 0-{ns_Uhi:.2f}")

# --------------------------------------------------- C8-C10 pressure specification
Pref = P["P_ref"]["value"]
Pi_max, Pi_in = P["Pi_max"]["value"], P["Pi_inlet_max"]["value"]
PiP0 = lambda P0: ((P0 / Pref) ** 2 - 1) / 2
dp_bed = AUD["US4397768_bed_dP_m_water"] * AUD["m_water_Pa"]
gauge = [b * 1e5 for b in AUD["US7557223_inlet_gauge_bar"]]
anchors = [PiP0(Pref + dp_bed)] + [PiP0(AUD["standard_atmosphere_Pa"] + gg) for gg in gauge]
p0_max = math.sqrt(1 + 2 * Pi_max)
prof = ", ".join(f"{math.sqrt(1 + 2 * Pi_max * (1 - z)):.3f}" for z in (0, 0.25, 0.5, 0.75, 1.0))
record("C8", "PASS" if Pi_max >= max(anchors) - 1e-3 else "FAIL",
       f"outlet-specified Ergun (P(L) = {Pref / 1e3:.1f} kPa): documented anchors give Pi = "
       f"{anchors[0]:.3f} (2 m water-gauge bed dP) and {anchors[1]:.3f}-{anchors[2]:.3f} (0.35-0.55 bar gauge "
       f"inlet); study range Pi = 0-{Pi_max:.2f} -> isothermal P0 = 1.000-{p0_max:.3f} P_out "
       f"({Pref * p0_max / 1e3:.1f} kPa max); P(z)/P_out at z/L = 0,.25,.5,.75,1: {prof}; P >= P_out > 0 everywhere")
pin = math.sqrt(1 - 2 * Pi_in)
record("C9", "PASS" if pin > 3 * P["P_min_frac"]["value"] else "FAIL",
       f"inlet-specified comparison (P(0) = P_ref) limited to Pi <= {Pi_in:.2f}: isothermal P_out/P0 = {pin:.3f}; "
       f"1-atm inlet at Pi = {Pi_max:.2f} would be infeasible (1 - 2 Pi = {1 - 2 * Pi_max:.2f} < 0)",
       "runs that trip the P < P_min guard in deep runaway are discarded as invalid, never read as runaway")
record("C10", "FLAG",
       f"pressure level effect: outlet-specified Pi = {Pi_max:.2f} raises inlet rate x{p0_max ** 2:.2f} (rate ~ P^2), "
       f"inlet N/S ~ {NS / p0_max ** 2:.2f} vs {NS:.2f} without pressure drop",
       "expected physics, not an error: M3 separates pressure-level and pressure-decline effects "
       "(Ergun/outlet spec vs constant P = P0 vs constant P = P_ref)")

# ----------------------------------------------------- C11-C14 documented limits
lam = g["lam"].v
record("C11", "FLAG", f"O2 held at inlet mole fraction (benchmark pseudo-first-order): lambda = {lam:.3f} -> "
       f"O2 depletion up to {100 * lam:.0f}% at X = 1, ~{50 * lam:.0f}% at X = 0.5",
       "documented benchmark assumption (RE_CODE, ZB2020); M8 model-form check with O2 depletion on")
record("C12", "FLAG", f"adiabatic ceiling T0 + dT_ad = {q['T0'].v + dTad:.0f} K, far above documented operation "
       f"(salt 618-658 K; observed hot spot ~713 K); kinetic validity range not documented",
       "criterion locates the boundary before the runaway branch; runaway-branch temperatures reported "
       "only qualitatively")
record("C13", "FLAG", "no documented particle/voidage data for the adopted dataset: U attainability and the "
       "viscous Ergun term cannot be evaluated",
       "U above the documented 373 W/(m2 K) labelled 'hypothetical cooling capacity'; pressure drop studied "
       "through Pi (inertia-dominated Ergun form) anchored to documented industrial pressures")
diffs = AUD["zb2020_differences"]
record("C14", "FLAG", "the two accessed reproductions differ in L, U, c_p, mass flux and feed: " +
       ", ".join(f"{k}: RE {v[0]} vs ZB {v[1]}" for k, v in diffs.items()),
       "coherent RE_CODE set adopted; ZB2020 used only for dimensionally unambiguous values that agree "
       "(E/R, dH, D_t, T_ref); c_p spanned 992-1046 J/(kg K) in M8")

# ----------------------------------------------------------- S source audit
t = AUD["re_table_6_5_printed"]
dTad_tab = t["dH_kJ_per_kmol"] * P["y_x0"]["value"] / (P["M_wf"]["value"] * 1e3 * P["cp"]["value"] / 1e3)
Da_tab = Da * t["km_per_s"] / P["k_m"]["value"]
record("S1", "PASS", f"R&E Table 6.5 misprints: printed dH gives dT_ad = {abs(dTad_tab):.2f} K and printed "
       f"k_m gives Da = {Da_tab:.0f} (complete conversion in the first cm) - both contradict the published profiles",
       "values from the code that generated the figures (k_m = 2.0822 1/s, dH = -1.284e6 kJ/kmol) adopted")
z = AUD["zb2020_table3_literal"]
k1 = z["P_O2_kPa"] * math.exp(z["lnA"] - P["E_over_R"]["value"] / P["T_m"]["value"])  # kmol/(h kg kPa), as printed
kv = z["rho_b"] * (k1 * 1000 / 3600 / 1000) * P["R_gas"]["value"] * P["T_m"]["value"]  # 1/s, volumetric
rhoL_geo = z["rho_b"] * math.pi * (z["D_t"] / 2) ** 2
record("S2", "PASS", f"ZB2020 Table 3 read literally (k in kPa units) gives a volumetric k = {kv:.3g} 1/s at 625 K "
       f"({kv / P['k_m']['value']:.0f}x RE_CODE); its linear density {z['rho_L']} kg/m != rho_b*A_t = "
       f"{rhoL_geo:.4f} kg/m",
       "ZB2020 pre-exponential not adopted; only E/R and dH (unit-unambiguous, agreeing with RE_CODE) used")
p_imp = t["km_per_s"] / (z["rho_b"] * math.exp(z["lnA"] - 13636 / 625) * 0.0820574 * 625)
record("S3", "INFO", f"printed R&E k_m = 1922.6 equals rho_b*exp(19.837-13636/625)*p_O2*R*T on an atm, per-hour "
       f"basis with p_O2 = {p_imp:.4f} atm (air O2 fraction) - consistent with an atm-based original constant and a "
       f"unit misprint; recorded, not used")

# ------------------------------------------------ C15-C17 criteria and decisions
C, M = cfg.get("criteria", {}), cfg.get("model", {})
need = ("hot_spot", "runaway_boundary", "Ucrit", "light_off", "engineering_T_limit", "dynamic_stability")
miss = [k for k in need if k not in C or not C[k].get("definition")]
rb = C.get("runaway_boundary", {})
rb_ok = any(ref_info(s)[1] == "full_text" for s in rb.get("sources", []))
record("C15", "PASS" if not miss and rb_ok else "DECIDE", "criteria: " + (
    f"missing {miss}" if miss else "hot spot, runaway boundary, U_crit, light-off, engineering limit and dynamic "
    "stability defined separately; runaway boundary taken from an accessed peer-reviewed full text"))
lim = str(C.get("engineering_T_limit", {}).get("definition", "")).lower()
record("C16", "PASS" if lim.startswith("none") else "DECIDE",
       "engineering temperature limit: none used (no sourced catalyst/material limit; 440 C is an observed hot spot)")
dec = [k for k in ("reaction", "kinetic_law", "rate_basis", "coordinate", "pressure_drop",
                   "coolant_configurations") if k not in M]
record("C17", "PASS" if not dec else "DECIDE",
       "model decisions recorded: " + ("all" if not dec else f"missing {dec}") +
       " (kinetic law, basis, coordinate, pressure spec, both coolant configurations)")
finish()
