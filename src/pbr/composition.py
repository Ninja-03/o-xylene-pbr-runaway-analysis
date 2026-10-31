"""(4) Composition / speciation.

Reaction (Phase 2): C8H10 + nu_O2 O2 -> C8H4O3 + 3 H2O (nu_O2 from the file; the H2O coefficient
follows from the hydrogen balance and the full C/H/O balance is checked). Total moles unchanged.
Conversion X = 1 - N_x/N_x0 (molar-flow basis); with N_total = N_f: y_x = y_x0 (1 - X).
FIXED-O2 ASSUMPTION (Phase 2, model.oxygen_depletion = false): the O2 dependence enters the rate only
through P/P_ref with y_O2 held at its inlet value y_O2_0. O2 depletion is NOT modelled; the helpers
marked 'diagnostic' quantify the neglected effect and are never used in the rate.
"""
from . import properties

FORMULA = {"o-xylene": (8, 10, 0), "O2": (0, 0, 2), "PA": (8, 4, 3), "H2O": (0, 2, 1)}  # (C, H, O)


def stoichiometry(p):
    nu_h2o = (FORMULA["o-xylene"][1] - FORMULA["PA"][1]) / FORMULA["H2O"][1]
    nu = {"o-xylene": -1.0, "O2": -p.nu_O2, "PA": 1.0, "H2O": nu_h2o}
    for e, sym in enumerate("CHO"):
        bal = sum(n * FORMULA[k][e] for k, n in nu.items())
        if abs(bal) > 1e-12:
            raise ValueError(f"element balance for {sym} does not close: {bal}")
    return nu


def y_x(X, p):
    return p.y_x0 * (1.0 - X)


def c_x(X, T, P, p):
    """o-xylene concentration y_x P/(R T) [mol/m3]."""
    return y_x(X, p) * properties.molar_concentration(T, P, p)


def y_O2_in_rate(p):
    """O2 mole fraction implied by the rate law: fixed at its inlet value (lumped into k_m at P_ref)."""
    return p.y_O2_0


def molar_flows(X, p):
    """Per-tube molar flows [mol/s] at conversion X. 'O2_if_consumed' is DIAGNOSTIC only."""
    nu, ext = stoichiometry(p), p.N_x0 * X
    return {"o-xylene": p.N_x0 + nu["o-xylene"] * ext, "PA": nu["PA"] * ext, "H2O": nu["H2O"] * ext,
            "O2_if_consumed": p.y_O2_0 * p.N_f + nu["O2"] * ext, "total": p.N_f}


def o2_depletion_fraction(X, p):
    """DIAGNOSTIC only: fraction of inlet O2 that stoichiometric consumption would remove at X."""
    return p.nu_O2 * p.y_x0 * X / p.y_O2_0
