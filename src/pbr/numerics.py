"""Numerical settings (NOT physical parameters).

Solver designation: BDF is PRIMARY; Radau and LSODA are CROSS-CHECK solvers (Phase 3 instruction). The
parameter file's numerics block was revised in Phase 3 to say the same (it previously named Radau) and
params.check_spec() stops if the two ever disagree.
rtol and atol_scaled come from the numerics block (revised in Phase 3 from 1e-8/1e-10 to 1e-10/1e-12 after
V5/V6/V7/V9 exposed error amplification near runaway); the absolute tolerance vector
is atol_scaled times the natural scale of each state. Guard and shooting tolerances are tied to rtol
so that they always sit just above the integration accuracy.
"""
import numpy as np

PRIMARY = "BDF"
CROSSCHECK = ("Radau", "LSODA")
N_GRID = 1501          # dense-output evaluation points for stored profiles
X_BOUND_REL = 10.0     # 0 <= X <= 1 is checked with allowance X_BOUND_REL * rtol
SHOOT_SCAN = 9         # coarse scan: one sign change and monotonic P(L)(P0) before Brent
SHOOT_XTOL_REL = 0.1   # Brent tolerance on P0 = SHOOT_XTOL_REL * rtol * P_ref
SHOOT_RES_REL = 10.0   # accepted |P(L)/P_ref - 1| <= SHOOT_RES_REL * rtol
BRACKET_MARGIN_REL = 100.0  # analytic shooting bracket widened by BRACKET_MARGIN_REL * rtol (integration error at the ends)
BRENT_MAXITER = 200


def tolerances(p, rtol=None, atol_scaled=None):
    r = float(p.numerics["rtol"]) if rtol is None else float(rtol)
    a = float(p.numerics["atol_scaled"]) if atol_scaled is None else float(atol_scaled)
    return r, a


def x_bound_tol(rtol):
    return X_BOUND_REL * rtol


def state_scales(p):
    """Natural scales of [X, T, P, Q]: 1, T_c, P_ref, Qrho c_p T_c."""
    return np.array([1.0, p.Tc, p.P_ref, p.Qrho * p.cp * p.Tc])
