"""(3) Kinetics - the settled Phase 2 law, implemented exactly (no reduction to a generic 1st-order form):

    r = k_m exp[-(E/R)(1/T - 1/T_m)] (P/P_ref) c_x,   c_x = y_x P/(R T),   y_x = y_x0 (1 - X)

Units: k_m [1/s], c_x [mol/m3] -> r [mol/(m3 bed s)]. Basis: per unit bed volume (catalyst lumped
into k_m; no catalyst-mass coordinate).
SIGN CONVENTION (used everywhere): r >= 0 is the rate of CONSUMPTION of o-xylene; the production rate
of species i is nu_i r with nu = (o-xylene -1, O2 -nu_O2, PA +1, H2O +3).
PRESSURE DEPENDENCE created by this law: at fixed T and X
    r = k(T) y_x0 (1 - X) P^2 / (P_ref R T)   ->   r ~ P^2
one power from c_x ~ P (first order in o-xylene) and one from P/P_ref (first order in the O2 partial
pressure with y_O2 fixed). At P = P_ref the law is identical to the RE_CODE law r = k c_x.
"""
import numpy as np

from . import composition


def k_arrhenius(T, p):
    return p.k_m * np.exp(-p.E_over_R * (1.0 / T - 1.0 / p.T_m))


def rate(X, T, P, p):
    return k_arrhenius(T, p) * (P / p.P_ref) * composition.c_x(X, T, P, p)


def rate_and_partials(X, T, P, p):
    """r and its exact partial derivatives (dr/dX, dr/dT, dr/dP)."""
    k = k_arrhenius(T, p)
    fP = P / p.P_ref
    r = k * fP * composition.c_x(X, T, P, p)
    dr_dX = -k * fP * p.y_x0 * P / (p.R_gas * T)
    dr_dT = r * (p.E_over_R / T ** 2 - 1.0 / T)
    dr_dP = 2.0 * r / P
    return r, dr_dX, dr_dT, dr_dP
