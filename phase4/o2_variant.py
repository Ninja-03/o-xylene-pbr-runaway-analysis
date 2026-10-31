"""Phase 4 MODEL-FORM CHECK: oxygen-depletion variant. This is NOT the base model.

Base (Phase 2/3):  r = k(T) (P/P_ref) c_x                      (O2 mole fraction held at y_O2_0)
Variant:           r = k(T) [y_O2(X)/y_O2_0] (P/P_ref) c_x,   y_O2(X) = y_O2_0 - nu_O2 y_x0 X   (Delta n = 0)
i.e. the same first-order O2 partial-pressure dependence (Phase 2 kinetic_law_note: p_O2 = y_O2 P) with
stoichiometric O2 consumption, lambda = nu_O2 y_x0 / y_O2_0. The frozen Phase 3 files are NOT edited: inside
`with o2_depletion():` the two kinetics functions are replaced in memory, so reactor / shooting / sensitivity /
hotspot run the variant unchanged; the originals are restored on exit.
"""
import contextlib
import math

import numpy as np

from phase4 import common  # noqa: F401  (sets sys.path)
from pbr import kinetics, reactor


def lam(p):
    return p.nu_O2 * p.y_x0 / p.y_O2_0


@contextlib.contextmanager
def o2_depletion():
    base_rate, base_part = kinetics.rate, kinetics.rate_and_partials

    def rate(X, T, P, p):
        return base_rate(X, T, P, p) * (1.0 - lam(p) * X)

    def rate_and_partials(X, T, P, p):
        r, rX, rT, rP = base_part(X, T, P, p)
        f = 1.0 - lam(p) * X
        return r * f, rX * f - r * lam(p), rT * f, rP * f

    kinetics.rate, kinetics.rate_and_partials = rate, rate_and_partials
    try:
        yield
    finally:
        kinetics.rate, kinetics.rate_and_partials = base_rate, base_part


def verify(p):
    """(i) isothermal analytic solution of dX/dz = kappa (1 - X)(1 - lam X):
           X = (1 - E)/(1 - lam E), E = exp[-(1 - lam) kappa z];
       (ii) analytic partial derivatives vs central differences; (iii) identical to base at X = 0;
       (iv) frozen functions restored afterwards."""
    q = p.with_(dH_rxn=0.0, U=0.0)
    k = q.k_m * math.exp(-q.E_over_R * (1.0 / q.T0 - 1.0 / q.T_m))
    kappa = q.A_c * k * q.P_ref / (q.N_f * q.R_gas * q.T0)
    l = lam(q)
    X, T, P = 0.4, 700.0, 1.2 * p.P_ref
    base_r0 = kinetics.rate(0.0, T, P, p)
    with o2_depletion():
        s = reactor.integrate(q)
        E = np.exp(-(1.0 - l) * kappa * s.z)
        Xa = (1.0 - E) / (1.0 - l * E)
        err = float(np.max(np.abs(s.X - Xa)))
        _, rX, rT, rP = kinetics.rate_and_partials(X, T, P, p)
        cd = lambda f, x, h: (f(x + h) - f(x - h)) / (2.0 * h)
        dev = max(abs(cd(lambda x: kinetics.rate(x, T, P, p), X, 1e-6) - rX) / abs(rX),
                  abs(cd(lambda t: kinetics.rate(X, t, P, p), T, 1e-3) - rT) / abs(rT),
                  abs(cd(lambda pp: kinetics.rate(X, T, pp, p), P, 1e-1) - rP) / abs(rP))
        same0 = abs(kinetics.rate(0.0, T, P, p) - base_r0) <= 1e-14 * abs(base_r0)
    restored = kinetics.rate.__module__ == "pbr.kinetics" and kinetics.rate_and_partials.__module__ == "pbr.kinetics"
    ok = s.ok and err <= 1e-6 and dev <= 1e-6 and same0 and restored
    return {"lambda": l, "status": s.status.value, "isothermal_X_L_numerical": float(s.X[-1]),
            "isothermal_X_L_analytic": float(Xa[-1]), "isothermal_max_abs_err": err, "partials_max_rel_dev": dev,
            "identical_to_base_at_X0": bool(same0), "frozen_functions_restored": bool(restored), "pass": bool(ok)}
