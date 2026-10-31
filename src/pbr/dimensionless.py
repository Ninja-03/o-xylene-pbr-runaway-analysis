"""(12) Independent dimensionless implementation of the Phase 2 model.

Coded separately from reactor.py / kinetics.py / shooting.py (none of them is imported): the groups
are computed here from the raw parameter-file values, the equations are coded here, and the outlet
shooting uses its own bracket and root finding.
  gamma = (E/R)/T_c        dT_ad = (-dH) y_x0/(M_wf c_p)        B = gamma dT_ad/T_c
  Da = A_c k(T_c) P_ref L/(N_f R T_c),  k(T_c) = k_m exp[(E/R)(1/T_m - 1/T_c)],  N_f = Qrho/M_wf
  St = 2 pi R_tube U L/(Qrho c_p)       theta = gamma (T - T_c)/T_c,  p = P/P_ref,  zeta = z/L
  dX/dzeta     = Da p^2 (1 - X) exp[theta/(1 + theta/gamma)] / (1 + theta/gamma)
  dtheta/dzeta = B dX/dzeta - St theta
  dp/dzeta     = -Pi (1 + theta/gamma)/p      (pressure-drop modes only)
  dI/dzeta     = theta                         (auxiliary: Q_cool = Qrho c_p (T_c/gamma) St I)
Initial conditions: X = 0, theta = theta0 = gamma (T0 - T_c)/T_c, p = p0.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import numerics
from .guards import RunStatus


def groups(p):
    A_c = math.pi * p.R_tube ** 2
    N_f = p.Qrho / p.M_wf
    gamma = p.E_over_R / p.Tc
    dT_ad = (-p.dH_rxn) * p.y_x0 / (p.M_wf * p.cp)
    k_c = p.k_m * math.exp(p.E_over_R * (1.0 / p.T_m - 1.0 / p.Tc))
    return {"gamma": gamma, "dT_ad": dT_ad, "B": gamma * dT_ad / p.Tc,
            "Da": A_c * k_c * p.P_ref * p.L / (N_f * p.R_gas * p.Tc),
            "St": 2.0 * math.pi * p.R_tube * p.U * p.L / (p.Qrho * p.cp),
            "theta0": gamma * (p.T0 - p.Tc) / p.Tc}


def rhs(zeta, y, g, Pi, pd_on):
    X, th, pp = y[0], y[1], y[2]
    s = 1.0 + th / g["gamma"]
    dX = g["Da"] * pp * pp * (1.0 - X) * np.exp(th / s) / s
    return np.array([dX, g["B"] * dX - g["St"] * th, -Pi * s / pp if pd_on else 0.0, th])


def integrate(p, mode="none", Pi=0.0, p0=None, solver=None, rtol=None, atol_scaled=None, n_grid=numerics.N_GRID):
    g = groups(p)
    solver = solver or numerics.PRIMARY
    rt, a_s = numerics.tolerances(p, rtol, atol_scaled)
    pd_on = mode in ("outlet", "inlet")
    p_in = 1.0 if mode in ("none", "inlet") else float(p0)
    th_scale = max(1.0, g["B"], abs(g["theta0"]))
    atol = a_s * np.array([1.0, th_scale, 1.0, th_scale])
    x_tol = numerics.x_bound_tol(rt)

    def fun(zeta, y):
        return rhs(zeta, y, g, Pi, pd_on)

    def ev_p(zeta, y):
        return y[2] - p.P_min_frac

    def ev_xhi(zeta, y):
        return y[0] - (1.0 + x_tol)

    def ev_xlo(zeta, y):
        return y[0] + x_tol

    def ev_hs(zeta, y):
        return fun(zeta, y)[1]

    for ev, d, term in ((ev_p, -1, True), (ev_xhi, 1, True), (ev_xlo, -1, True), (ev_hs, -1, False)):
        ev.terminal, ev.direction = term, d
    out = {"groups": g, "mode": mode, "Pi": Pi, "solver": solver, "rtol": rt}
    try:
        sol = solve_ivp(fun, (0.0, 1.0), np.array([0.0, g["theta0"], p_in, 0.0]), method=solver, rtol=rt,
                        atol=atol, dense_output=True, events=[ev_p, ev_xhi, ev_xlo, ev_hs])
    except (FloatingPointError, ValueError, ZeroDivisionError) as e:
        out.update(status=RunStatus.SOLVER_FAILURE, detail=str(e))
        return out
    if sol.status == -1:
        out.update(status=RunStatus.SOLVER_FAILURE, detail=sol.message)
        return out
    for i, name in enumerate(("pressure p <= P_min/P_ref", "X > 1", "X < 0")):
        if sol.t_events[i].size:
            out.update(status=RunStatus.INVALID_PHYSICAL_STATE, detail=name, zeta_guard=float(sol.t_events[i][0]))
            return out
    zeta = np.linspace(0.0, 1.0, n_grid)
    Y = sol.sol(zeta)
    Y[:, 0], Y[:, -1] = sol.y[:, 0], sol.y[:, -1]
    cands = [(0.0, float(Y[1, 0])), (1.0, float(Y[1, -1]))]
    cands += [(float(z), float(y[1])) for z, y in zip(sol.t_events[3], sol.y_events[3])]
    z_hs, th_hs = max(cands, key=lambda c: c[1])
    out.update(status=RunStatus.PASS, detail="ok", zeta=zeta, X=Y[0], theta=Y[1], p=Y[2], I=Y[3],
               z=zeta * p.L, T=p.Tc * (1.0 + Y[1] / g["gamma"]), P=Y[2] * p.P_ref,
               Q_cool=p.Qrho * p.cp * (p.Tc / g["gamma"]) * g["St"] * float(sol.y[3, -1]),
               zeta_hs=z_hs, z_hs=z_hs * p.L, theta_hs=th_hs, T_hs=p.Tc * (1.0 + th_hs / g["gamma"]),
               dense=sol.sol, nfev=int(sol.nfev))
    return out


def shoot_outlet(p, Pi, n_scan=numerics.SHOOT_SCAN, solver=None, rtol=None, atol_scaled=None, n_grid=numerics.N_GRID):
    """Independent outlet shooting on p0 = P0/P_ref (bracket from the theta bounds of the energy balance)."""
    g = groups(p)
    rt, _ = numerics.tolerances(p, rtol, atol_scaled)
    kw = dict(solver=solver, rtol=rtol, atol_scaled=atol_scaled, n_grid=n_grid)
    if Pi == 0.0:
        r = integrate(p, "outlet", 0.0, 1.0, **kw)
        return {"status": r["status"], "detail": "Pi = 0", "p0": 1.0, "residual": 0.0, "result": r}
    s_lo = 1.0 + min(g["theta0"], 0.0) / g["gamma"]
    s_hi = 1.0 + (max(g["theta0"], 0.0) + g["B"]) / g["gamma"]
    m = numerics.BRACKET_MARGIN_REL * rt          # same integration-error margin as the dimensional shooting
    lo, hi = math.sqrt(1.0 + 2.0 * Pi * s_lo) * (1.0 - m), math.sqrt(1.0 + 2.0 * Pi * s_hi) * (1.0 + m)

    def F(p0):
        r = integrate(p, "outlet", Pi, p0, **kw)
        if r["status"] == RunStatus.PASS:
            return float(r["p"][-1]) - 1.0, r
        if r["status"] == RunStatus.INVALID_PHYSICAL_STATE and r["detail"].startswith("pressure"):
            return (p.P_min_frac - 1.0) - (1.0 - r["zeta_guard"]), r
        return float("nan"), r

    ps = np.linspace(lo, hi, n_scan)
    Fs = np.array([F(x)[0] for x in ps])
    out = {"bracket": (lo, hi), "scan_p0": ps, "scan_F": Fs}
    if not np.all(np.isfinite(Fs)) or Fs[0] > 0.0 or Fs[-1] < 0.0 or not np.all(np.diff(Fs) > 0.0):
        out.update(status=RunStatus.ROOT_FAILURE, detail="invalid bracket or non-monotonic residual",
                   p0=float("nan"), residual=float("nan"), result=None)
        return out
    i = int(np.nonzero((Fs[:-1] <= 0.0) & (Fs[1:] >= 0.0))[0][0])
    try:
        root = ps[i] if Fs[i] == 0.0 else brentq(lambda x: F(x)[0], ps[i], ps[i + 1],
                                                  xtol=numerics.SHOOT_XTOL_REL * rt, maxiter=numerics.BRENT_MAXITER)
    except (ValueError, RuntimeError) as e:
        out.update(status=RunStatus.ROOT_FAILURE, detail=str(e), p0=float("nan"), residual=float("nan"), result=None)
        return out
    res, r = F(root)
    good = r["status"] == RunStatus.PASS and abs(res) <= numerics.SHOOT_RES_REL * rt
    out.update(status=RunStatus.PASS if good else RunStatus.ROOT_FAILURE, detail="ok" if good else f"residual {res:.3e}",
               p0=float(root), residual=float(res), result=r)
    return out
