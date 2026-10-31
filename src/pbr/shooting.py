"""(9) Outlet-pressure shooting: find the inlet pressure P0 with P(L) = P_ref (Phase 2 primary spec).

Rigorous bracket. The energy balance bounds the temperature by
    T_low = min(T0, T_c) <= T(z) <= T_up = max(T0, T_c) + dT_ad,
and the exact identity P(L)^2 = P0^2 - 2 Pi P_ref^2 int_0^1 (T/T_c) dzeta then gives
    P0_lo = P_ref sqrt(1 + 2 Pi T_low/T_c)  ->  P(L) <= P_ref,
    P0_hi = P_ref sqrt(1 + 2 Pi T_up /T_c)  ->  P(L) >= P_ref.
Residual F(P0) = P(L)/P_ref - 1. A trial stopped by the pressure guard at z_g < L is an
INVALID_PHYSICAL_STATE (never 'runaway'); for bracketing only, F is continued there by
F = (P_min - P_ref)/P_ref - (1 - z_g/L) (negative, continuous at z_g = L).
The bracket is widened by BRACKET_MARGIN_REL * rtol (numerics.py) to absorb integration error.
Before Brent's method, a scan over the bracket must show exactly one sign change and a strictly
increasing F (uniqueness of the shooting solution); otherwise ROOT_FAILURE.
"""
import dataclasses
import math

import numpy as np
from scipy.optimize import brentq

from . import numerics, reactor
from .guards import RunStatus


@dataclasses.dataclass
class ShootResult:
    status: RunStatus
    detail: str
    P0: float
    residual: float
    solution: object
    bracket: tuple
    scan_P0: np.ndarray
    scan_F: np.ndarray
    scan_status: list
    monotonic: bool
    sign_changes: int
    iterations: int
    function_calls: int


def analytic_bracket(p, Pi, rtol=None):
    """Rigorous bracket, widened by BRACKET_MARGIN_REL * rtol so that integration error at the ends cannot
    flip the residual sign (the exact bracket collapses to one point when T is known exactly, e.g. dH = U = 0)."""
    rt, _ = numerics.tolerances(p, rtol)
    m = numerics.BRACKET_MARGIN_REL * rt
    T_low = min(p.T0, p.Tc)
    T_up = max(p.T0, p.Tc) + p.dT_ad
    return (p.P_ref * math.sqrt(1.0 + 2.0 * Pi * T_low / p.Tc) * (1.0 - m),
            p.P_ref * math.sqrt(1.0 + 2.0 * Pi * T_up / p.Tc) * (1.0 + m))


def residual(P0, p, Pi, **kw):
    sol = reactor.integrate(p, mode="outlet", Pi=Pi, P0=P0, **kw)
    if sol.status == RunStatus.PASS:
        return sol.P[-1] / p.P_ref - 1.0, sol
    ge = sol.guard_event
    if sol.status == RunStatus.INVALID_PHYSICAL_STATE and ge and ge["guard"].startswith("pressure"):
        return (p.P_min - p.P_ref) / p.P_ref - (1.0 - ge["z"] / p.L), sol
    return float("nan"), sol


def solve_outlet(p, Pi, n_scan=numerics.SHOOT_SCAN, bracket=None, xtol_rel=None, res_tol=None, **kw):
    """Find P0 with P(L) = P_ref. kw go to reactor.integrate (solver, rtol, atol_scaled, config, ...)."""
    rtol, _ = numerics.tolerances(p, kw.get("rtol"), kw.get("atol_scaled"))
    xtol_rel = numerics.SHOOT_XTOL_REL * rtol if xtol_rel is None else xtol_rel
    res_tol = numerics.SHOOT_RES_REL * rtol if res_tol is None else res_tol
    if Pi == 0.0:
        sol = reactor.integrate(p, mode="outlet", Pi=0.0, P0=p.P_ref, **kw)
        return ShootResult(sol.status, "Pi = 0: P0 = P_ref (no pressure drop)", p.P_ref, 0.0, sol,
                           (p.P_ref, p.P_ref), np.array([p.P_ref]), np.array([0.0]), [sol.status], True, 0, 0, 1)
    lo, hi = bracket if bracket is not None else analytic_bracket(p, Pi, rtol)
    P0s = np.linspace(lo, hi, max(int(n_scan), 2))
    Fs, sts = [], []
    for P0 in P0s:
        F, s = residual(P0, p, Pi, **kw)
        Fs.append(F)
        sts.append(s.status)
    Fs = np.array(Fs)
    finite = bool(np.all(np.isfinite(Fs)))
    mono = finite and bool(np.all(np.diff(Fs) > 0.0))
    sgn = np.sign(Fs) if finite else np.zeros(1)
    nsc = int(np.sum(sgn[:-1] * sgn[1:] < 0) + np.sum(sgn == 0)) if finite else -1

    def result(status, detail, P0=float("nan"), res=float("nan"), sol=None, it=0, calls=0):
        return ShootResult(status, detail, P0, res, sol, (lo, hi), P0s, Fs, sts, mono, nsc, it, len(P0s) + calls)

    if not finite:
        return result(RunStatus.ROOT_FAILURE, "a trial integration failed for a reason other than the pressure guard")
    if Fs[0] > 0.0 or Fs[-1] < 0.0:
        return result(RunStatus.ROOT_FAILURE,
                      f"invalid bracket: F(lo) = {Fs[0]:.3e}, F(hi) = {Fs[-1]:.3e} (need F(lo) <= 0 <= F(hi))")
    if not mono or nsc != 1:
        return result(RunStatus.ROOT_FAILURE, f"P(L) not strictly increasing in P0 over the bracket "
                      f"(monotonic={mono}, sign changes={nsc}): uniqueness not established")
    i = int(np.nonzero((Fs[:-1] <= 0.0) & (Fs[1:] >= 0.0))[0][0])
    calls = [0]

    def fun(x):
        calls[0] += 1
        return residual(x, p, Pi, **kw)[0]

    try:
        if Fs[i] == 0.0 or Fs[i + 1] == 0.0:
            root, it, conv = (P0s[i] if Fs[i] == 0.0 else P0s[i + 1]), 0, True
        else:
            root, info = brentq(fun, P0s[i], P0s[i + 1], xtol=xtol_rel * p.P_ref,
                                maxiter=numerics.BRENT_MAXITER, full_output=True, disp=False)
            it, conv = int(info.iterations), bool(info.converged)
    except (ValueError, RuntimeError) as e:
        return result(RunStatus.ROOT_FAILURE, f"Brent failure: {e}", calls=calls[0])
    F, sol = residual(root, p, Pi, **kw)
    if not conv or sol.status != RunStatus.PASS or not abs(F) <= res_tol:
        return result(RunStatus.ROOT_FAILURE, f"shooting not accepted: converged={conv}, status="
                      f"{sol.status.value}, residual={F:.3e} (tol {res_tol:.1e})", float(root), F, sol, it, calls[0] + 1)
    return result(RunStatus.PASS, "ok", float(root), float(F), sol, it, calls[0] + 1)
