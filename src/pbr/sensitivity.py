"""(11) Forward sensitivity equations and the SENSITIVITY-BASED PARAMETRIC RUNAWAY BOUNDARY (Phase 2, ZB2020).

S(z) = d[X, T, P](z)/dT0 with T_c held fixed:  dS/dz = J(y) S,  J = exact 3x3 Jacobian of
[dX/dz, dT/dz, dP/dz] (reactor.jac),  S(0) = (0, 1, dP0/dT0).
  modes 'none', 'inlet', 'constant' (P0 fixed):  dP0/dT0 = 0
  mode 'outlet' (P(L) = P_ref held fixed):        dP0/dT0 = -S^T0_P(L) / S^P0_P(L),
  where S^P0 solves the same equations with S(0) = (0, 0, 1).
Phase 2 criteria.runaway_boundary, with S_T the temperature component of S:
  if dS_T/dz(0) >= 0: S_int = max of S_T over (0, L]        (S rises from the inlet: runaway side)
  else:               S_int = max of S_T over [z_min1, L]   (z_min1 = first interior minimum of S_T;
                                                              S_int = S_T(L) if S_T never turns up)
  g = S_int - 1:  g < 0 safe side (inlet disturbances attenuated), g > 0 runaway side (amplified).
Name: "sensitivity-based parametric runaway boundary". It is NOT a dynamic-stability boundary, a
saddle-node point, a bifurcation point, or the exact Rawlings & Ekerdt runaway temperature (the documented
'about 631 K' statement is not reproduced by it; see verification V9). The Morbidelli-Varma generalized
criterion is NOT implemented (primary source not accessed).
U_crit (criteria.Ucrit): the U with g = 0 at fixed T0 (bracketed Brent); T0_crit: the T0 with g = 0
at fixed U (V12 consistency check).
"""
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import guards, numerics, pressure, reactor, shooting
from .guards import GuardError, RunStatus

T0_DIR, P0_DIR = "T0", "P0"


def _scales(p, kind):
    if kind == T0_DIR:
        return np.array([1.0 / p.Tc, 1.0, p.P_ref / p.Tc])
    return np.array([1.0 / p.P_ref, p.Tc / p.P_ref, 1.0])


def _aug_fun(p, Pi, pd_on, n_vec):
    def F(z, Y):
        y = Y[:4]
        J = reactor.jac(z, y, p, Pi, pd_on)[:3, :3]
        parts = [reactor.rhs(z, y, p, Pi, pd_on)]
        for k in range(n_vec):
            parts.append(J @ Y[4 + 3 * k: 7 + 3 * k])
        out = np.concatenate(parts)
        if not np.all(np.isfinite(out)):
            raise FloatingPointError(f"non-finite augmented right-hand side at z = {z:.6g} m")
        return out
    return F


def integrate_sensitivity(p, mode="none", Pi=0.0, P0=None, vectors=((0.0, 1.0, 0.0),), kinds=(T0_DIR,),
                          config="physical", solver=None, rtol=None, atol_scaled=None, check_ranges=True,
                          n_grid=numerics.N_GRID, extrema=True):
    t0 = time.perf_counter()
    solver = solver or numerics.PRIMARY
    rtol, a_s = numerics.tolerances(p, rtol, atol_scaled)
    out = {"mode": mode, "Pi": Pi, "solver": solver, "rtol": rtol, "atol_scaled": a_s}
    try:
        guards.check_inputs(p, mode, Pi, P0, config, check_ranges)
    except GuardError as e:
        out.update(status=e.status, detail=e.detail)
        return out
    pd_on = mode in pressure.PD_MODES
    P_in = reactor.inlet_pressure(p, mode, P0)
    n_vec = len(vectors)
    Y0 = np.concatenate([[0.0, p.T0, P_in, 0.0]] + [np.asarray(v, float) for v in vectors])
    atol = a_s * np.concatenate([numerics.state_scales(p)] + [_scales(p, k) for k in kinds])
    F = _aug_fun(p, Pi, pd_on, n_vec)
    x_tol = numerics.x_bound_tol(rtol)

    def ev_p(z, Y):
        return Y[2] - p.P_min

    def ev_xhi(z, Y):
        return Y[0] - (1.0 + x_tol)

    def ev_xlo(z, Y):
        return Y[0] + x_tol

    for ev, d in ((ev_p, -1), (ev_xhi, 1), (ev_xlo, -1)):
        ev.terminal, ev.direction = True, d
    events = [ev_p, ev_xhi, ev_xlo]
    if extrema:
        def ev_min(z, Y):
            return F(z, Y)[5]

        def ev_max(z, Y):
            return F(z, Y)[5]
        ev_min.terminal, ev_min.direction = False, 1
        ev_max.terminal, ev_max.direction = False, -1
        events += [ev_min, ev_max]
    try:
        sol = solve_ivp(F, (0.0, p.L), Y0, method=solver, rtol=rtol, atol=atol, dense_output=True, events=events)
    except (FloatingPointError, ValueError, ZeroDivisionError, np.linalg.LinAlgError) as e:
        out.update(status=RunStatus.SOLVER_FAILURE, detail=f"{type(e).__name__}: {e}")
        return out
    if sol.status == -1:
        out.update(status=RunStatus.SOLVER_FAILURE, detail=sol.message)
        return out
    for i, name in enumerate(reactor.GUARD_NAMES[:3]):
        if sol.t_events[i].size:
            out.update(status=RunStatus.INVALID_PHYSICAL_STATE,
                       detail=f"guard '{name}' at z = {sol.t_events[i][0]:.9g} m")
            return out
    zg = np.linspace(0.0, p.L, n_grid)
    Y = sol.sol(zg)
    Y[:, 0], Y[:, -1] = sol.y[:, 0], sol.y[:, -1]
    J0 = reactor.jac(0.0, Y0[:4], p, Pi, pd_on)[:3, :3]
    out.update(status=RunStatus.PASS, detail="ok", z=zg, y=Y[:4],
               S=[Y[4 + 3 * k: 7 + 3 * k] for k in range(n_vec)],
               S_L=[sol.y[4 + 3 * k: 7 + 3 * k, -1] for k in range(n_vec)],
               dS0=float((J0 @ np.asarray(vectors[0], float))[1]), dense=sol.sol, t_steps=sol.t,
               nfev=int(sol.nfev), wall_s=time.perf_counter() - t0, P0=P_in)
    if extrema:
        out.update(min_z=sol.t_events[3], max_z=sol.t_events[4],
                   min_S=np.asarray(sol.y_events[3])[:, 5] if sol.t_events[3].size else np.empty(0),
                   max_S=np.asarray(sol.y_events[4])[:, 5] if sol.t_events[4].size else np.empty(0))
    return out


def interior_max(res):
    """S_int per the Phase 2 criterion (see module docstring)."""
    L = float(res["z"][-1])
    ST_L = float(res["S_L"][0][1])
    if res["dS0"] >= 0.0:
        z_start, rule = 0.0, "dS_T/dz(0) >= 0: max of S_T over (0, L]"
    elif res["min_z"].size:
        z_start, rule = float(res["min_z"][0]), "max of S_T over [first interior minimum, L]"
    else:
        z_start, rule = L, "S_T decreases monotonically from the inlet: S_int = S_T(L)"
    cands = [(float(z), float(s)) for z, s in zip(res["max_z"], res["max_S"]) if z > z_start] + [(L, ST_L)]
    z_star, S_int = max(cands, key=lambda c: c[1])
    return S_int, z_star, z_start, rule


def criterion(p, mode="none", Pi=0.0, P0=None, config="physical", solver=None, rtol=None, atol_scaled=None,
              check_ranges=True):
    """g = S_int - 1 at one operating point (Phase 2 runaway-boundary criterion)."""
    kw = dict(config=config, solver=solver, rtol=rtol, atol_scaled=atol_scaled, check_ranges=check_ranges)
    out = {"mode": mode, "Pi": Pi, "T0": p.T0, "Tc": p.Tc, "U": p.U, "g": float("nan")}
    dP0dT0 = 0.0
    if mode == "outlet":
        sr = shooting.solve_outlet(p, Pi, **kw)
        if sr.status != RunStatus.PASS:
            out.update(status=RunStatus.ROOT_FAILURE, detail=f"shooting: {sr.detail}")
            return out
        P0 = sr.P0
        two = integrate_sensitivity(p, mode, Pi, P0, vectors=((0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
                                    kinds=(T0_DIR, P0_DIR), extrema=False, **kw)
        if two["status"] != RunStatus.PASS:
            out.update(status=two["status"], detail=two["detail"])
            return out
        dP0dT0 = -float(two["S_L"][0][2]) / float(two["S_L"][1][2])
    res = integrate_sensitivity(p, mode, Pi, P0, vectors=((0.0, 1.0, dP0dT0),), kinds=(T0_DIR,), extrema=True, **kw)
    if res["status"] != RunStatus.PASS:
        out.update(status=res["status"], detail=res["detail"])
        return out
    S_int, z_star, z_start, rule = interior_max(res)
    out.update(status=RunStatus.PASS, detail="ok", g=S_int - 1.0, S_int=S_int, z_Sint=z_star, z_start=z_start,
               rule=rule, dS0=res["dS0"], dP0dT0=dP0dT0, P0=res["P0"], S_T_L=float(res["S_L"][0][1]),
               S_P_L=float(res["S_L"][0][2]), n_min=int(res["min_z"].size), n_max=int(res["max_z"].size), res=res)
    return out


def _scan_and_brent(evalg, xs, expect, xtol):
    recs = [evalg(x) for x in xs]
    gs = np.array([r["g"] if r["status"] == RunStatus.PASS else np.nan for r in recs])
    info = {"scan_x": np.asarray(xs, float), "scan_g": gs, "scan_status": [r["status"].value for r in recs]}
    if not np.all(np.isfinite(gs)):
        return dict(status=RunStatus.ROOT_FAILURE, detail="criterion evaluation failed at a scan point", **info)
    sgn = np.sign(gs)
    changes = np.nonzero(sgn[:-1] * sgn[1:] < 0)[0]
    d = np.diff(gs)
    info["sign_changes"] = int(changes.size)
    info["globally_monotonic"] = bool(np.all(d < 0)) if expect == "decreasing" else bool(np.all(d > 0))
    lo_ok = gs[0] > 0 if expect == "decreasing" else gs[0] < 0
    hi_ok = gs[-1] < 0 if expect == "decreasing" else gs[-1] > 0
    if not (lo_ok and hi_ok):
        return dict(status=RunStatus.ROOT_FAILURE, detail=f"no valid bracket: g = {gs[0]:.4g} at {xs[0]:g}, "
                    f"{gs[-1]:.4g} at {xs[-1]:g} (no boundary in range)", **info)
    if changes.size != 1:
        return dict(status=RunStatus.ROOT_FAILURE, detail=f"{changes.size} sign changes: root not unique", **info)
    i = int(changes[0])
    # Uniqueness at finer resolution: a 9-point sub-scan of the root bracket must also change sign exactly once.
    # (Monotonicity is NOT required away from the root: g can be discontinuous on the g < 0 side where the
    # interior maximum of S_T disappears and S_int switches to S_T(L); that region contains no root.)
    sub_x = np.linspace(xs[i], xs[i + 1], 9)
    sub = [evalg(x) for x in sub_x]
    sub_g = np.array([r["g"] if r["status"] == RunStatus.PASS else np.nan for r in sub])
    s2 = np.sign(sub_g)
    info.update(subscan_x=sub_x, subscan_g=sub_g, subscan_sign_changes=int(np.sum(s2[:-1] * s2[1:] < 0)))
    if not np.all(np.isfinite(sub_g)) or info["subscan_sign_changes"] != 1:
        return dict(status=RunStatus.ROOT_FAILURE, detail="root not unique inside the bracket (sub-scan)", **info)
    j = int(np.nonzero(s2[:-1] * s2[1:] < 0)[0][0])
    lo_b, hi_b = sub_x[j], sub_x[j + 1]

    def g_only(x):
        r = evalg(x)
        if r["status"] != RunStatus.PASS:
            raise RuntimeError(f"criterion failed at {x}: {r['detail']}")
        return r["g"]

    try:
        root, rr = brentq(g_only, lo_b, hi_b, xtol=xtol, maxiter=numerics.BRENT_MAXITER,
                          full_output=True, disp=False)
    except (ValueError, RuntimeError) as e:
        return dict(status=RunStatus.ROOT_FAILURE, detail=f"Brent failure: {e}", **info)
    at = evalg(root)
    # Monotonicity WHERE REQUIRED: strictly monotonic through the root (root -/+ xtol and -/+ 10 xtol).
    # probe both sides of the root inside the full scan range (NOT the sub-bracket: the root may sit on its edge)
    pts = np.unique(np.clip(root + np.array([-10.0, -1.0, 1.0, 10.0]) * xtol, xs[0], xs[-1]))
    if pts.size < 3:
        return dict(status=RunStatus.ROOT_FAILURE, detail="root at the edge of the scan range: cannot probe both sides",
                    root=float(root), **info)
    loc_g = np.array([evalg(x)["g"] for x in pts])
    dl = np.diff(loc_g)
    info.update(local_x=pts, local_g=loc_g,
                local_monotonic_at_root=bool(np.all(dl < 0)) if expect == "decreasing" else bool(np.all(dl > 0)))
    good = bool(rr.converged) and at["status"] == RunStatus.PASS and info["local_monotonic_at_root"]
    info.update(bracket=(float(lo_b), float(hi_b)), coarse_bracket=(float(xs[i]), float(xs[i + 1])),
                iterations=int(rr.iterations),
                function_calls=int(rr.function_calls), converged=bool(rr.converged),
                g_at_root=float(at.get("g", np.nan)),
                at_root={k: at[k] for k in ("S_int", "z_Sint", "z_start", "rule", "dS0", "dP0dT0", "n_min", "n_max")
                         if k in at})
    return dict(status=RunStatus.PASS if good else RunStatus.ROOT_FAILURE,
                detail="ok" if good else "Brent not converged or g not monotonic through the root",
                root=float(root), xtol=xtol, **info)


def find_Ucrit(p, config="physical", mode="none", Pi=0.0, P0=None, U_bracket=None, n_scan=9, xtol=None,
               solver=None, rtol=None, atol_scaled=None):
    lo, hi = U_bracket if U_bracket is not None else p.ranges["U"]
    xtol = float(p.numerics["root_tolerance_U"]) if xtol is None else xtol

    def evalg(U):
        return criterion(p.with_(U=U), mode=mode, Pi=Pi, P0=P0, config=config, solver=solver, rtol=rtol,
                         atol_scaled=atol_scaled)
    r = _scan_and_brent(evalg, np.linspace(lo, hi, n_scan), "decreasing", xtol)
    r["U_crit"] = r.get("root", float("nan"))
    return r


def find_T0crit(p, U, config="physical", mode="none", Pi=0.0, T0_bracket=None, n_scan=9, xtol=None,
                solver=None, rtol=None, atol_scaled=None):
    lo, hi = T0_bracket if T0_bracket is not None else p.ranges["T0"]
    xtol = float(p.numerics["boundary_resolution_T0"]) if xtol is None else xtol

    def point(T0):
        return p.with_(U=U, T0=T0, Tc=T0) if config == "normalized" else p.with_(U=U, T0=T0)

    def evalg(T0):
        return criterion(point(T0), mode=mode, Pi=Pi, config=config, solver=solver, rtol=rtol,
                         atol_scaled=atol_scaled)
    r = _scan_and_brent(evalg, np.linspace(lo, hi, n_scan), "increasing", xtol)
    r["T0_crit"] = r.get("root", float("nan"))
    return r
