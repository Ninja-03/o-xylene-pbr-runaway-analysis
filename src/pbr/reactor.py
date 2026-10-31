"""(8) Reactor ODE system and guarded integration (dimensional Phase 2 model).

State y = [X, T, P, Q] (conversion [-], temperature [K], pressure [Pa], cooling duty [W]):
    dX/dz = A_c r / N_x0                                          (material.py)
    dT/dz = [(-dH) A_c r - 2 pi R_tube U (T - T_c)] / (Qrho c_p)  (energy.py)
    dP/dz = -(Pi P_ref^2/L)(T/T_c)/P   in modes 'outlet'/'inlet', otherwise 0   (pressure.py)
    dQ/dz = 2 pi R_tube U (T - T_c)    auxiliary quadrature (does not feed back)
r from kinetics.py (r >= 0 = o-xylene consumption). The exact analytic Jacobian is supplied to the
stiff solvers. integrate() never raises: every outcome carries a RunStatus (guards.py).
"""
import dataclasses
import time

import numpy as np
from scipy.integrate import solve_ivp

from . import energy, guards, kinetics, material, numerics, pressure
from .guards import GuardError, RunStatus

IX, IT, IP, IQ = 0, 1, 2, 3
GUARD_NAMES = ("pressure P <= P_min", "conversion X > 1", "conversion X < 0", "temperature T <= 0")
_FAILURE_ORDER = (RunStatus.SOLVER_FAILURE, RunStatus.INVALID_PHYSICAL_STATE)


@dataclasses.dataclass
class Solution:
    status: RunStatus
    detail: str
    mode: str
    Pi: float
    P0: float
    config: str
    solver: str
    rtol: float
    atol: np.ndarray
    z: np.ndarray
    y: np.ndarray
    dense: object
    t_steps: np.ndarray
    y_steps: np.ndarray
    hs_z: np.ndarray
    hs_y: np.ndarray
    guard_event: object
    nfev: int
    njev: int
    nlu: int
    wall_s: float
    params_sha256: str
    overrides: tuple

    @property
    def X(self):
        return self.y[IX]

    @property
    def T(self):
        return self.y[IT]

    @property
    def P(self):
        return self.y[IP]

    @property
    def Q(self):
        return self.y[IQ]

    @property
    def ok(self):
        return self.status == RunStatus.PASS


def rhs(z, y, p, Pi, pd_on):
    X, T, P = y[IX], y[IT], y[IP]
    r = kinetics.rate(X, T, P, p)
    return np.array([material.dXdz(r, p), energy.dTdz(r, T, p),
                     pressure.dPdz(T, P, p, Pi) if pd_on else 0.0, energy.q_cool(T, p)])


def jac(z, y, p, Pi, pd_on):
    X, T, P = y[IX], y[IT], y[IP]
    _, rX, rT, rP = kinetics.rate_and_partials(X, T, P, p)
    a = material.coeff(p)
    h = (-p.dH_rxn) * p.A_c / energy.heat_capacity_flow(p)
    b = p.perimeter * p.U / energy.heat_capacity_flow(p)
    J = np.zeros((4, 4))
    J[IX, IX], J[IX, IT], J[IX, IP] = a * rX, a * rT, a * rP
    J[IT, IX], J[IT, IT], J[IT, IP] = h * rX, h * rT - b, h * rP
    if pd_on:
        J[IP, IT], J[IP, IP] = pressure.dPdz_partials(T, P, p, Pi)
    J[IQ, IT] = p.perimeter * p.U
    return J


def inlet_pressure(p, mode, P0):
    return p.P_ref if mode in ("none", "inlet") else float(P0)


def _meta(p, mode, Pi, P0, config, solver, rtol, atol):
    return dict(mode=mode, Pi=float(Pi), P0=float("nan") if P0 is None else float(P0), config=config,
                solver=solver, rtol=rtol, atol=atol, params_sha256=p.sha256, overrides=p.overrides)


def _failed(status, detail, meta, t0, stats=None, guard_event=None):
    stats = stats or dict(nfev=0, njev=0, nlu=0)
    return Solution(status=status, detail=detail, z=np.empty(0), y=np.empty((4, 0)), dense=None,
                    t_steps=np.empty(0), y_steps=np.empty((4, 0)), hs_z=np.empty(0), hs_y=np.empty((0, 4)),
                    guard_event=guard_event, wall_s=time.perf_counter() - t0, **stats, **meta)


def integrate(p, mode="none", Pi=0.0, P0=None, config="physical", solver=None, rtol=None,
              atol_scaled=None, n_grid=numerics.N_GRID, check_ranges=True, X0=0.0, fault_after_z=None):
    """Integrate the dimensional model on [0, L]. fault_after_z injects NaN (V8 guard test only)."""
    t0 = time.perf_counter()
    solver = solver or numerics.PRIMARY
    rtol, a_s = numerics.tolerances(p, rtol, atol_scaled)
    atol = a_s * numerics.state_scales(p)
    meta = _meta(p, mode, Pi, P0, config, solver, rtol, atol)
    try:
        guards.check_inputs(p, mode, Pi, P0, config, check_ranges, X0)
    except GuardError as e:
        return _failed(e.status, e.detail, meta, t0)
    pd_on = mode in pressure.PD_MODES
    meta["P0"] = inlet_pressure(p, mode, P0)
    y0 = np.array([X0, p.T0, meta["P0"], 0.0])
    x_tol = numerics.x_bound_tol(rtol)

    def f(z, y):
        out = np.full(4, np.nan) if (fault_after_z is not None and z >= fault_after_z) else rhs(z, y, p, Pi, pd_on)
        if not np.all(np.isfinite(out)):
            raise FloatingPointError(f"non-finite right-hand side at z = {z:.6g} m")
        return out

    def ev_p(z, y):
        return y[IP] - p.P_min

    def ev_xhi(z, y):
        return y[IX] - (1.0 + x_tol)

    def ev_xlo(z, y):
        return y[IX] + x_tol

    def ev_t(z, y):
        return y[IT]

    def ev_hs(z, y):
        return f(z, y)[IT]

    for ev, d in ((ev_p, -1), (ev_xhi, 1), (ev_xlo, -1), (ev_t, -1)):
        ev.terminal, ev.direction = True, d
    ev_hs.terminal, ev_hs.direction = False, -1
    try:
        sol = solve_ivp(f, (0.0, p.L), y0, method=solver, rtol=rtol, atol=atol,
                        jac=lambda z, y: jac(z, y, p, Pi, pd_on), dense_output=True,
                        events=[ev_p, ev_xhi, ev_xlo, ev_t, ev_hs])
    except (FloatingPointError, ValueError, ZeroDivisionError, np.linalg.LinAlgError) as e:
        return _failed(RunStatus.SOLVER_FAILURE, f"{type(e).__name__}: {e}", meta, t0)
    stats = dict(nfev=int(sol.nfev), njev=int(sol.njev), nlu=int(sol.nlu))
    if sol.status == -1:
        return _failed(RunStatus.SOLVER_FAILURE, f"solver failure: {sol.message}", meta, t0, stats)
    guard_event = None
    for i, name in enumerate(GUARD_NAMES):
        if sol.t_events[i].size:
            guard_event = {"guard": name, "z": float(sol.t_events[i][0])}
            break
    zg = np.linspace(0.0, float(sol.t[-1]), n_grid)
    try:
        Y = sol.sol(zg)
    except Exception as e:  # noqa: BLE001 - any interpolation error is a solver failure
        return _failed(RunStatus.SOLVER_FAILURE, f"interpolation failure: {e}", meta, t0, stats, guard_event)
    Y[:, 0], Y[:, -1] = sol.y[:, 0], sol.y[:, -1]
    dP = pressure.dPdz(Y[IT], Y[IP], p, Pi) if pd_on else np.zeros(zg.size)
    problems = []
    if guard_event is not None:
        problems.append((RunStatus.INVALID_PHYSICAL_STATE,
                         f"guard '{guard_event['guard']}' at z = {guard_event['z']:.9g} m"))
    problems += guards.check_solution(Y[IX], Y[IT], Y[IP], dP, p, pd_on, x_tol)
    if not problems and sol.t[-1] < p.L:
        problems.append((RunStatus.SOLVER_FAILURE, f"integration ended at z = {sol.t[-1]} < L"))
    if problems:
        status = min((s for s, _ in problems), key=_FAILURE_ORDER.index)
        detail = "; ".join(d for _, d in problems)
    else:
        status, detail = RunStatus.PASS, "ok"
    hs_y = np.asarray(sol.y_events[4]) if len(sol.y_events[4]) else np.empty((0, 4))
    return Solution(status=status, detail=detail, z=zg, y=Y, dense=sol.sol, t_steps=sol.t, y_steps=sol.y,
                    hs_z=sol.t_events[4], hs_y=hs_y, guard_event=guard_event,
                    wall_s=time.perf_counter() - t0, **stats, **meta)


def evaluate(sol, z):
    """Dense-output evaluation with an explicit range guard (interpolation failure -> SOLVER_FAILURE)."""
    z = np.atleast_1d(np.asarray(z, float))
    z_end = float(sol.z[-1]) if sol.z.size else 0.0
    if sol.dense is None or z.min() < 0.0 or z.max() > z_end:
        raise GuardError(RunStatus.SOLVER_FAILURE, f"interpolation outside the integrated range [0, {z_end}] m")
    out = sol.dense(z)
    if not np.all(np.isfinite(out)):
        raise GuardError(RunStatus.SOLVER_FAILURE, "non-finite interpolated state")
    return out


def step_quadrature(sol, fun, n_gauss=5):
    """int_0^{z_end} fun(z, y) dz with Gauss-Legendre on every solver step of the dense output
    (independent of the auxiliary state Q)."""
    xg, wg = np.polynomial.legendre.leggauss(n_gauss)
    total = 0.0
    for a, b in zip(sol.t_steps[:-1], sol.t_steps[1:]):
        zz = 0.5 * (b - a) * xg + 0.5 * (a + b)
        total += 0.5 * (b - a) * float(np.dot(wg, fun(zz, sol.dense(zz))))
    return total
