"""(13) Numerical and physical guards; machine-readable run statuses.

Before integration: INVALID_INPUT (non-finite or non-physical input) or OUT_OF_MODEL_RANGE (outside the
Phase 2 ranges or the Phase 2 feasible range of a pressure configuration).
During/after integration: INVALID_PHYSICAL_STATE (P <= P_min, X outside [0, 1], T <= 0, dP/dz >= 0 in
a pressure-drop mode), SOLVER_FAILURE (solver error, NaN/Inf, interpolation failure),
ROOT_FAILURE (invalid bracket, non-monotonic residual, non-convergence).
There is deliberately NO 'runaway' status: guard violations and numerical failures are never
interpreted as thermal runaway. Runaway is decided only by the Phase 2 sensitivity criterion.
"""
import enum
import math

import numpy as np

from . import params as _params
from . import pressure


class RunStatus(str, enum.Enum):
    PASS = "PASS"
    INVALID_INPUT = "INVALID_INPUT"
    OUT_OF_MODEL_RANGE = "OUT_OF_MODEL_RANGE"
    INVALID_PHYSICAL_STATE = "INVALID_PHYSICAL_STATE"
    SOLVER_FAILURE = "SOLVER_FAILURE"
    ROOT_FAILURE = "ROOT_FAILURE"


class GuardError(Exception):
    def __init__(self, status, detail):
        super().__init__(f"{status.value}: {detail}")
        self.status = status
        self.detail = detail


CONFIGS = ("physical", "normalized", "verification")
_POSITIVE = ("k_m", "T_m", "E_over_R", "y_x0", "M_wf", "P_ref", "T0", "Tc", "Qrho", "R_tube", "L",
             "cp", "R_gas", "P_min_frac")


def check_inputs(p, mode, Pi, P0=None, config="physical", check_ranges=True, X0=0.0):
    bad_in, bad_rng = RunStatus.INVALID_INPUT, RunStatus.OUT_OF_MODEL_RANGE
    for f in _params.FIELDS:
        v = getattr(p, f)
        if not math.isfinite(v):
            raise GuardError(bad_in, f"{f} = {v} is not finite")
    for f in _POSITIVE:
        if getattr(p, f) <= 0.0:
            raise GuardError(bad_in, f"{f} = {getattr(p, f)} must be > 0")
    if p.U < 0.0:
        raise GuardError(bad_in, f"U = {p.U} must be >= 0")
    if mode not in pressure.MODES:
        raise GuardError(bad_in, f"unknown pressure mode {mode!r}")
    if config not in CONFIGS:
        raise GuardError(bad_in, f"unknown configuration {config!r}")
    if not math.isfinite(Pi):
        raise GuardError(bad_in, f"Pi = {Pi} is not finite")
    if mode == "none" and Pi != 0.0:
        raise GuardError(bad_in, "baseline mode 'none' has no pressure drop (Pi must be 0)")
    if mode in ("outlet", "constant") and (P0 is None or not math.isfinite(P0) or P0 <= 0.0):
        raise GuardError(bad_in, f"mode {mode!r} needs a finite inlet pressure P0 > 0 (got {P0})")
    if not (0.0 <= X0 <= 1.0):
        raise GuardError(RunStatus.INVALID_PHYSICAL_STATE, f"initial conversion X0 = {X0} outside [0, 1]")
    if not check_ranges:
        return
    if config == "verification":
        raise GuardError(bad_in, "configuration 'verification' is only allowed with check_ranges=False")
    if Pi < 0.0 or Pi > p.Pi_max:
        raise GuardError(bad_rng, f"Pi = {Pi} outside the Phase 2 study range [0, {p.Pi_max}]")
    if mode == "inlet" and Pi > p.Pi_inlet_max:
        raise GuardError(bad_rng, f"inlet-specified comparison only for Pi <= {p.Pi_inlet_max} (Phase 2); got {Pi}")
    for f in ("U", "y_x0", "cp", "T0"):
        lo, hi = p.ranges[f]
        v = getattr(p, f)
        if not (lo <= v <= hi):
            raise GuardError(bad_rng, f"{f} = {v} outside the Phase 2 range [{lo}, {hi}]")
    if config == "physical":
        lo, hi = p.ranges["Tc"]
        if not (lo <= p.Tc <= hi):
            raise GuardError(bad_rng, f"Tc = {p.Tc} outside the Phase 2 salt range [{lo}, {hi}]")
    elif config == "normalized" and abs(p.Tc - p.T0) > 1e-12 * p.T0:
        raise GuardError(bad_rng, "normalized configuration requires Tc = T0")


def check_solution(X, T, P, dPdz, p, pd_on, x_tol):
    """Post-integration state checks on the evaluated profile; returns [(status, detail), ...]."""
    if not (np.all(np.isfinite(X)) and np.all(np.isfinite(T)) and np.all(np.isfinite(P))):
        return [(RunStatus.SOLVER_FAILURE, "non-finite state (NaN/Inf) in the solution")]
    out = []
    if X.min() < -x_tol or X.max() > 1.0 + x_tol:
        out.append((RunStatus.INVALID_PHYSICAL_STATE,
                    f"conversion outside [0, 1] (allowance {x_tol:.1e}): min {X.min():.3e}, max {X.max():.12f}"))
    if T.min() <= 0.0:
        out.append((RunStatus.INVALID_PHYSICAL_STATE, f"temperature <= 0 (min {T.min():.6g} K)"))
    if pd_on:
        if P.min() <= p.P_min:
            out.append((RunStatus.INVALID_PHYSICAL_STATE, f"P <= P_min = {p.P_min:.6g} Pa (min {P.min():.6g} Pa)"))
        if np.any(dPdz >= 0.0):
            out.append((RunStatus.INVALID_PHYSICAL_STATE, f"dP/dz >= 0 at {int(np.sum(dPdz >= 0.0))} points"))
        if np.any(np.diff(P) >= 0.0):
            out.append((RunStatus.INVALID_PHYSICAL_STATE, "P(z) not strictly decreasing"))
    elif P.min() <= 0.0:
        out.append((RunStatus.INVALID_PHYSICAL_STATE, "P <= 0"))
    return out
