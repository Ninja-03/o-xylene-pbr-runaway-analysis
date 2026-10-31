"""(10) Hot-spot detection (Phase 2 criteria.hot_spot): T_hs = max over 0 <= z <= L of T(z), z_hs its
location, dT_hs = T_hs - T_c. Interior maxima come from the integrator's event root finding on
dT/dz = 0 (direction -1, i.e. a maximum); inlet and outlet boundary values are candidates as well; the
dense-profile maximum is reported as a cross-check. A calculated hot spot is a MODEL OUTPUT: the
kinetic validity range is undocumented (Phase 2 C12), so it is not an engineering allowable
temperature and no temperature limit is applied here.
"""
import numpy as np


def detect(sol, p):
    cands = [(0.0, float(sol.T[0]), "inlet"), (float(sol.z[-1]), float(sol.T[-1]), "outlet")]
    for zz, yy in zip(sol.hs_z, sol.hs_y):
        cands.append((float(zz), float(yy[1]), "interior"))
    z_hs, T_hs, kind = max(cands, key=lambda c: c[1])
    i = int(np.argmax(sol.T))
    return {"z_hs": z_hs, "T_hs": T_hs, "dT_hs": T_hs - p.Tc, "kind": kind,
            "n_interior_maxima": int(len(sol.hs_z)), "grid_T_max": float(sol.T[i]),
            "grid_z_max": float(sol.z[i]), "event_minus_grid_K": T_hs - float(sol.T[i])}


def metrics(sol, p):
    """Scalar outputs used by verification and the sanity run."""
    hs = detect(sol, p)
    return {"X_L": float(sol.X[-1]), "T_L": float(sol.T[-1]), "T_max": hs["T_hs"], "z_hs": hs["z_hs"],
            "hs_kind": hs["kind"], "P0": float(sol.P[0]), "P_L": float(sol.P[-1]),
            "dP": float(sol.P[0] - sol.P[-1]), "Q_cool": float(sol.Q[-1])}
