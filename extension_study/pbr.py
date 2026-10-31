"""
1-D pseudo-homogeneous PBR, o-xylene (R&E dataset), rebuilt from manuscript Eqs (1)-(5), Table 1.
Extensions:
  n      : total-pressure order of the rate, r = k(T) (P/Pref)^(n-1) c_x   (n=2 base, n=1 alternative)
  kappa  : co-current coolant capacity ratio  m*cp/(mc*cpc);  kappa=0 -> constant Tc
Criterion: g = Sint - 1 with S = dT/dT0 (coolant inlet temperature held fixed).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

R = 8.314462618
P = dict(km=2.0822, Tm=625.0, ER=13636.0, dH=-1.284e6, yx0=0.019, Mwf=0.0303866,
         Pref=101000.0, mdot=2.6371e-3, Rt=0.0125, L=1.5, cp=992.0)
Nf = P['mdot'] / P['Mwf']
Nx0 = P['yx0'] * Nf
Ac = np.pi * P['Rt'] ** 2
NGRID = 3001
ZG = np.linspace(0.0, P['L'], NGRID)


def make_rhs(U, Tc0, Pi, n, kappa):
    Pref, L = P['Pref'], P['L']
    twopiR = 2 * np.pi * P['Rt']
    mcp = P['mdot'] * P['cp']

    def rhs(z, y):
        X, T, p, Tc = y            # p = P/Pref
        k = P['km'] * np.exp(-P['ER'] * (1.0 / T - 1.0 / P['Tm']))
        cx = P['yx0'] * (1.0 - X) * p * Pref / (R * T)
        r = k * p ** (n - 1.0) * cx
        dX = Ac * r / Nx0
        dT = (-P['dH'] * Ac * r - twopiR * U * (T - Tc)) / mcp
        dp = -Pi * (T / Tc0) / p / L          # Eq.(5), reference Tc0 (as in manuscript: Tc)
        dTc = kappa * twopiR * U * (T - Tc) / mcp
        return [dX, dT, dp, dTc]
    return rhs


def integrate(T0, U, Tc0, p0, Pi, n, kappa, rtol=1e-10, dense=True):
    rhs = make_rhs(U, Tc0, Pi, n, kappa)
    sol = solve_ivp(rhs, (0, P['L']), [0.0, T0, p0, Tc0], method='BDF', rtol=rtol,
                    atol=[1e-12, 1e-9, 1e-12, 1e-9], t_eval=ZG if dense else None)
    return sol


def inlet_p(T0, U, Tc0, Pi, n, kappa, mode, plevel):
    """mode: 'const' -> p0=plevel everywhere (Pi forced 0);  'outlet' -> shoot p(L)=1."""
    if mode == 'const':
        return plevel, 0.0
    p_iso = np.sqrt(1 + 2 * Pi)

    def res(p0):
        s = integrate(T0, U, Tc0, p0, Pi, n, kappa, rtol=1e-11, dense=False)
        if s.status != 0 or s.y[2, -1] <= 0:
            return np.nan
        return s.y[2, -1] - 1.0
    lo, hi = 0.97 * p_iso, 1.08 * p_iso
    flo, fhi = res(lo), res(hi)
    if not (np.isfinite(flo) and np.isfinite(fhi) and flo * fhi < 0):
        lo, hi = 0.9 * p_iso, 1.3 * p_iso
        flo, fhi = res(lo), res(hi)
    return brentq(res, lo, hi, xtol=1e-13, rtol=1e-13), Pi


def profile(T0, U, Tc0, Pi, n, kappa, mode, plevel=1.0):
    p0, Pi_eff = inlet_p(T0, U, Tc0, Pi, n, kappa, mode, plevel)
    return integrate(T0, U, Tc0, p0, Pi_eff, n, kappa, rtol=1e-11), p0


def g_value(T0, U, Tc0=625.0, Pi=0.0, n=2.0, kappa=0.0, mode='const', plevel=1.0, h=0.02, full=False):
    """g = Sint - 1, S by central difference in T0 (re-shooting P0 in outlet mode)."""
    sp, p0p = profile(T0 + h, U, Tc0, Pi, n, kappa, mode, plevel)
    sm, p0m = profile(T0 - h, U, Tc0, Pi, n, kappa, mode, plevel)
    s0, p00 = profile(T0, U, Tc0, Pi, n, kappa, mode, plevel)
    if min(sp.status, sm.status, s0.status) != 0 or sp.y.shape[1] != NGRID or sm.y.shape[1] != NGRID:
        return np.nan
    S = (sp.y[1] - sm.y[1]) / (2 * h)
    S[0] = 1.0
    # Sint per Sec. 3.7
    dS0 = S[1] - S[0]
    if dS0 >= 0:
        idx = np.arange(1, NGRID)
        i = idx[np.argmax(S[1:])]
        Sint = _refine(S, i)
    else:
        d = np.diff(S)
        up = np.where((d[:-1] < 0) & (d[1:] >= 0))[0]
        if len(up) == 0:
            Sint = S[-1]
        else:
            imin = up[0] + 1
            i = imin + np.argmax(S[imin:])
            Sint = _refine(S, i)
    if full:
        return Sint - 1.0, dict(S=S, T=s0.y[1], X=s0.y[0], p=s0.y[2], Tc=s0.y[3], p0=p00)
    return Sint - 1.0


def _refine(S, i):
    if 0 < i < len(S) - 1:
        a, b, c = S[i - 1], S[i], S[i + 1]
        den = a - 2 * b + c
        if den < 0:
            return b - 0.125 * (c - a) ** 2 / den
    return S[i]


def ucrit(T0=625.0, Tc0=625.0, Pi=0.0, n=2.0, kappa=0.0, mode='const', plevel=1.0,
          lo=330.0, hi=750.0, tol=0.005):
    f = lambda U: g_value(T0, U, Tc0, Pi, n, kappa, mode, plevel)
    # coarse scan upward from lo to find first +->- sign change (amplified at low U, attenuated at high U)
    Us = np.linspace(lo, hi, 9)
    vals = [f(u) for u in Us]
    for a, b, fa, fb in zip(Us[:-1], Us[1:], vals[:-1], vals[1:]):
        if np.isfinite(fa) and np.isfinite(fb) and fa > 0 > fb:
            # sub-scan
            sub = np.linspace(a, b, 9)
            sv = [f(u) for u in sub]
            for a2, b2, fa2, fb2 in zip(sub[:-1], sub[1:], sv[:-1], sv[1:]):
                if fa2 > 0 > fb2:
                    return brentq(f, a2, b2, xtol=tol)
    return np.nan
