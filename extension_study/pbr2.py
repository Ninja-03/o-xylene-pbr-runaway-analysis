"""
Dataset 2: van Welsenaere-Froment single reaction as tabulated in Zuluaga-Botero et al. (2020) Tables 2-3.
Per tube; units kcal, kmol, h, m, K, kPa.  Rate per unit length:
   r_L = rhoL * exp(19.837 - a/T) * (PO2/101.325) * (pA/101.325)       [kmol/(h m)]
   (printed k1 = 21.07 e^19.837 with kPa pressures is inconsistent by 101.325^2; corrected here, see text)
   PO2 = 21.07 * p^(n-1) kPa   (n=1: fixed O2 as in source; n=2: O2 scales with total pressure)
   pA  = yA0 (1-X) P  (kPa), P = p*101.325
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

Pref = 101.325
D = dict(Dt=0.025, u=3600., rho=1.293, cp=0.2498, M=29.48, rhoL=0.5203, a=13636., dH=-306707.,
         U=82.7, L=2.0, wc=72000/3000, cpc=0.3105)
A = np.pi/4*D['Dt']**2
mdot = D['u']*D['rho']*A
Ft = mdot/D['M']
mcp = mdot*D['cp']
KAPPA_DOC = mcp/(D['wc']*D['cpc'])
NGRID = 3001
ZG = np.linspace(0, D['L'], NGRID)


def rhs_factory(U, Tc0, Pi, n, kappa, yA0):
    FA0 = Ft*yA0
    pip = np.pi*D['Dt']
    def rhs(z, y):
        X, T, p, Tc = y
        k = np.exp(19.837 - D['a']/T)
        pA = yA0*(1-X)*p*Pref
        PO2 = 21.07*p**(n-1)
        rL = D['rhoL']*k*(PO2/Pref)*(pA/Pref)
        dX = rL/FA0
        dT = (-D['dH']*rL - pip*U*(T-Tc))/mcp
        dp = -Pi*(T/Tc0)/p/D['L']
        dTc = kappa*pip*U*(T-Tc)/mcp
        return [dX, dT, dp, dTc]
    return rhs


def integ(T0, U, Tc0, p0, Pi, n, kappa, yA0, rtol=1e-11, dense=True):
    return solve_ivp(rhs_factory(U, Tc0, Pi, n, kappa, yA0), (0, D['L']), [0., T0, p0, Tc0], method='BDF',
                     rtol=rtol, atol=[1e-12, 1e-9, 1e-12, 1e-9], t_eval=ZG if dense else None)


def inlet_p(T0, U, Tc0, Pi, n, kappa, yA0, mode, plevel):
    if mode == 'const':
        return plevel, 0.0
    piso = np.sqrt(1+2*Pi)
    def res(p0):
        s = integ(T0, U, Tc0, p0, Pi, n, kappa, yA0, dense=False)
        return np.nan if (s.status != 0 or s.y[2, -1] <= 0) else s.y[2, -1]-1.0
    lo, hi = 0.97*piso, 1.08*piso
    if not (res(lo)*res(hi) < 0):
        lo, hi = 0.9*piso, 1.3*piso
    return brentq(res, lo, hi, xtol=1e-13, rtol=1e-13), Pi


def prof(T0, U, Tc0, Pi, n, kappa, yA0, mode, plevel):
    p0, Pie = inlet_p(T0, U, Tc0, Pi, n, kappa, yA0, mode, plevel)
    return integ(T0, U, Tc0, p0, Pie, n, kappa, yA0), p0


def _refine(S, i):
    if 0 < i < len(S)-1:
        a, b, c = S[i-1], S[i], S[i+1]; den = a-2*b+c
        if den < 0: return b-0.125*(c-a)**2/den
    return S[i]


def g_value(T0, U, yA0, Tc0=625., Pi=0., n=1., kappa=0., mode='const', plevel=1., h=0.02, full=False):
    sp, _ = prof(T0+h, U, Tc0, Pi, n, kappa, yA0, mode, plevel)
    sm, _ = prof(T0-h, U, Tc0, Pi, n, kappa, yA0, mode, plevel)
    s0, p00 = prof(T0, U, Tc0, Pi, n, kappa, yA0, mode, plevel)
    if min(sp.status, sm.status, s0.status) != 0 or sp.y.shape[1] != NGRID or sm.y.shape[1] != NGRID:
        return np.nan
    S = (sp.y[1]-sm.y[1])/(2*h); S[0] = 1.
    if S[1]-S[0] >= 0:
        i = 1+np.argmax(S[1:]); Sint = _refine(S, i)
    else:
        d = np.diff(S); up = np.where((d[:-1] < 0) & (d[1:] >= 0))[0]
        if len(up) == 0: Sint = S[-1]
        else:
            im = up[0]+1; i = im+np.argmax(S[im:]); Sint = _refine(S, i)
    if full: return Sint-1., dict(S=S, T=s0.y[1], X=s0.y[0], p=s0.y[2], p0=p00)
    return Sint-1.


def ucrit(yA0, T0=625., Tc0=625., Pi=0., n=1., kappa=0., mode='const', plevel=1., lo=10., hi=400., tol=0.002):
    f = lambda U: g_value(T0, U, yA0, Tc0, Pi, n, kappa, mode, plevel)
    Us = np.linspace(lo, hi, 9); v = [f(u) for u in Us]
    for a, b, fa, fb in zip(Us[:-1], Us[1:], v[:-1], v[1:]):
        if np.isfinite(fa) and np.isfinite(fb) and fa > 0 > fb:
            sub = np.linspace(a, b, 9); sv = [f(u) for u in sub]
            for a2, b2, fa2, fb2 in zip(sub[:-1], sub[1:], sv[:-1], sv[1:]):
                if fa2 > 0 > fb2: return brentq(f, a2, b2, xtol=tol)
    return np.nan


def crit_pA0(U=82.7, kappa=0., n=1., lo=0.8, hi=3.0, T0=625., Tc0=625.):
    """critical inlet o-xylene partial pressure (kPa) at fixed U: g changes - -> + as pA0 rises"""
    f = lambda pa: g_value(T0, U, pa/Pref, Tc0, 0., n, kappa)
    return brentq(f, lo, hi, xtol=1e-5)
