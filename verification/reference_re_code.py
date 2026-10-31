"""V9 external reference: transcription of the published Rawlings & Ekerdt companion code
(main.py for Figs 6.35/6.36, 'Converted from xyleneT.m', GNU GPL), accessed from
https://sites.engineering.ucsb.edu/~jbraw/chemreacfun/fig-html/ch6/fig-6-36.html

SOURCE-DERIVED: every literal below is the published script's own value, in the script's units
(kmol, kJ, kN/m^2, s, m). This file is used ONLY by test V9 as the independent comparator; no model
module imports it. Only the file-output line (save_ascii) is omitted.
"""
import numpy as np
from scipy.integrate import solve_ivp

R_tube = 1.25e-2      # tube radius, m
Ac     = np.pi*R_tube**2
Qrho   = 0.0026371    # mass flow, kg/sec
Pf     = 1.01e2       # feed pressure, kN/m^2
Mwf    = 0.98*(0.79*28+0.21*32)+0.02*106.17  # mol wt feed, kg/kmol
Tf     = 625.         # feed temperature, K
Rg     = 8.314        # kJ/(K kmol)
Nf     = Qrho / Mwf   # kmol/sec
E      = 13636.       # K
Tm     = 625.         # K
km     = 2.0822       # 1/sec
Ta     = 625.         # K
Cp     = 0.992        # kJ/kg K
delH   = -1.284e6     # kJ/kmol
U      = 0.373        # kJ/(m^2 sec K)
beta   = delH * Ac / (Qrho*Cp)
Gamma  = 2*np.pi*R_tube*U / (Qrho*Cp)
l      = 1.5          # m

p = dict(Nf=Nf, Pf=Pf, Rg=Rg, E=E, Tm=Tm, km=km, Ta=Ta, beta=beta, Gamma=Gamma, Ac=Ac)


def pfr(t, x, p):
    Nx = x[0]; T = x[1]
    Q  = p['Nf'] / (p['Pf']/(p['Rg']*T))
    cx = Nx / Q
    k  = p['km'] * np.exp(-p['E']*(1./T - 1./p['Tm']))
    rate = k * cx
    return [-rate*p['Ac'],
            -p['beta']*rate + p['Gamma']*(p['Ta']-T)]


npts   = 200
z      = np.linspace(0., l, npts)
yxfeed = 0.019
Nxf    = yxfeed * Nf
Tfeed  = [615., 620., 625., 630.]


def run(Tf_i, rtol=None, atol=None, method="Radau"):
    """Published settings by default: Radau, rtol = atol = sqrt(eps), t_eval = z (200 points)."""
    eps_sq = np.sqrt(np.finfo(float).eps)
    sol = solve_ivp(lambda t, x: pfr(t, x, p), [0., l], [Nxf, Tf_i], method=method, t_eval=z,
                    rtol=eps_sq if rtol is None else rtol, atol=eps_sq if atol is None else atol)
    return z, sol.y[0, :], sol.y[1, :], sol
