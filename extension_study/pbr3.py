"""
Multi-step networks of Zuluaga-Botero et al. (2020), Table 3, per tube, units kcal-kmol-h-m-K-kPa.
  MODEL='delasa': A->B (k21), B->COx (k22), A->COx (k23)           [3 reactions]
  MODEL='cc'    : Chandrasekharan-Calderbank 5 reactions (k31..k35), site coverage theta_v
Assumptions: total molar flow constant (dilute feed), P_i = y_i * P, y_i = yA0 * (F_i/F_A0);
             oxygen absorbed in the rate constants (fixed, as in the source); n=2 variant multiplies all rates by p.
State: [X (conversion of A), T, p, Tc, fB, fD, fE]  (f_i = F_i / F_A0)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import pbr2
from pbr2 import D, A, mdot, Ft, mcp, KAPPA_DOC, Pref, NGRID, ZG

MODEL = 'delasa'
GSCALE = 1.0          # multiplies alpha*rho_b*F_t of the 5-reaction rate (interpretation test)
DH22_DELASA = -783700.0          # Hess closure (dH23 - dH21); printed value -7,837,000 is exactly 10x this
DH = dict(d21=-307122., d22=DH22_DELASA, d23=-1090822.,
          c31=-108999., c32=-295999., c33=-550000., c34=-99001., c35=-87800.)
ALPHA, RHOB = 3.2832e-3, 1300.0
# --- five-reaction source audit (Chandrasekharan & Calderbank 1979, Table 3, CO case) ---
# as published by Zuluaga-Botero et al. (2020): pre-exponentials = 36 x C&C values (3600 s/h over 100 kPa/atm), dH35 = -87800, kaPO2 = 0.026
# C&C-consistent variant: KSCALE = (3600/101.325)/36, DH['c35'] = -88000, KA = 7.224e-6*3600 = 0.0260064
KSCALE = 1.0
KA = 0.026


def rhs_factory(U, Tc0, Pi, n, kappa, yA0):
    FA0 = Ft*yA0
    pip = np.pi*D['Dt']
    def rhs(z, y):
        X, T, p, Tc, fB, fD, fE = y
        pf = p**(n-1.0)
        P = p*Pref
        PA = yA0*(1-X)*P; PB = yA0*fB*P; PD = yA0*fD*P; PE = yA0*fE*P
        if MODEL == 'delasa':
            k21 = 21.07*np.exp(10.59-13587.68/T); k22 = 21.07*np.exp(11.62-15801.97/T); k23 = 21.07*np.exp(9.73-14392.88/T)
            r21 = D['rhoL']*pf*k21*PA; r22 = D['rhoL']*pf*k22*PB; r23 = D['rhoL']*pf*k23*PA
            dX = (r21+r23)/FA0; dB = (r21-r22)/FA0; dD = 0.; dE = 0.
            heat = -(DH['d21']*r21 + DH['d22']*r22 + DH['d23']*r23)
        else:
            k31 = KSCALE*137.808*np.exp(-7381.11/T); k32 = KSCALE*46.62*np.exp(-6552.43/T); k33 = KSCALE*12.888*np.exp(-6153.13/T)
            k34 = KSCALE*20.016*np.exp(-5585.44/T); k35 = KSCALE*114.84*np.exp(-6962.56/T)
            th = KA/(KA + (k31+3*k32+6.5*k33)*PA + k34*PD + k35*PE)
            G = th*ALPHA*RHOB*Ft*pf*GSCALE
            r31 = G*k31*PA; r32 = G*k32*PA; r33 = G*k33*PA; r34 = G*k34*PD; r35 = G*k35*PE
            dX = (r31+r32+r33)/FA0; dB = (r32+r35)/FA0; dD = (r31-r34)/FA0; dE = (r34-r35)/FA0
            heat = -(DH['c31']*r31 + DH['c32']*r32 + DH['c33']*r33 + DH['c34']*r34 + DH['c35']*r35)
        dT = (heat - pip*U*(T-Tc))/mcp
        dp = -Pi*(T/Tc0)/p/D['L']
        dTc = kappa*pip*U*(T-Tc)/mcp
        return [dX, dT, dp, dTc, dB, dD, dE]
    return rhs


def integ(T0, U, Tc0, p0, Pi, n, kappa, yA0, rtol=1e-11, dense=True):
    return solve_ivp(rhs_factory(U, Tc0, Pi, n, kappa, yA0), (0, D['L']), [0., T0, p0, Tc0, 0., 0., 0.], method='BDF',
                     rtol=rtol, atol=[1e-12, 1e-9, 1e-12, 1e-9, 1e-12, 1e-12, 1e-12], t_eval=ZG if dense else None)

# re-bind the pbr2 machinery to the network integrator (same criterion, same procedures)
pbr2.integ = integ
g_value = pbr2.g_value
ucrit = pbr2.ucrit
crit_pA0 = pbr2.crit_pA0

# ---- faster settings for the network runs (validated against the rtol=1e-11 reference, see validate_fast.py)
RTOL = 1e-10
_integ_raw = integ
def integ(T0, U, Tc0, p0, Pi, n, kappa, yA0, rtol=None, dense=True):
    return _integ_raw(T0, U, Tc0, p0, Pi, n, kappa, yA0, rtol=RTOL if rtol is None else rtol, dense=dense)

def inlet_p_fast(T0, U, Tc0, Pi, n, kappa, yA0, mode, plevel):
    if mode == 'const':
        return plevel, 0.0
    piso = np.sqrt(1+2*Pi)
    def res(p0):
        s = integ(T0, U, Tc0, p0, Pi, n, kappa, yA0, rtol=1e-9, dense=False)
        return np.nan if (s.status != 0 or s.y[2, -1] <= 0) else s.y[2, -1]-1.0
    lo, hi = 0.995*piso, 1.03*piso
    flo, fhi = res(lo), res(hi)
    if not (np.isfinite(flo) and np.isfinite(fhi) and flo*fhi < 0):
        lo, hi = 0.9*piso, 1.3*piso
    return brentq(res, lo, hi, xtol=1e-10, rtol=1e-10), Pi

pbr2.integ = integ
pbr2.inlet_p = inlet_p_fast
