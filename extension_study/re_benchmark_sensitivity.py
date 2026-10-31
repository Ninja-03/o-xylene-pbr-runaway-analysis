# Rawlings-Ekerdt benchmark: which input would have to change, and by how much, to move the excursion from ~627 K to the printed ~631 K?
# Excursion marker: feed temperature T0 at which Tmax reaches 800 K (Tc = 625 K, U = 373 W m-2 K-1, no pressure drop).
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
R = 8.314462618
BASE = dict(km=2.0822, Tm=625.0, ER=13636.0, dH=-1.284e6, yx0=0.019, Mwf=0.0303866, Pref=101000.0,
            mdot=2.6371e-3, Rt=0.0125, L=1.5, cp=992.0, U=373.0)

def tmax(T0, q, Tc=625.0):
    Nf = q['mdot']/q['Mwf']; Nx0 = q['yx0']*Nf; Ac = np.pi*q['Rt']**2
    def rhs(z, y):
        X, T = y
        k = q['km']*np.exp(-q['ER']*(1/T - 1/q['Tm']))
        cx = q['yx0']*(1-X)*q['Pref']/(R*T)
        r = k*cx
        return [Ac*r/Nx0, (-q['dH']*Ac*r - 2*np.pi*q['Rt']*q['U']*(T-Tc))/(q['mdot']*q['cp'])]
    s = solve_ivp(rhs, (0, q['L']), [0., T0], method='BDF', rtol=1e-10, atol=[1e-12, 1e-9])
    return s.y[1].max() if s.status == 0 else np.inf

def t0_exc(q, level=800.0):
    return brentq(lambda T0: tmax(T0, q) - level, 560.0, 700.0, xtol=1e-4)

print("check vs manuscript Table 4:", [round(tmax(T, BASE), 1) for T in (625., 627., 628.)], "(manuscript 718.9, 828.8, 917.7)")
t_base = t0_exc(BASE)
print(f"base excursion marker T0(Tmax = 800 K) = {t_base:.2f} K ; target 631 K -> shift {631 - t_base:+.2f} K\n")
print(f"{'parameter':10s} {'dT0*/dln(x) (K)':>18s} {'factor needed for 631 K':>26s} {'relative change':>16s}")
rows = []
for name in ('km', 'ER', 'dH', 'yx0', 'mdot', 'Rt', 'L', 'cp', 'U'):
    def f(fac):
        q = dict(BASE); q[name] = BASE[name]*fac
        return t0_exc(q) - 631.0
    q1 = dict(BASE); q1[name] = BASE[name]*1.002; q2 = dict(BASE); q2[name] = BASE[name]*0.998
    sens = (t0_exc(q1) - t0_exc(q2))/(np.log(1.002) - np.log(0.998))
    try:
        fac = brentq(f, 0.985, 1.04 if name == 'cp' else 1.015, xtol=1e-7); need = f"{fac:.4f}"; rel = f"{100*(fac-1):+.1f} %"
    except ValueError:
        need = "no effect (hot spot interior)"; rel = "-"
    rows.append((name, sens, need, rel)); print(f"{name:10s} {sens:18.2f} {need:>26s} {rel:>16s}")
