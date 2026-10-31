# Second dataset against the hot-spot excursion reported in van Welsenaere & Froment (1970), Fig. 3:
# excursion between T0 = Tw = 627 and 628 K at p0(o-xylene) = 0.017 atm.  Constant Tc = T0, no pressure drop, fixed oxygen (n = 1).  (~1 min)
import numpy as np, pbr2
from pbr2 import *
pA0 = 0.017*101.325            # kPa
yA0 = pA0/Pref
print(f"pA0 = {pA0:.4f} kPa (0.017 atm)")
for label, rho in (("rho_L = 0.6381", 1300*A), ("rho_L = 0.5203", 0.5203)):
    pbr2.D['rhoL'] = rho
    print(label)
    prev = None
    for T0 in np.arange(624.0, 633.01, 1.0):
        s = pbr2.integ(T0, 82.7, T0, 1.0, 0.0, 1.0, 0.0, yA0)       # Tc0 = T0 (Tw = T0)
        if s.status != 0:
            print(f"   T0 = Tw = {T0:5.1f} K: integration failed"); continue
        Tm = s.y[1].max(); X = s.y[0][-1]
        jump = '' if prev is None else f"  (change in Tmax {Tm-prev:+7.1f} K)"
        print(f"   T0 = Tw = {T0:5.1f} K: Tmax = {Tm:8.1f} K   X(L) = {X:.3f}{jump}")
        prev = Tm

# fine scan of the excursion with rho_L = 0.6381 (location of the steepest rise of Tmax with T0)
pbr2.D['rhoL'] = 1300*A
Ts = np.arange(626.5, 629.01, 0.1)
Tm = []
for T0 in Ts:
    s = pbr2.integ(T0, 82.7, T0, 1.0, 0.0, 1.0, 0.0, yA0)
    Tm.append(s.y[1].max())
Tm = np.array(Tm)
slope = np.diff(Tm)/np.diff(Ts)
k = int(np.argmax(slope))
print(f"rho_L = 0.6381: steepest rise of Tmax with T0 at T0 = Tw = {0.5*(Ts[k]+Ts[k+1]):.2f} K (dTmax/dT0 = {slope[k]:.0f});")
for T0, T in zip(Ts[::5], Tm[::5]): print(f"   T0 = {T0:6.2f} K  Tmax = {T:7.1f} K")
# first T0 at which Tmax exceeds Tw + 100 K (a simple excursion marker)
idx = np.where(Tm > Ts + 100)[0]
print("first T0 with Tmax > T0 + 100 K:", f"{Ts[idx[0]]:.1f} K" if len(idx) else "none in 626.5-629")
