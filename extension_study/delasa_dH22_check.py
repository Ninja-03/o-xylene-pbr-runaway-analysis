# de Lasa network: Hess-law check of the printed dH22 and its effect
# (feed 1.65 kPa, T0 = Tc = 625 K, no pressure drop, constant Tc, rho_L = rho_b*pi*Dt^2/4 = 0.6381 kg/m)   (~1 min)
import numpy as np, pbr2, pbr3
from pbr3 import *
pbr2.D['rhoL'] = 1300*pbr2.A
pbr3.MODEL = 'delasa'
yA0 = 1.65/Pref
print("rho_L =", pbr2.D['rhoL'])
print("Hess closure: dH23 - dH21 =", -1090822 - (-307122), " printed dH22 = -7837000  (ratio %.1f)" % (-7837000/(-1090822 + 307122)))
for label, d22 in (("Hess value -783700", -783700.), ("printed value -7837000", -7837000.)):
    pbr3.DH['d22'] = d22
    print(label)
    for U in (120, 125, 150, 175, 200, 225, 250):
        s = pbr3.integ(625., U, 625., 1.0, 0., 1., 0., yA0)
        if s.status != 0:
            print(f"   U = {U} kcal m-2 h-1 K-1: integration fails (status {s.status}) at z = {s.t[-1]:.3f} m, Tmax so far {s.y[1].max():.0f} K")
        else:
            print(f"   U = {U} kcal m-2 h-1 K-1: g = {g_value(625, U, yA0):.3f}")
