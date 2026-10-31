# Critical inlet o-xylene partial pressure of the second dataset at the tabulated U (82.7 kcal m-2 h-1 K-1), for both catalyst loadings
# and both coolant treatments.   Compare: van Welsenaere & Froment 1.651 kPa; Henning & Perez 1.677-1.815 kPa.  (~1 min)
import numpy as np, pbr2
from pbr2 import *
for label, rho in (("rho_L = 0.6381 (rho_b*pi*Dt^2/4)", 1300*A), ("rho_L = 0.5203 (Zuluaga-Botero Table 3)", 0.5203)):
    pbr2.D['rhoL'] = rho
    for name, kap, lo, hi in (("constant Tc", 0.0, 1.0, 2.6), ("co-current, kappa=0.0766", KAPPA_DOC, 0.8, 2.3)):
        pc = crit_pA0(lo=lo, hi=hi, kappa=kap)
        print(f"{label:42s} {name:26s} critical pA0 = {pc:.3f} kPa", flush=True)
