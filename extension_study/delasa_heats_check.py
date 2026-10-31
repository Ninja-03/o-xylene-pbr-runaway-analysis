# de Lasa network: Zuluaga-Botero heats (Hess-corrected dH22) versus the original de Lasa (1983) heats.
# feed 1.65 kPa, T0 = Tc = 625 K, no pressure drop, constant and documented co-current coolant; rho_L = 0.6381 kg/m.  (~1 min)
import pbr2, pbr3
from pbr3 import *
pbr2.D['rhoL'] = 1300*pbr2.A
pbr3.MODEL = 'delasa'
print("rho_L =", pbr2.D['rhoL'], "kg/m   (expect 0.6381 for the primary loading)")
yA0 = 1.65/Pref; W = 4184/3600
sets = {"Zuluaga-Botero (-307122, -783700, -1090822)": dict(d21=-307122., d22=-783700., d23=-1090822.),
        "de Lasa 1983     (-307000, -783000, -1090000)": dict(d21=-307000., d22=-783000., d23=-1090000.)}
res = {}
for label, h in sets.items():
    pbr3.DH.update(h)
    for kappa in (0.0, KAPPA_DOC):
        U = ucrit(yA0, kappa=kappa, hi=1500)
        res[(label, kappa)] = U
        print(f"{label}  kappa={kappa:.4f}  Ucrit,0 = {U:.4f} kcal = {U*W:.2f} W m-2 K-1", flush=True)
labels = list(sets)
for kappa in (0.0, KAPPA_DOC):
    a, b = res[(labels[0], kappa)], res[(labels[1], kappa)]
    print(f"kappa={kappa:.4f}: change {100*(b-a)/a:+.3f} %")
