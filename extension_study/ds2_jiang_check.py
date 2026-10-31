# Second dataset vs Jiang, Jiang & Pan (2009), CIESC Journal 60(10), 2490-2496: Table 2 (critical o-xylene partial pressure, kPa,
# T0 = Tc = 623, 625, 628, 630 K; van Welsenaere-Froment, Morbidelli-Varma and divergence criteria) and Fig. 4 groups
# (gamma = 21.82, Da = 0.47, St = 22.8, tau = 0.0317).  Our threshold criterion g = 0, constant Tc and co-current coolant.  (~2 min)
import numpy as np, pbr2
from pbr2 import *
jiang = {623: (1.77, 1.94, 1.79), 625: (1.68, 1.85, 1.69), 628: (1.56, 1.73, 1.57), 630: (1.48, 1.65, 1.49)}   # VF, MV, div
print("groups at pA0 = 1.65 kPa, T0 = Tc = 625 K")
for label, rho in (("rho_L = 0.6381", 1300*A), ("rho_L = 0.5203", 0.5203)):
    pbr2.D['rhoL'] = rho
    yA0 = 1.65/Pref
    k = np.exp(19.837 - pbr2.D['a']/625.); rL = rho*k*(21.07/Pref)*(1.65/Pref); Da = pbr2.D['L']*rL/(Ft*yA0)
    print(f"  {label}: Da = {Da:.3f}  (Jiang 0.47)   St(82.7) = {np.pi*pbr2.D['Dt']*82.7*pbr2.D['L']/mcp:.2f} (22.8)   gamma = {pbr2.D['a']/625.:.2f} (21.82)   tau (Jiang, W_co = 15 kg/s) = {np.pi*0.025**2/4*3000*1.0*1.293*1.044/(4.18*15):.4f} (0.0317 as printed)")
print("\ncritical inlet partial pressure of o-xylene (kPa), T0 = Tc, U = 82.7 kcal m-2 h-1 K-1, n = 1")
print(f"{'T0 (K)':>7s} | {'Jiang: VF / MV / div':>22s} | {'ours rho_L=0.6381':>18s} | {'ours rho_L=0.5203':>18s}")
res = {}
for T0, ref in jiang.items():
    row = []
    for rho in (1300*A, 0.5203):
        pbr2.D['rhoL'] = rho
        row.append(crit_pA0(U=82.7, kappa=0.0, n=1.0, lo=1.0, hi=3.0, T0=float(T0), Tc0=float(T0)))
    res[T0] = row
    print(f"{T0:7d} | {ref[0]:6.2f} / {ref[1]:4.2f} / {ref[2]:4.2f}   | {row[0]:18.3f} | {row[1]:18.3f}")
# co-current coolant with Jiang's capacity ratio (tau = 0.0317) at T0 = 625 K; Jiang: div criterion with coolant 1.64 kPa (1.69 without)
pbr2.D['rhoL'] = 1300*A
kj = 0.0317
pc = crit_pA0(U=82.7, kappa=kj, n=1.0, lo=1.0, hi=3.0, T0=625., Tc0=625.)
print(f"\nco-current coolant, kappa = tau = {kj}: ours {pc:.3f} kPa  (Jiang modified-div 1.64 kPa; unmodified div 1.69 kPa)")
