# Re-runs the coolant T0-scan quoted in Section 4.12
# (co-current coolant, kappa = 0.02, i.e. 16.2 K coolant rise at full conversion)
from pbr import ucrit

KAPPA = 0.02
manuscript = {610.0: 389.27, 625.0: 400.10, 640.0: 423.06}   # W m-2 K-1

for T0, ref in manuscript.items():
    U = ucrit(T0=T0, Tc0=625.0, kappa=KAPPA, lo=330, hi=1500)
    status = "OK" if abs(U - ref) < 0.01 else "CHECK"
    print(f"T0 = {T0:.0f} K   Ucrit = {U:.3f}   manuscript = {ref:.2f}   diff = {U-ref:+.3f}   {status}")
