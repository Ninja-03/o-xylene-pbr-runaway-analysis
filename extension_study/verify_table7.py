# Checks that pbr.py reproduces the manuscript's own Table 7 and base point
from pbr import ucrit, g_value

U0 = ucrit()
print(f"Ucrit,0 = {U0:.3f}   manuscript 380.38")
g, d = g_value(625, 373, full=True)
print(f"Tmax = {d['T'].max():.3f}   manuscript 718.879")

# (n, Pi, Ucrit(C), Ucrit(A'))  from manuscript Table 7
table7 = [(2.0, 0.10, 446.60, 459.37), (2.0, 0.20, 515.93, 537.63), (2.0, 0.35, 623.13, 654.25),
          (1.0, 0.10, 411.66, 418.14), (1.0, 0.20, 441.49, 452.57), (1.0, 0.35, 483.66, 499.59)]

for n, Pi, refC, refA in table7:
    Uc = ucrit(Pi=Pi, n=n, mode='outlet')
    _, dd = g_value(625, Uc, Pi=Pi, n=n, mode='outlet', full=True)
    Ua = ucrit(n=n, mode='const', plevel=dd['p0'])
    ok = abs(Uc - refC) < 0.03 and abs(Ua - refA) < 0.03
    print(f"n={n:g} Pi={Pi}:  C = {Uc:.2f} (ms {refC})   A' = {Ua:.2f} (ms {refA})   {'OK' if ok else 'CHECK'}", flush=True)
