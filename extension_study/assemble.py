import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

U0 = 380.3823690189159          # kappa=0, n-independent (p=1 -> same boundary)
dTad = 809.329

# --- load results -------------------------------------------------------------------
rows = []
# n = 2 and 1 from t2 (Table-7 reproduction), recompute from printed values (precise values via json below)
base = json.load(open('base_n12.json'))        # created below from t2 rerun values
for r in base: rows.append(r)
for r in json.load(open('sweep_n.json')):
    rows.append(dict(kappa=0.0, n=r['n'], Pi=r['Pi'], U0=r['U0'], p0=r['p0'], UC=r['UC'], UA=r['UA']))
for r in json.load(open('coolant_pressure.json')):
    if r['kind'] == 'decomp':
        rows.append(dict(kappa=r['kappa'], n=r['n'], Pi=r['Pi'], U0=r['U0'], p0=r['p0'], UC=r['UC'], UA=r['UA']))
for r in rows:
    r['level'] = r['UA'] - r['U0']
    r['decl'] = r['UC'] - r['UA']
    r['comb'] = r['UC'] - r['U0']
    r['peff'] = (r['UC'] / r['U0']) ** (1.0 / r['n'])
    r['peff_frac'] = (r['peff'] - 1) / (r['p0'] - 1)
    r['scale_err'] = r['UA'] - r['U0'] * r['p0'] ** r['n']
    r['decl_over_level'] = r['decl'] / r['level']
json.dump(rows, open('all_rows.json', 'w'), indent=1)
for r in sorted(rows, key=lambda r: (r['kappa'], r['n'], r['Pi'])):
    print(f"k={r['kappa']:<5} n={r['n']:<4} Pi={r['Pi']:<4} p0={r['p0']:.4f} lev={r['level']:+8.2f} dec={r['decl']:+7.2f} comb={r['comb']:+8.2f} "
          f"lev/|dec|={-r['level']/r['decl']:5.2f} peff={r['peff']:.4f} frac={r['peff_frac']:.3f} scaleErr={r['scale_err']:+.4f}")
