# Interpretation test for the 5-reaction network: scale alpha*rho_b*F_t by 0.5 and 2 (kappa = 0, n = 2)
import json, time
import pbr3
from pbr3 import *
pbr3.MODEL = 'cc'
yA0 = 1.65/Pref
res = []
for sc in (0.5, 2.0):
    pbr3.GSCALE = sc
    U0 = ucrit(yA0, hi=1500)
    print('scale', sc, 'U0 kcal', U0, 'W', U0*4184/3600, flush=True)
    for Pi in (0.011, 0.2):
        Uc = ucrit(yA0, Pi=Pi, n=2.0, mode='outlet', hi=1500)
        g, d = g_value(625, Uc, yA0, Pi=Pi, n=2.0, mode='outlet', full=True); p0 = d['p0']
        Ua = ucrit(yA0, n=2.0, mode='const', plevel=p0, hi=1500)
        r = dict(scale=sc, kappa=0.0, n=2.0, Pi=Pi, U0=U0, p0=p0, UC=Uc, UA=Ua); res.append(r); print(r, flush=True)
        json.dump(res, open('cc_scale_results.json', 'w'), indent=1)
