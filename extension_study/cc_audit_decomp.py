# Pressure decomposition of the five-reaction network with C&C-consistent source parameters (resumable)
#   exact atm->kPa conversion of the pre-exponentials, dH35 = -88000 kcal/kmol, kaPO2 = 7.224e-6*3600
import json, time, os
import pbr3
from pbr3 import *
pbr3.MODEL = 'cc'; pbr3.KSCALE = (3600/101.325)/36; pbr3.KA = 7.224e-6*3600; pbr3.DH['c35'] = -88000.
yA0 = 1.65/Pref
out = 'cc_audit_decomp.json'
res = json.load(open(out)) if os.path.exists(out) else []
def save(): json.dump(res, open(out, 'w'), indent=1)
def have(kind, kappa, n=None, Pi=None):
    return any(r['kind'] == kind and abs(r['kappa']-kappa) < 1e-12 and (n is None or r['n'] == n) and (Pi is None or r['Pi'] == Pi) for r in res)
t0 = time.time()
for kappa in (0.0, KAPPA_DOC):
    if have('U0', kappa): U0 = [r['U0'] for r in res if r['kind'] == 'U0' and abs(r['kappa']-kappa) < 1e-12][0]
    else:
        U0 = ucrit(yA0, kappa=kappa, hi=1500); res.append(dict(kind='U0', kappa=kappa, U0=U0)); save(); print(res[-1], flush=True)
    for Pi in (0.011, 0.2):
        for n in (1.0, 2.0):
            if have('decomp', kappa, n, Pi): continue
            Uc = ucrit(yA0, Pi=Pi, n=n, kappa=kappa, mode='outlet', hi=1500)
            g, d = g_value(625, Uc, yA0, Pi=Pi, n=n, kappa=kappa, mode='outlet', full=True); p0 = d['p0']
            Ua = ucrit(yA0, n=n, kappa=kappa, mode='const', plevel=p0, hi=1500)
            r = dict(kind='decomp', kappa=kappa, n=n, Pi=Pi, U0=U0, p0=p0, UC=Uc, UA=Ua, scale_check=U0*p0**n)
            res.append(r); save(); print(r, f"{time.time()-t0:.0f}s", flush=True)
print('DONE', flush=True)
