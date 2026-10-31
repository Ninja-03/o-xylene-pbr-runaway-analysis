# Second-dataset single reaction and the two networks at Pi = 0.121 (Ergun reconstruction with the 2 mm particle diameter of the original sources, standard-condition density basis)
# usage: python pi121_run.py ds2 | delasa | cc     (resumable; writes pi121_<model>.json)
import sys, json, time, os
model = sys.argv[1]
PI = 0.121
if model == 'ds2':
    import pbr2
    from pbr2 import *
    pbr2.D['rhoL'] = 1300*A
    HI = 800
else:
    import pbr2, pbr3
    from pbr3 import *
    pbr3.MODEL = model
    if model == 'delasa': pbr2.D['rhoL'] = 1300*pbr2.A
    HI = 1500
yA0 = 1.65/Pref
out = f'pi121_{model}.json'
res = json.load(open(out)) if os.path.exists(out) else []
def save(): json.dump(res, open(out, 'w'), indent=1)
def have(kind, kappa, n=None):
    return any(r['kind'] == kind and abs(r['kappa']-kappa) < 1e-12 and (n is None or r['n'] == n) for r in res)
t0 = time.time()
for kappa in (0.0, KAPPA_DOC):
    if have('U0', kappa): U0 = [r['U0'] for r in res if r['kind'] == 'U0' and abs(r['kappa']-kappa) < 1e-12][0]
    else:
        U0 = ucrit(yA0, kappa=kappa, hi=HI); res.append(dict(kind='U0', kappa=kappa, U0=U0)); save(); print(res[-1], flush=True)
    for n in (1.0, 2.0):
        if have('decomp', kappa, n): continue
        Uc = ucrit(yA0, Pi=PI, n=n, kappa=kappa, mode='outlet', hi=HI)
        g, d = g_value(625, Uc, yA0, Pi=PI, n=n, kappa=kappa, mode='outlet', full=True); p0 = d['p0']
        Ua = ucrit(yA0, n=n, kappa=kappa, mode='const', plevel=p0, hi=HI)
        r = dict(kind='decomp', kappa=kappa, n=n, Pi=PI, U0=U0, p0=p0, UC=Uc, UA=Ua, scale_check=U0*p0**n)
        res.append(r); save(); print(r, f"{time.time()-t0:.0f}s", flush=True)
print('DONE', flush=True)
