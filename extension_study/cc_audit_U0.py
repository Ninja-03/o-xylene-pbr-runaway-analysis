# U_crit,0 of the five-reaction network (kappa = 0 and documented kappa) for the published and corrected source parameters
import json, time
import pbr3
from pbr3 import *
pbr3.MODEL = 'cc'
yA0 = 1.65/Pref
cases = {'published (36, dH5 -87800, ka 0.026)': dict(KSCALE=1.0, dH5=-87800., KA=0.026),
         'dH5 = -88000 only': dict(KSCALE=1.0, dH5=-88000., KA=0.026),
         'exact conversion only': dict(KSCALE=(3600/101.325)/36, dH5=-87800., KA=0.026),
         'all C&C-consistent': dict(KSCALE=(3600/101.325)/36, dH5=-88000., KA=7.224e-6*3600)}
out = []
for name, c in cases.items():
    pbr3.KSCALE, pbr3.KA = c['KSCALE'], c['KA']; pbr3.DH['c35'] = c['dH5']
    for kappa in (0.0, KAPPA_DOC):
        t = time.time(); U0 = ucrit(yA0, kappa=kappa, hi=1500)
        r = dict(case=name, kappa=kappa, U0_kcal=U0, U0_W=U0*4184/3600); out.append(r); print(r, f'{time.time()-t:.0f}s', flush=True)
        json.dump(out, open('cc_audit_U0.json', 'w'), indent=1)
