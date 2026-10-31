import json, numpy as np
from pbr2 import *
yA0=1.65/Pref
res=[]
def save(): json.dump(res,open('ds2_results.json','w'),indent=1)
for kappa in (0.0,KAPPA_DOC):
    U0=ucrit(yA0,kappa=kappa,hi=800)
    res.append(dict(kind='U0',kappa=kappa,U0=U0)); save(); print(res[-1],flush=True)
    for n in (1.0,2.0):
        for Pi in (0.011,0.1,0.2,0.35):
            Uc=ucrit(yA0,Pi=Pi,n=n,kappa=kappa,mode='outlet',hi=800)
            g,d=g_value(625,Uc,yA0,Pi=Pi,n=n,kappa=kappa,mode='outlet',full=True); p0=d['p0']
            Ua=ucrit(yA0,n=n,kappa=kappa,mode='const',plevel=p0,hi=800)
            r=dict(kind='decomp',kappa=kappa,n=n,Pi=Pi,U0=U0,p0=p0,UC=Uc,UA=Ua,scale_check=U0*p0**n)
            res.append(r); save(); print(r,flush=True)
