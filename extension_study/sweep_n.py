import json, numpy as np
from pbr import *
res=[]
for n in (1.25,1.5,1.75):
    U0=ucrit(n=n)
    for Pi in (0.1,0.2,0.35):
        Uc=ucrit(Pi=Pi,n=n,mode='outlet')
        g,d=g_value(625,Uc,Pi=Pi,n=n,mode='outlet',full=True); p0=d['p0']
        Ua=ucrit(n=n,mode='const',plevel=p0)
        r=dict(n=n,Pi=Pi,U0=U0,p0=p0,UC=Uc,UA=Ua,scale_check=U0*p0**n)
        res.append(r); print(r,flush=True)
        json.dump(res,open('sweep_n.json','w'),indent=1)
