import json, numpy as np
from pbr import *
res=[]
dTad=809.329
# boundary vs T0 under coolant warming
for kappa in ():
    for T0 in (610.0,625.0,640.0):
        res.append(dict(kind='T0scan',kappa=kappa,T0=T0,Ucrit=ucrit(T0=T0,kappa=kappa,lo=330,hi=1500)))
        print(res[-1],flush=True)
json.dump(res,open('coolant_pressure.json','w'),indent=1)
for kappa in (0.02,0.05):
    for n in (2.0,1.0):
        U0=ucrit(n=n,kappa=kappa,lo=330,hi=1500)
        for Pi in (0.1,0.2,0.35):
            Uc=ucrit(Pi=Pi,n=n,kappa=kappa,mode='outlet',lo=330,hi=1500)
            g,d=g_value(625,Uc,Pi=Pi,n=n,kappa=kappa,mode='outlet',full=True); p0=d['p0']
            Ua=ucrit(n=n,kappa=kappa,mode='const',plevel=p0,lo=330,hi=1500)
            r=dict(kind='decomp',kappa=kappa,dTc=kappa*dTad,n=n,Pi=Pi,U0=U0,p0=p0,UC=Uc,UA=Ua,scale_check=U0*p0**n)
            res.append(r); print(r,flush=True)
            json.dump(res,open('coolant_pressure.json','w'),indent=1)
