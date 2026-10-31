import json
from pbr import *
out=[]
for n in (2.0,1.0):
    U0=ucrit(n=n)
    for Pi in (0.1,0.2,0.35):
        Uc=ucrit(Pi=Pi,n=n,mode='outlet')
        g,d=g_value(625,Uc,Pi=Pi,n=n,mode='outlet',full=True); p0=d['p0']
        Ua=ucrit(n=n,mode='const',plevel=p0)
        out.append(dict(kappa=0.0,n=n,Pi=Pi,U0=U0,p0=p0,UC=Uc,UA=Ua)); print(out[-1],flush=True)
        json.dump(out,open('base_n12.json','w'),indent=1)
