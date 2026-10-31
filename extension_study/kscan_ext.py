import json
from pbr import *
out=[]
for kappa in (0.0,0.0025,0.005,0.01,0.02,0.03,0.04,0.05,0.075,0.1):
    U=ucrit(kappa=kappa,lo=330,hi=2000)
    out.append(dict(kappa=kappa,Ucrit=U)); print(kappa,U,flush=True)
json.dump(out,open('kappa_scan_ext.json','w'),indent=1)
