import time, json, pbr3
from pbr3 import *
yA0=1.65/Pref
for m in ('delasa','cc'):
    pbr3.MODEL=m
    ref=json.load(open(f'ref_slow_{m}.json'))
    t=time.time(); U0=ucrit(yA0,hi=1500); print(m,"U0 fast",U0,"slow",ref[0]['U0'],"diff",U0-ref[0]['U0'],f"{time.time()-t:.0f}s",flush=True)
    r=ref[1]; t=time.time()
    Uc=ucrit(yA0,Pi=r['Pi'],n=r['n'],mode='outlet',hi=1500)
    print("  UC fast",Uc,"slow",r['UC'],"diff",Uc-r['UC'],f"{time.time()-t:.0f}s",flush=True)
