# derived quantities for the Pi = 0.121 reruns (same definitions as the earlier row files)
import json, math, sys
f = 4184/3600
def rows(model):
    res = json.load(open(f'pi121_{model}.json'))
    out = []
    for r in res:
        if r['kind'] != 'decomp': continue
        U0, UC, UA, p0, n = r['U0'], r['UC'], r['UA'], r['p0'], r['n']
        level = (UA-U0)*f; decl = (UC-UA)*f; comb = (UC-U0)*f
        peff = (UC/U0)**(1/n); phi = (peff-1)/(p0-1); neff = math.log(UA/U0)/math.log(p0)
        out.append(dict(r, U0_W=U0*f, UC_W=UC*f, UA_W=UA*f, scale_check_W=r['scale_check']*f, level=level, decl=decl, comb=comb,
                        ratio=(level/abs(decl) if decl < 0 else float('nan')), peff=peff, phi=phi, scale_dev_W=(UA-r['scale_check'])*f, scale_err_W=(UA-r['scale_check'])*f, neff=neff))
    return out
if __name__ == '__main__':
    for m in sys.argv[1:]:
        R = rows(m); json.dump(R, open(f'pi121_{m}_rows.json', 'w'), indent=1)
        print('==', m)
        for r in R: print(f"kappa={r['kappa']:.4f} n={r['n']:.0f} Pi={r['Pi']} p0={r['p0']:.4f} U0={r['U0_W']:.2f} level={r['level']:+.3f} decl={r['decl']:+.3f} comb={r['comb']:+.3f} ratio={r['ratio']:.2f} phi={r['phi']:.3f} neff={r['neff']:.2f} scale_dev={r['scale_dev_W']:+.5f}")
