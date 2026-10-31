import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.labelsize':9,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'font.family':'DejaVu Sans'})
rows=json.load(open('all_rows.json'))
ks=json.load(open('kappa_scan_ext.json'))
dTad=809.329
col={0.1:'#1f77b4',0.2:'#d62728',0.35:'#2ca02c'}

# ---------------- Figure 13: co-current coolant -----------------
fig,ax=plt.subplots(1,2,figsize=(7.2,3.1))
x=[k['kappa']*dTad for k in ks]; y=[k['Ucrit'] for k in ks]
ax[0].plot(x,y,'o-',color='#1f4e79',ms=4)
ax[0].axhline(ks[0]['Ucrit'],color='grey',ls=':',lw=0.8)
ax[0].set_xlabel('Coolant temperature rise at full conversion, $\\kappa\\,\\Delta T_{ad}$ (K)')
ax[0].set_ylabel('$U_{crit}$ (W m$^{-2}$ K$^{-1}$)')
ax[0].set_title('(a) boundary vs coolant warming ($T_0=T_{c,in}=625$ K)',fontsize=8.5)
ax[0].text(30,ks[0]['Ucrit']-12,'constant-$T_c$ value',fontsize=7.5,color='grey')
# panel b: level / |decline| ratio vs Pi
cc={0.0:'#1f77b4',0.02:'#ff7f0e',0.05:'#2ca02c'}
for kap in (0.0,0.02,0.05):
    for n,ls,m in ((2.0,'-','o'),(1.0,'--','s')):
        sel=sorted([r for r in rows if abs(r['kappa']-kap)<1e-9 and r['n']==n],key=lambda r:r['Pi'])
        ax[1].plot([r['Pi'] for r in sel],[-r['level']/r['decl'] for r in sel],ls,marker=m,color=cc[kap],ms=4,
                   label=f'{kap*dTad:.0f} K, $n$={n:g}')
ax[1].axhline(1,color='k',lw=0.6)
ax[1].set_ylim(0,10)
ax[1].set_xlabel('$\\Pi$ (outlet-fixed)'); ax[1].set_ylabel('level / |decline|')
ax[1].set_title('(b) level-to-decline ratio',fontsize=8.5)
ax[1].legend(title='coolant rise, rate order',ncol=2,frameon=False,fontsize=6.5,title_fontsize=6.5,loc='lower right')
ax[1].text(0.1,1.15,'level = decline',fontsize=6.5)
plt.tight_layout(); plt.savefig('fig14_coolant.png',dpi=300); plt.close()

# ---------------- Figure 14: pressure-order continuum -----------------
fig,ax=plt.subplots(1,3,figsize=(7.2,2.8))
for Pi,c in col.items():
    sel=sorted([r for r in rows if r['kappa']==0.0 and r['Pi']==Pi],key=lambda r:r['n'])
    n=[r['n'] for r in sel]
    ax[0].plot(n,[r['level'] for r in sel],'o-',color=c,ms=3.5,label=f'$\\Pi$={Pi}')
    ax[1].plot(n,[-r['decl'] for r in sel],'o-',color=c,ms=3.5)
    ax[2].plot(n,[r['peff_frac'] for r in sel],'o-',color=c,ms=3.5)
    selc=sorted([r for r in rows if r['kappa']==0.02 and r['Pi']==Pi],key=lambda r:r['n'])
    if selc:
        ax[2].plot([r['n'] for r in selc],[r['peff_frac'] for r in selc],'s',mfc='none',color=c,ms=5)
ax[0].set_title('(a) level contribution',fontsize=8.5); ax[0].set_ylabel('W m$^{-2}$ K$^{-1}$')
ax[1].set_title('(b) decline magnitude',fontsize=8.5)
ax[2].set_title('(c) $(p_{eff}-1)/(p_0-1)$',fontsize=8.5)
ax[2].set_ylim(0.7,1.0)
for a in ax: a.set_xlabel('total-pressure order $n$')
ax[0].legend(frameon=False,fontsize=7)
ax[2].text(1.02,0.715,'open squares: co-current\ncoolant, 16 K rise',fontsize=6.5)
plt.tight_layout(); plt.savefig('fig13_order.png',dpi=300); plt.close()
print('figs done')
