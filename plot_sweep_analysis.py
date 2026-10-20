#!/usr/bin/env python3
"""
plot_sweep_analysis.py -- visual summary of the 7x7 (alpha_s, B) sweep.

Figure 1  sweep_eta_comparison   the eta = 1,4,8,12 question Wilcox raised
Figure 2  sweep_energy_vs_N      T, U, V, Ep per quark vs number of quarks

    python3 plot_sweep_analysis.py 20260728_pm30run.csv
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap, BoundaryNorm
from matplotlib.ticker import AutoMinorLocator

ap = argparse.ArgumentParser()
ap.add_argument('csv', nargs='?',
                default=os.path.join(_HERE, 'alphas_B_sweep_results.csv'))
ap.add_argument('--outdir', default=_HERE)
a = ap.parse_args()

rows = list(csv.DictReader(open(a.csv)))
for r in rows:
    r['_x'], r['_y'], r['_e'] = int(r['xx_pct']), int(r['yy_pct']), int(r['eta'])
XS = sorted({r['_x'] for r in rows})
YS = sorted({r['_y'] for r in rows})
ES = sorted({r['_e'] for r in rows})
V = {(r['_x'], r['_y'], r['_e']): r for r in rows}
XX0 = {r['_x']: float(r['xx']) for r in rows}
YY0 = {r['_y']: float(r['yy']) for r in rows}


def val(x, y, e, col='Ep_per_quark'):
    return float(V[(x, y, e)][col])


Ns = [5 * e for e in ES]


plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'dejavuserif', 'axes.linewidth': 1.1,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 5, 'ytick.major.size': 5, 'font.size': 10.5,
})

# ===========================================================================
# Figure 1 -- the eta = 1, 4, 8, 12 comparison
# ===========================================================================
SEL = [e for e in (1, 4, 8, 12) if e in ES]
cols = {1: 'black', 4: 'crimson', 8: 'royalblue', 12: 'forestgreen'}
marks = {1: 'o', 4: 'p', 8: '^', 12: 's'}

fig = plt.figure(figsize=(14.2, 4.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.12], wspace=0.28)
axA, axB, axC = (fig.add_subplot(gs[0]), fig.add_subplot(gs[1]),
                 fig.add_subplot(gs[2]))

# (a) vs bag constant, at baseline alpha_s
for e in SEL:
    axA.plot([YY0[y] for y in YS], [val(0, y, e) for y in YS],
             '-', marker=marks[e], ms=6, color=cols[e], lw=1.3,
             label=rf'$\eta={e}$  ($N={5*e}$)')
axA.set_xlabel(r'$B^{1/4}$ (MeV)')
axA.set_ylabel('Total energy/quark (MeV)')
axA.set_title(r'(a)  varying $B$ at baseline $\alpha_s$', fontsize=11)
axA.legend(fontsize=8.5, loc='upper left')
axA.xaxis.set_minor_locator(AutoMinorLocator(2))
axA.yaxis.set_minor_locator(AutoMinorLocator(2))

# (b) vs alpha_s, at baseline bag constant
for e in SEL:
    axB.plot([XX0[x] for x in XS], [val(x, 0, e) for x in XS],
             '-', marker=marks[e], ms=6, color=cols[e], lw=1.3,
             label=rf'$\eta={e}$')
axB.set_xlabel(r'$\alpha_s$')
axB.set_ylabel('Total energy/quark (MeV)')
axB.set_title(r'(b)  varying $\alpha_s$ at baseline $B$', fontsize=11)
axB.legend(fontsize=8.5, loc='upper left')
axB.xaxis.set_minor_locator(AutoMinorLocator(2))
axB.yaxis.set_minor_locator(AutoMinorLocator(2))

# (c) which eta is the global minimum
# ===========================================================================
# Figure 3 -- energy per quark RELATIVE to eta = 1
#   removes the overall vertical offset so the crossings become visible
# ===========================================================================
def relative_to_eta1(vary, fname, suptitle, legend_fmt):
    keys = YS if vary == 'yy' else XS
    shades = plt.cm.viridis(np.linspace(0.05, 0.88, len(keys)))
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 4.9))
    for axx, (col, ylab) in zip(axs, [('Ep_per_quark', 'Total'),
                                      ('T_per_quark', 'Kinetic')]):
        for k, colr in zip(keys, shades):
            get = (lambda e: val(0, k, e, col)) if vary == 'yy' \
                else (lambda e: val(k, 0, e, col))
            ref = get(1)
            axx.plot(Ns, [get(e) - ref for e in ES], '-', marker='o', ms=4.5,
                     color=colr, lw=1.2, label=legend_fmt(k))
        axx.axhline(0, color='black', lw=1.0, ls='--', zorder=0)
        axx.set_xlabel('Number of Quarks')
        axx.set_ylabel(rf'{ylab} energy/quark $-$ value at $\eta=1$  (MeV)')
        axx.xaxis.set_minor_locator(AutoMinorLocator(2))
        axx.yaxis.set_minor_locator(AutoMinorLocator(2))
    axs[0].set_title('(a)  total energy', fontsize=11)
    axs[1].set_title('(b)  kinetic energy', fontsize=11)
    axs[0].legend(fontsize=7.5, ncol=2, handlelength=1.8, loc='best')
    fig.suptitle(suptitle, fontsize=12, y=0.985)
    fig.tight_layout(pad=1.2, rect=[0, 0, 1, 0.95])
    fig.savefig(f'{a.outdir}/{fname}.pdf', bbox_inches='tight')
    fig.savefig(f'{a.outdir}/{fname}.png', dpi=200, bbox_inches='tight')
    plt.close(fig)


relative_to_eta1('xx', 'sweep_relative_to_eta1_alphas',
                 r'energy per quark relative to $\eta=1$, varying $\alpha_s$'
                 rf'  ($B^{{1/4}} = {YY0[0]:g}$ MeV).'
                 '  below zero = configuration favoured over $\\eta=1$',
                 lambda x: rf'$\alpha_s={XX0[x]:.4f}$')

relative_to_eta1('yy', 'sweep_relative_to_eta1_B',
                 r'energy per quark relative to $\eta=1$, varying $B$'
                 rf'  ($\alpha_s = {XX0[0]:g}$).'
                 '  below zero = configuration favoured over $\\eta=1$',
                 lambda y: rf'$B^{{1/4}}={YY0[y]:.1f}$')

W = np.array([[min(ES, key=lambda e: val(x, y, e)) for x in XS] for y in YS])
uniq = sorted(set(W.ravel()))
cmap = ListedColormap([cols.get(u, 'grey') for u in uniq])
idx = np.vectorize(lambda v: uniq.index(v))(W)
axC.pcolormesh(np.arange(len(XS) + 1), np.arange(len(YS) + 1), idx,
               cmap=cmap, norm=BoundaryNorm(range(len(uniq) + 1), len(uniq)))
for i in range(len(YS)):
    for j in range(len(XS)):
        axC.text(j + .5, i + .5, rf'$\eta$={W[i, j]}', ha='center',
                 va='center', fontsize=8.5, color='white')
axC.set_xticks(np.arange(len(XS)) + .5, [f'{x:+d}' for x in XS])
axC.set_yticks(np.arange(len(YS)) + .5, [f'{y:+d}' for y in YS])
axC.set_xlabel(r'change in $\alpha_s$ (%)')
axC.set_ylabel(r'change in $B^{1/4}$ (%)')
axC.set_title('(c)  configuration with the lowest energy/quark', fontsize=11)
axC.tick_params(top=False, right=False)

fig.tight_layout(pad=1.1)
fig.savefig(f'{a.outdir}/sweep_eta_comparison.pdf', bbox_inches='tight')
fig.savefig(f'{a.outdir}/sweep_eta_comparison.png', dpi=200, bbox_inches='tight')

# ===========================================================================
# Figure 2 -- T, U, V, Ep per quark vs number of quarks
#   2a: curves = bag constant, at baseline alpha_s
#   2b: curves = alpha_s,      at baseline bag constant
# ===========================================================================
panels = [('T_per_quark', 'Kinetic energy/quark (MeV)', 'o'),
          ('U_per_quark', 'Potential energy/quark (MeV)', 's'),
          ('V_per_quark', 'Volume energy/quark (MeV)', '^'),
          ('Ep_per_quark', 'Total energy/quark (MeV)', 'p')]
tags = ['(a)', '(b)', '(c)', '(d)']


def energy_vs_N(vary, fname, suptitle, legend_fmt):
    """vary='yy' -> one curve per bag constant; vary='xx' -> one per alpha_s."""
    keys = YS if vary == 'yy' else XS
    shades = plt.cm.viridis(np.linspace(0.05, 0.88, len(keys)))
    fig, axs = plt.subplots(2, 2, figsize=(9.6, 7.6))
    for axx, (col, ylab, mk), tag in zip(axs.ravel(), panels, tags):
        for k, colr in zip(keys, shades):
            ys = [val(0, k, e, col) if vary == 'yy' else val(k, 0, e, col)
                  for e in ES]
            axx.plot(Ns, ys, '-', marker=mk, ms=5, color=colr, lw=1.1,
                     label=legend_fmt(k))
        axx.set_xlabel('Number of Quarks')
        axx.set_ylabel(ylab)
        axx.xaxis.set_minor_locator(AutoMinorLocator(2))
        axx.yaxis.set_minor_locator(AutoMinorLocator(2))
        axx.text(0.5, 0.95, tag, transform=axx.transAxes, ha='center',
                 va='top', fontsize=12)
    axs[0, 0].legend(fontsize=7, loc='upper right', ncol=2, handlelength=1.8)
    fig.suptitle(suptitle, fontsize=12, y=0.985)
    fig.tight_layout(pad=1.2, rect=[0, 0, 1, 0.965])
    fig.savefig(f'{a.outdir}/{fname}.pdf', bbox_inches='tight')
    fig.savefig(f'{a.outdir}/{fname}.png', dpi=200, bbox_inches='tight')
    plt.close(fig)


energy_vs_N('yy', 'sweep_energy_vs_N',
            r'varying the bag constant at baseline $\alpha_s$'
            rf'  ($\alpha_s = {XX0[0]:g}$)',
            lambda y: rf'$B^{{1/4}}={YY0[y]:.1f}$')

energy_vs_N('xx', 'sweep_energy_vs_N_alphas',
            r'varying the strong coupling constant at baseline $B$'
            rf'  ($B^{{1/4}} = {YY0[0]:g}$ MeV)',
            lambda x: rf'$\alpha_s={XX0[x]:.4f}$')

W = np.array([[min(ES, key=lambda e: val(x, y, e)) for x in XS] for y in YS])
bad = [(x, y) for x in XS for y in YS if val(x, y, 4) >= val(x, y, 1)]
print('eta=4 not favoured over eta=1 in', len(bad), 'of', len(XS) * len(YS),
      'combos:', bad)
import collections
print('global minimum counts:', dict(collections.Counter(W.ravel().tolist())))
print('saved sweep_eta_comparison, sweep_energy_vs_N, '
      'sweep_energy_vs_N_alphas, sweep_relative_to_eta1_* (pdf + png)')
