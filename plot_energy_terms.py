#!/usr/bin/env python3
"""
plot_energy_terms.py -- the individual energy contributions (kinetic, potential,
volume) alongside the total, for every parameter sweep.

    energy_terms_vs_B.pdf        vs B^(1/4) at baseline alpha_s   (extends Fig. 5)
    energy_terms_vs_alphas.pdf   vs alpha_s at baseline B         (extends Fig. 6)
    energy_terms_vs_mass.pdf     vs the quark masses, three sweep modes

All per quark, without rest mass:   E_p = T + U + V_vol
Curves are eta = 1, 4, 8, 12.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

PM30 = os.path.join(_HERE, 'alphas_B_sweep_results.csv')
MASS = os.path.join(_HERE, 'mass_sweep_results.csv')
OUT  = _HERE
SEL  = [1, 4, 8, 12]
PCTS = [-30, -20, -10, 0, 10, 20, 30]

# ---------------------------------------------------------------- style
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'cm', 'axes.linewidth': 1.0,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 6, 'ytick.major.size': 6,
    'xtick.minor.size': 3.2, 'ytick.minor.size': 3.2,
    'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'font.size': 11.5, 'legend.fontsize': 9.5,
    'legend.framealpha': 1.0, 'legend.edgecolor': 'black',
    'legend.fancybox': False,
})

# colour + line style + marker shape + fill, so no single channel carries meaning
STYLE = {
    1:  dict(color='black',       ls='-',  marker='o', fill=True),
    4:  dict(color='crimson',     ls='--', marker='p', fill=False),
    8:  dict(color='royalblue',   ls='-.', marker='^', fill=True),
    12: dict(color='forestgreen', ls=':',  marker='s', fill=False),
}

TERMS = [('T_per_quark',  'Kinetic energy/quark (MeV)'),
         ('U_per_quark',  'Potential energy/quark (MeV)'),
         ('V_per_quark',  'Volume energy/quark (MeV)'),
         ('Ep_per_quark', 'Total energy/quark (MeV)')]


TAGPOS = {'ul': (0.055, 0.945, 'left', 'top'),
          'll': (0.055, 0.055, 'left', 'bottom'),
          'ur': (0.945, 0.945, 'right', 'top'),
          'lr': (0.945, 0.055, 'right', 'bottom')}


def draw(ax, xs, series, xlabel, ylabel, tag, xmaj, legend=None, tagpos='ul'):
    for e in SEL:
        s = STYLE[e]
        ax.plot(xs, series[e], ls=s['ls'], color=s['color'], lw=1.4,
                marker=s['marker'], ms=6.5, mew=1.5,
                mfc=(s['color'] if s['fill'] else 'white'), mec=s['color'],
                label=rf'$\eta = {e}$  ($N = {5*e}$)')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.xaxis.set_major_locator(MultipleLocator(xmaj))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='both', top=True, right=True, direction='in')
    tx, ty, tha, tva = TAGPOS[tagpos]
    ax.text(tx, ty, f'({tag})', transform=ax.transAxes,
            ha=tha, va=tva, fontsize=12.5)
    if legend:
        ax.legend(loc=legend, handlelength=2.6, borderpad=0.5, labelspacing=0.35)


# ================================================================ alpha_s / B
rows = list(csv.DictReader(open(PM30)))
G = {(int(r['xx_pct']), int(r['yy_pct']), int(r['eta'])): r for r in rows}
XXV = {p: float(G[(p, 0, 1)]['xx']) for p in PCTS}
YYV = {p: float(G[(0, p, 1)]['yy']) for p in PCTS}


def sweep_fig(vary, fname, suptitle, xlabel, xvals, xmaj):
    """vary='yy' -> x is the bag constant; vary='xx' -> x is alpha_s."""
    fig, axs = plt.subplots(2, 2, figsize=(9.6, 7.6))
    tags = ['a', 'b', 'c', 'd']
    # rising curves -> legend upper left; falling -> upper right
    legs = ['upper left', None, None, None]
    tagp = ['lr', 'll', 'ul', 'ul']
    for ax, (col, ylab), tag, lg, tp in zip(axs.ravel(), TERMS, tags, legs, tagp):
        ser = {}
        for e in SEL:
            ser[e] = [float(G[((0, p, e) if vary == 'yy' else (p, 0, e))][col])
                      for p in PCTS]
        draw(ax, xvals, ser, xlabel, ylab, tag, xmaj, legend=lg, tagpos=tp)
    fig.suptitle(suptitle, fontsize=12, y=0.985)
    fig.tight_layout(pad=1.2, rect=[0, 0, 1, 0.962])
    fig.savefig(f'{OUT}/{fname}.pdf', bbox_inches='tight')
    fig.savefig(f'{OUT}/{fname}.png', dpi=210, bbox_inches='tight')
    plt.close(fig)
    print('wrote', fname)


sweep_fig('yy', 'energy_terms_vs_B',
          r'$\pm30\%$ variation of $B^{1/4}$ at baseline $\alpha_s = 0.346$',
          r'$B^{1/4}$ (MeV)', [YYV[p] for p in PCTS], 20)

sweep_fig('xx', 'energy_terms_vs_alphas',
          r'$\pm30\%$ variation of $\alpha_s$ at baseline $B^{1/4} = 107.6$ MeV',
          r'$\alpha_s$', [XXV[p] for p in PCTS], 0.05)

# ================================================================ quark masses
mrows = list(csv.DictReader(open(MASS)))
D = {}
for r in mrows:
    p, e = int(r['pct']), int(r['eta'])
    if p == 0:
        for m in ('both', 'light', 'heavy'):
            D[(m, p, e)] = r
    else:
        D[(r['mode'], p, e)] = r

MODES = [('both',  'Change in quark mass (%)',          PCTS, 10),
         ('light', r'$m_{\mathrm{light}}$ (MeV)',      None, 50),
         ('heavy', r'$m_{\mathrm{heavy}}$ (MeV)',      None, 200)]

fig, axs = plt.subplots(3, 3, figsize=(13.0, 10.2))
tags = [['a', 'b', 'c'], ['d', 'e', 'f'], ['g', 'h', 'i']]
for i, (col, ylab) in enumerate(TERMS[:3]):
    for j, (mode, xlab, xv, xmaj) in enumerate(MODES):
        ax = axs[i, j]
        if xv is None:
            key = 'm_light' if mode == 'light' else 'm_heavy'
            xs = [float(D[(mode, p, SEL[0])][key]) for p in PCTS]
        else:
            xs = xv
        ser = {e: [float(D[(mode, p, e)][col]) for p in PCTS] for e in SEL}
        draw(ax, xs, ser, xlab, ylab, tags[i][j], xmaj,
             legend=('upper right' if (i == 0 and j == 0) else None),
             tagpos='ll')
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo - 0.13 * (hi - lo), hi)
    axs[i, 0].set_ylabel(ylab)

for j, (mode, _, _, _) in enumerate(MODES):
    title = {'both': 'both masses scaled together',
             'light': 'light quark mass only',
             'heavy': 'heavy quark mass only'}[mode]
    axs[0, j].set_title(title, fontsize=11.5, pad=8)

fig.suptitle(r'$\pm30\%$ variation of the quark masses at '
             r'$\alpha_s = 0.346$ and $B^{1/4} = 107.6$ MeV', fontsize=12, y=0.992)
fig.tight_layout(pad=1.1, rect=[0, 0, 1, 0.972])
fig.savefig(f'{OUT}/energy_terms_vs_mass.pdf', bbox_inches='tight')
fig.savefig(f'{OUT}/energy_terms_vs_mass.png', dpi=190, bbox_inches='tight')
plt.close(fig)
print('wrote energy_terms_vs_mass')

# ---------------------------------------------------------------- ranges table
print()
print('Range of each term over the full sweep (eta = 1):')
print('%-26s %10s %10s %10s %10s' % ('sweep', 'T', 'U', 'V_vol', 'E_p'))
for lbl, getter in [
        (r'B^(1/4)  75.3->139.9', lambda c: (float(G[(0, -30, 1)][c]), float(G[(0, 30, 1)][c]))),
        ('alpha_s  0.242->0.450', lambda c: (float(G[(-30, 0, 1)][c]), float(G[(30, 0, 1)][c]))),
        ('both masses -30->+30%', lambda c: (float(D[('both', -30, 1)][c]), float(D[('both', 30, 1)][c]))),
        ('m_light     -30->+30%', lambda c: (float(D[('light', -30, 1)][c]), float(D[('light', 30, 1)][c]))),
        ('m_heavy     -30->+30%', lambda c: (float(D[('heavy', -30, 1)][c]), float(D[('heavy', 30, 1)][c])))]:
    cells = []
    for col, _ in TERMS:
        a, b = getter(col)
        cells.append(f'{100*(b-a)/abs(a):+9.1f}%')
    print('%-26s %s' % (lbl, ' '.join(cells)))
