#!/usr/bin/env python3
"""
plot_mass_sweep.py -- figures in the style of Figures 5 and 6, for the
+/- 30% quark-mass variation.

Figure A  mass_sweep_energy    total energy/quark without rest mass
Figure B  mass_sweep_mass      total mass per quark (rest mass included)

Each has three panels: overall mass scale, light quark only, heavy quark only.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

ap = argparse.ArgumentParser()
ap.add_argument('csv', nargs='?', default='mass_sweep_results.csv')
ap.add_argument('--outdir', default=_HERE)
a = ap.parse_args()

rows = list(csv.DictReader(open(a.csv)))
for r in rows:
    r['_e'], r['_p'] = int(r['eta']), int(r['pct'])
ES = sorted({r['_e'] for r in rows})
PCTS = sorted({r['_p'] for r in rows})
MODES = ['both', 'light', 'heavy']

# at 0% every mode coincides; that row is stored once under 'both'
D = {}
for r in rows:
    if r['_p'] == 0:
        for m in MODES:
            D[(m, r['_p'], r['_e'])] = r
    else:
        D[(r['mode'], r['_p'], r['_e'])] = r


def val(mode, p, e, col):
    return float(D[(mode, p, e)][col])


plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'dejavuserif', 'axes.linewidth': 1.1,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 5, 'ytick.major.size': 5,
    'xtick.minor.size': 2.8, 'ytick.minor.size': 2.8, 'font.size': 10.5,
})

cols = {1: 'black', 4: 'crimson', 8: 'royalblue', 12: 'forestgreen'}
marks = {1: 'o', 4: 'p', 8: '^', 12: 's'}
XLAB = {'both': r'quark mass scale (% change)',
        'light': r'$m_{\rm light}$ (MeV)',
        'heavy': r'$m_{\rm heavy}$ (MeV)'}
TITLE = {'both': r'(a)  both masses scaled together',
         'light': r'(b)  light quark mass only',
         'heavy': r'(c)  heavy quark mass only'}


def xaxis(mode):
    if mode == 'both':
        return PCTS
    col = 'm_light' if mode == 'light' else 'm_heavy'
    return [val(mode, p, ES[0], col) for p in PCTS]


def make(col, ylab, fname, suptitle):
    fig, axs = plt.subplots(1, 3, figsize=(14.2, 4.6))
    for axx, mode in zip(axs, MODES):
        xs = xaxis(mode)
        for e in ES:
            axx.plot(xs, [val(mode, p, e, col) for p in PCTS], '-',
                     marker=marks[e], ms=6, color=cols[e], lw=1.3,
                     label=rf'$\eta={e}$  ($N={5*e}$)')
        axx.set_xlabel(XLAB[mode])
        axx.set_ylabel(ylab)
        axx.set_title(TITLE[mode], fontsize=11)
        axx.xaxis.set_minor_locator(AutoMinorLocator(2))
        axx.yaxis.set_minor_locator(AutoMinorLocator(2))
    axs[0].legend(fontsize=8.5, loc='best')
    fig.suptitle(suptitle, fontsize=12, y=0.99)
    fig.tight_layout(pad=1.2, rect=[0, 0, 1, 0.95])
    fig.savefig(f'{a.outdir}/{fname}.pdf', bbox_inches='tight')
    fig.savefig(f'{a.outdir}/{fname}.png', dpi=200, bbox_inches='tight')
    plt.close(fig)


make('Ep_per_quark', 'Total energy/quark (MeV)', 'mass_sweep_energy',
     r'$\pm30\%$ variation of the quark masses, '
     r'at $\alpha_s = 0.346$ and $B^{1/4} = 107.6$ MeV'
     '\n(total energy without rest mass per quark)')

make('M_per_quark', 'Total mass/quark (MeV)', 'mass_sweep_mass',
     r'$\pm30\%$ variation of the quark masses, '
     r'at $\alpha_s = 0.346$ and $B^{1/4} = 107.6$ MeV'
     '\n(total mass per quark, rest mass included)')

# ---- how much does each quantity move? ----
print(f'{"mode":<7}{"quantity":<14}' + ''.join(f'{"eta"+str(e):>12}' for e in ES))
for mode in MODES:
    for col, lab in [('Ep_per_quark', 'Ep/quark'), ('M_per_quark', 'M/quark')]:
        line = f'{mode:<7}{lab:<14}'
        for e in ES:
            lo, hi = val(mode, PCTS[0], e, col), val(mode, PCTS[-1], e, col)
            line += f'{100*(hi/lo-1):11.1f}%'
        print(line)
print('\nsaved mass_sweep_energy and mass_sweep_mass (pdf + png)')
