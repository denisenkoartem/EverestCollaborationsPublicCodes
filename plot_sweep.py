#!/usr/bin/env python3
"""
plot_sweep.py -- four-panel figure from sweep_results.csv.

x-axis  : yy = B^(1/4) in MeV
curves  : one per alpha_s value
panels  : (a) kinetic (b) potential (c) volume (d) total, all per quark

If the sweep covers more than one eta, pass --eta N to pick which one.
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
ap.add_argument('--csv', default=os.path.join(_HERE, 'sweep_results.csv'))
ap.add_argument('--eta', type=int, default=None)
ap.add_argument('--out', default=os.path.join(_HERE, 'sweep_energies'))
a = ap.parse_args()

rows = list(csv.DictReader(open(a.csv)))
etas = sorted({int(r['eta']) for r in rows})
eta = a.eta if a.eta is not None else etas[0]
rows = [r for r in rows if int(r['eta']) == eta]
N = int(rows[0]['N'])

alphas = sorted({float(r['xx']) for r in rows})
bags = sorted({float(r['yy']) for r in rows})
grid = {(float(r['xx']), float(r['yy'])): r for r in rows}

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'dejavuserif', 'axes.linewidth': 1.1,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 5, 'ytick.major.size': 5,
    'xtick.minor.size': 2.8, 'ytick.minor.size': 2.8, 'font.size': 11,
})

panels = [
    ('T_per_quark',  'Kinetic energy/quark (MeV)',   'o'),
    ('U_per_quark',  'Potential energy/quark (MeV)', 's'),
    ('V_per_quark',  'Volume energy/quark (MeV)',    '^'),
    ('Ep_per_quark', 'Total energy/quark (MeV)',     'p'),
]
tags = ['(a)', '(b)', '(c)', '(d)']
colours = plt.cm.viridis(np.linspace(0.05, 0.85, len(alphas)))

fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4))
for ax, (col, ylab, mark), tag in zip(axes.ravel(), panels, tags):
    for al, colr in zip(alphas, colours):
        y = [float(grid[(al, b)][col]) for b in bags]
        ax.plot(bags, y, '-', marker=mark, ms=6, color=colr, lw=1.3,
                mec=colr, mfc=colr, label=rf'$\alpha_s = {al:.4f}$')
    ax.set_xlabel(r'$B^{1/4}$ (MeV)')
    ax.set_ylabel(ylab)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.text(0.5, 0.95, tag, transform=ax.transAxes, ha='center', va='top', fontsize=12)
    ax.axvline(107.6, color='0.6', lw=0.8, ls=':', zorder=0)

axes[0, 0].legend(loc='upper left', fontsize=8, framealpha=0.9, handlelength=2.2)
fig.suptitle(rf'$\eta = {eta}$   ($N = {N}$ quarks)', fontsize=13, y=0.985)
fig.tight_layout(pad=1.2, rect=[0, 0, 1, 0.965])
fig.savefig(a.out + '.pdf', bbox_inches='tight')
fig.savefig(a.out + '.png', dpi=200, bbox_inches='tight')
print('eta =', eta, ' N =', N)
print('alpha_s values :', alphas)
print('B^(1/4) values :', bags)
print('saved', a.out + '.pdf / .png')
