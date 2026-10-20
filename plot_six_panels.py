#!/usr/bin/env python3
"""
plot_six_panels.py -- the six panels from the slide figures, as standalone plots.

    (a) E_p/quark vs overall quark mass scale
    (b) E_p/quark vs m_light
    (c) E_p/quark vs m_heavy
    (d) M/quark   vs overall quark mass scale
    (e) M/quark   vs m_light
    (f) M/quark   vs m_heavy

Each is written as its own PDF and PNG.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

CSV = 'mass_sweep_results.csv'
OUT = os.path.join(_HERE, 'panels')
os.makedirs(OUT, exist_ok=True)

rows = list(csv.DictReader(open(CSV)))
for r in rows:
    r['_e'], r['_p'] = int(r['eta']), int(r['pct'])
ES = sorted({r['_e'] for r in rows})
PCTS = sorted({r['_p'] for r in rows})
MODES = ['both', 'light', 'heavy']

D = {}
for r in rows:
    if r['_p'] == 0:
        for m in MODES:
            D[(m, r['_p'], r['_e'])] = r
    else:
        D[(r['mode'], r['_p'], r['_e'])] = r

val = lambda m, p, e, c: float(D[(m, p, e)][c])

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'axes.linewidth': 1.0,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 6, 'ytick.major.size': 6,
    'xtick.minor.size': 3.2, 'ytick.minor.size': 3.2,
    'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'font.size': 12, 'legend.fontsize': 10.5,
    'legend.framealpha': 1.0, 'legend.edgecolor': 'black',
    'legend.fancybox': False,
})

# colour + line style + marker + fill, so no single channel carries the meaning
STYLE = {
    1:  dict(color='black',       ls='-',   marker='o', fill=True),
    4:  dict(color='crimson',     ls='--',  marker='p', fill=False),
    8:  dict(color='royalblue',   ls='-.',  marker='^', fill=True),
    12: dict(color='forestgreen', ls=':',   marker='s', fill=False),
}

XLAB = {'both':  'Change in quark mass (%)',
        'light': r'$m_{\mathrm{light}}$ (MeV)',
        'heavy': r'$m_{\mathrm{heavy}}$ (MeV)'}
XTICK = {'both': 10, 'light': 50, 'heavy': 200}


def xaxis(mode):
    if mode == 'both':
        return PCTS
    col = 'm_light' if mode == 'light' else 'm_heavy'
    return [val(mode, p, ES[0], col) for p in PCTS]


PANELS = [
    ('a', 'both',  'Ep_per_quark', 'Total energy/quark (MeV)',  'panel_a_energy_scale'),
    ('b', 'light', 'Ep_per_quark', 'Total energy/quark (MeV)',  'panel_b_energy_mlight'),
    ('c', 'heavy', 'Ep_per_quark', 'Total energy/quark (MeV)',  'panel_c_energy_mheavy'),
    ('d', 'both',  'M_per_quark',  'Total mass/quark (MeV)',    'panel_d_mass_scale'),
    ('e', 'light', 'M_per_quark',  'Total mass/quark (MeV)',    'panel_e_mass_mlight'),
    ('f', 'heavy', 'M_per_quark',  'Total mass/quark (MeV)',    'panel_f_mass_mheavy'),
]

# energy panels fall left-to-right -> legend upper right, tag lower left
# mass panels rise    left-to-right -> legend upper left,  tag lower right
LEGLOC = {'Ep_per_quark': 'upper right', 'M_per_quark': 'upper left'}
TAGPOS = {'Ep_per_quark': (0.055, 0.075, 'left', 'bottom'),
          'M_per_quark':  (0.945, 0.075, 'right', 'bottom')}

for tag, mode, col, ylab, fname in PANELS:
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    xs = xaxis(mode)
    for e in ES:
        st = STYLE[e]
        ax.plot(xs, [val(mode, p, e, col) for p in PCTS],
                ls=st['ls'], color=st['color'], lw=1.4,
                marker=st['marker'], ms=7, mew=1.5,
                mfc=(st['color'] if st['fill'] else 'white'),
                mec=st['color'],
                label=rf'$\eta = {e}$  ($N = {5*e}$)')

    ax.set_xlabel(XLAB[mode])
    ax.set_ylabel(ylab)
    ax.xaxis.set_major_locator(MultipleLocator(XTICK[mode]))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='both', top=True, right=True, direction='in')

    tx, ty, tha, tva = TAGPOS[col]
    ax.text(tx, ty, f'({tag})', transform=ax.transAxes,
            ha=tha, va=tva, fontsize=13)
    ax.legend(loc=LEGLOC[col], handlelength=2.6, borderpad=0.5,
              labelspacing=0.35)

    fig.tight_layout(pad=0.6)
    fig.savefig(f'{OUT}/{fname}.pdf', bbox_inches='tight')
    fig.savefig(f'{OUT}/{fname}.png', dpi=220, bbox_inches='tight')
    plt.close(fig)
    print('wrote', fname)

print('\nall six panels in', OUT)
