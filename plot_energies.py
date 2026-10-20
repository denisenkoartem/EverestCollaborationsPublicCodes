import os
_HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

UP = _HERE
OUT = _HERE

# columns: N, r1, r2, r3, T, U, V, Ep
lo = np.genfromtxt(f'{UP}/data_paper_B_96_84.dat')   # B = 96.84
hi = np.genfromtxt(f'{UP}/data_paper_B_107_6.dat')   # B = 107.6

N_lo, N_hi = lo[:, 0], hi[:, 0]

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif'],
    'mathtext.fontset': 'dejavuserif',
    'axes.linewidth': 1.1,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'xtick.minor.size': 2.8,
    'ytick.minor.size': 2.8,
    'font.size': 11,
})

# panel: (column index, y-label, colour, marker, linestyle for B=107.6)
panels = [
    (4, 'Kinetic energy/quark (MeV)',   'black',            'o', '-'),
    (5, 'Potential energy/quark (MeV)', 'forestgreen',      's', ':'),
    (6, 'Volume energy/quark (MeV)',    'royalblue',        '^', '--'),
    (7, 'Total energy/quark (MeV)',     'crimson',          'p', '-.'),
]

fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.25))
tags = ['(a)', '(b)', '(c)', '(d)']

for ax, (col, ylab, colour, mark, ls), tag in zip(axes.ravel(), panels, tags):
    # B = 107.6  -> filled markers, panel's signature linestyle
    ax.plot(N_hi, hi[:, col], ls=ls, color=colour, marker=mark,
            ms=7, mfc=colour, mec=colour, lw=1.2,
            label=r'$B = (107.6)^4$ MeV$^4$')
    # B = 96.84  -> open markers, solid thin line
    ax.plot(N_lo, lo[:, col], ls='-', color=colour, marker=mark,
            ms=7, mfc='white', mec=colour, mew=1.4, lw=1.0, alpha=0.95,
            label=r'$B = (96.84)^4$ MeV$^4$')

    ax.set_xlabel('Number of Quarks')
    ax.set_ylabel(ylab)
    ax.set_xlim(0, 75)
    ax.set_xticks([20, 40, 60])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.text(0.5, 0.94, tag, transform=ax.transAxes,
            ha='center', va='top', fontsize=12)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9,
              handlelength=2.6, borderpad=0.4)

fig.suptitle(r'$\alpha_s = 0.3114$', fontsize=13, y=0.985)
fig.tight_layout(pad=1.2, rect=[0, 0, 1, 0.965])
fig.savefig(f'{OUT}/energy_per_quark_two_bag_constants.pdf', bbox_inches='tight')
fig.savefig(f'{OUT}/energy_per_quark_two_bag_constants.png', dpi=200, bbox_inches='tight')

# ---- numbers worth reporting alongside the figure ----
print(f'{"quantity":<10} {"B=107.6 range":>26} {"B=96.84 range":>26} {"mean ratio":>11}')
for col, name in [(4, 'T'), (5, 'U'), (6, 'V'), (7, 'Ep')]:
    a, b = lo[:, col], hi[:, col]
    print(f'{name:<10} {b.min():11.4f} to {b.max():11.4f} '
          f'{a.min():11.4f} to {a.max():11.4f} {np.mean(a / b):11.4f}')
print()
print('B ratio 96.84/107.6 =', 96.84 / 107.6)
print('saved PDF and PNG to', OUT)
