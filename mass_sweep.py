#!/usr/bin/env python3
"""
mass_sweep.py -- vary the quark masses by +/- 30%, holding alpha_s and B fixed.

Three modes, run as separate curve sets:

    both    m_light and m_heavy scaled together  (mass ratio aa unchanged)
    light   m_light only                         (aa scales by the factor)
    heavy   m_heavy only                         (aa scales by 1/factor)

The light mass is not just an additive rest-mass term: it sets the length
scale through a1 = 1/zz, so it enters R, c, d and the volume energy as well.
The heavy mass enters dynamically only through aa = m_light/m_heavy.

    M     = E_wom + n1 g1 m_light + (n2 g2 + n3 g3) m_heavy
    R     = (3a)/(8 alpha_s) (3 pi eta/2)^(2/3),      a = 1/m_light
    V_vol = (4 pi/3) R^3 x1^3 B

Usage
-----
    python3 mass_sweep.py                     etas 1,4,8,12, all three modes
    python3 mass_sweep.py --etas 1,4          a subset
    python3 mass_sweep.py --modes light       one mode
    python3 mass_sweep.py --jobs 8
"""

import argparse, csv, os, sys, time
from dataclasses import replace
from multiprocessing import Pool

from tf_core import Model
from configs import CONFIGS

# hold the couplings at the values used for Figures 5 and 6
BASE_XX = 0.346      # alpha_s
BASE_YY = 107.6      # B^(1/4) in MeV

FIELDS = ['mode', 'pct', 'factor', 'm_light', 'm_heavy', 'aa', 'eta', 'N',
          'x1', 'x2', 'x3', 'R', 'm2',
          'T1', 'V1', 'Volumeenergy', 'Energywithoutmass', 'Energy',
          'r1', 'r2', 'r3', 'T_per_quark', 'U_per_quark', 'V_per_quark',
          'Ep_per_quark', 'M_per_quark', 'seconds']


def make_cfg(eta_key, mode, pct):
    """Config with alpha_s and B pinned, and the masses scaled by mode/pct."""
    base = CONFIGS[eta_key]
    f = 1 + pct / 100.0
    fl = f if mode in ('both', 'light') else 1.0
    fh = f if mode in ('both', 'heavy') else 1.0
    return replace(base,
                   xx=BASE_XX,
                   yy=BASE_YY,
                   zz=base.zz * fl,             # light mass -> a1 = 1/zz
                   aa=base.aa * fl / fh,        # alfa2 = m_light / m_heavy
                   m_light=base.m_light * fl,
                   m_heavy=base.m_heavy * fh)


def one_run(job):
    eta_key, mode, pct, solver = job
    cfg = make_cfg(eta_key, mode, pct)
    t0 = time.time()
    try:
        res = Model(cfg, solver=solver).run()
    except Exception as e:
        sys.stderr.write(f'  FAILED eta={eta_key} {mode} {pct:+d}%: {e}\n')
        return None
    dt = time.time() - t0

    N = cfg.N
    row = dict(mode=mode, pct=pct, factor=1 + pct / 100.0,
               m_light=cfg.m_light, m_heavy=cfg.m_heavy, aa=cfg.aa,
               eta=eta_key, N=N, seconds=round(dt, 1))
    for k in ['x1', 'x2', 'x3', 'R', 'm2', 'T1', 'V1',
              'Volumeenergy', 'Energywithoutmass', 'Energy']:
        row[k] = res[k]
    row['r1'] = res['R'] * res['x1'] * 1e3
    row['r2'] = res['R'] * res['x2'] * 1e3
    row['r3'] = res['R'] * res['x3'] * 1e3
    row['T_per_quark'] = res['T1'] / N
    row['U_per_quark'] = res['V1'] / N
    row['V_per_quark'] = res['Volumeenergy'] / N
    row['Ep_per_quark'] = res['Energywithoutmass'] / N
    row['M_per_quark'] = res['Energy'] / N
    return row


def already_done(path):
    if not os.path.exists(path):
        return set()
    with open(path) as f:
        return {(int(r['eta']), r['mode'], int(r['pct']))
                for r in csv.DictReader(f)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--etas', default='1,4,8,12')
    ap.add_argument('--modes', default='both,light,heavy')
    ap.add_argument('--span', type=float, default=30.0)
    ap.add_argument('--step', type=int, default=10)
    ap.add_argument('--jobs', type=int, default=max(1, os.cpu_count() or 1))
    ap.add_argument('--root', action='store_true')
    ap.add_argument('--out', default='mass_sweep_results.csv')
    a = ap.parse_args()

    etas = [int(x) for x in a.etas.split(',')]
    modes = a.modes.split(',')
    span, step = int(a.span), a.step
    pcts = list(range(-span, span + 1, step))
    if 0 not in pcts:
        pcts = sorted(pcts + [0])
    solver = 'root' if a.root else 'minimize'

    done = already_done(a.out)
    jobs = []
    for e in etas:
        for m in modes:
            for p in pcts:
                # at 0% all modes coincide; run it once under 'both'
                if p == 0 and m != modes[0]:
                    continue
                if (e, m, p) not in done:
                    jobs.append((e, m, p, solver))

    print(f'output        : {a.out}')
    print(f'alpha_s, B1/4 : {BASE_XX}, {BASE_YY}  (held fixed)')
    print(f'etas          : {etas}')
    print(f'modes         : {modes}')
    print(f'percentages   : {pcts}')
    print(f'm_light range : {[round(306*(1+p/100), 1) for p in pcts]}')
    print(f'm_heavy range : {[round(1553*(1+p/100), 1) for p in pcts]}')
    print(f'already done  : {len(done)}   to do now: {len(jobs)}')
    print(f'workers       : {a.jobs}\n')
    if not jobs:
        print('nothing to do')
        return

    new = not os.path.exists(a.out)
    fh = open(a.out, 'a', newline='')
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
        fh.flush()

    t0, n = time.time(), 0
    with Pool(a.jobs) as pool:
        for row in pool.imap_unordered(one_run, jobs):
            n += 1
            if row is None:
                continue
            w.writerow(row)
            fh.flush()
            left = (time.time() - t0) / n * (len(jobs) - n)
            print(f'[{n:3d}/{len(jobs)}] eta={row["eta"]:2d} {row["mode"]:<5s} '
                  f'{row["pct"]:+3d}%  Ep/q={row["Ep_per_quark"]:8.4f}  '
                  f'M/q={row["M_per_quark"]:9.3f}  ({row["seconds"]:.0f}s, '
                  f'~{left/60:.0f} min left)')
    fh.close()
    print(f'\ndone in {(time.time()-t0)/60:.1f} min -> {a.out}')


if __name__ == '__main__':
    main()
