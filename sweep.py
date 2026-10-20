#!/usr/bin/env python3
"""
sweep.py -- scan alfas (xx) and the bag constant root (yy) over +/- 30%.

Default grid is 0, +/-10, +/-20, +/-30 percent on each axis:

    7 xx values  x  7 yy values  =  49 parameter sets
    each parameter set is run for every eta you ask for (default: all 14)

Set the range with --span and the increment with --step.

Usage
-----
    python3 sweep.py                  all 14 etas x 49 sets  (686 runs)
    python3 sweep.py --etas 2         one eta  x 49 sets      (49 runs)
    python3 sweep.py --etas 1,2,5     a subset
    python3 sweep.py --jobs 8         use 8 worker processes
    python3 sweep.py --span 20        0, +/-10, +/-20  (5 values per axis)
    python3 sweep.py --step 15        0, +/-15, +/-30  (5 values per axis)
    python3 sweep.py --root           Brent root-finding for the boundaries
    python3 sweep.py --out my.csv     choose the output file

Output file defaults to sweep_results_pm<span>.csv, so a +/-30% scan is not
mixed into a +/-20% file by accident. Results append as they finish; re-running
skips combinations already present, so an interrupted sweep resumes.

NOTE: percentages are relative to the xx and yy in configs.py. If you edit
those base values, start a new output file -- a row labelled "0%" means "the
base at the time it ran", and mixing bases in one file is silent.

Timing: roughly 15-40 s per run on one core. The default 686 runs is ~4-7 h
serial, well under an hour on 8 cores. Progress is checkpointed every row.
"""

import argparse, csv, os, sys, time
from dataclasses import replace
from multiprocessing import Pool

from tf_core import Model
from configs import CONFIGS

FIELDS = ['xx_pct', 'yy_pct', 'xx', 'yy', 'eta', 'N',
          'x1', 'x2', 'x3', 'R', 'm2',
          'T1', 'V1', 'Volumeenergy', 'Energywithoutmass', 'Energy',
          'r1', 'r2', 'r3', 'T_per_quark', 'U_per_quark', 'V_per_quark', 'Ep_per_quark',
          'seconds']


def percent_steps(span, step):
    """Percentages from -span to +span in increments of step, always incl. 0.

    span=30, step=10  ->  [-30, -20, -10, 0, 10, 20, 30]
    """
    span, step = int(round(span)), int(round(step))
    if step <= 0:
        return [0]
    out = list(range(-span, span + 1, step))
    if 0 not in out:
        out = sorted(out + [0])
    return out


def one_run(job):
    eta_key, xx_pct, yy_pct, solver = job
    base = CONFIGS[eta_key]
    xx = base.xx * (1 + xx_pct / 100.0)
    yy = base.yy * (1 + yy_pct / 100.0)
    cfg = replace(base, xx=xx, yy=yy)

    t0 = time.time()
    try:
        res = Model(cfg, solver=solver).run()
    except Exception as e:
        sys.stderr.write(f'  FAILED eta={eta_key} xx{xx_pct:+d}% yy{yy_pct:+d}%: {e}\n')
        return None
    dt = time.time() - t0

    N = cfg.N
    row = dict(xx_pct=xx_pct, yy_pct=yy_pct, xx=xx, yy=yy, eta=eta_key, N=N,
               seconds=round(dt, 1))
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
    return row


def already_done(path):
    if not os.path.exists(path):
        return set()
    done = set()
    with open(path) as f:
        for r in csv.DictReader(f):
            done.add((int(r['eta']), int(r['xx_pct']), int(r['yy_pct'])))
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--etas', default='all')
    ap.add_argument('--span', type=float, default=30.0,
                    help='percentage span either side of the base (default 30)')
    ap.add_argument('--step', type=int, default=10,
                    help='percentage increment between grid points (default 10)')
    ap.add_argument('--jobs', type=int, default=max(1, (os.cpu_count() or 1)))
    ap.add_argument('--root', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    out = a.out or f'sweep_results_pm{int(a.span)}.csv'
    etas = sorted(CONFIGS) if a.etas == 'all' else [int(x) for x in a.etas.split(',')]
    steps = percent_steps(a.span, a.step)
    solver = 'root' if a.root else 'minimize'

    base = CONFIGS[etas[0]]
    done = already_done(out)
    jobs = [(e, xp, yp, solver) for e in etas for xp in steps for yp in steps
            if (e, xp, yp) not in done]

    total = len(etas) * len(steps) ** 2
    print(f'output file   : {out}')
    print(f'grid (%)      : {steps}')
    print(f'span / step   : +/- {a.span:g}%  in steps of {a.step}%'
          f'   ({len(steps)} values per axis)')
    print(f'etas          : {etas}')
    print(f'base xx, yy   : {base.xx}, {base.yy}   (from configs.py)')
    print(f'xx values     : {[round(base.xx*(1+s/100), 6) for s in steps]}')
    print(f'yy values     : {[round(base.yy*(1+s/100), 4) for s in steps]}')
    print(f'combinations  : {len(steps)**2} parameter sets x {len(etas)} eta = {total} runs')
    print(f'already done  : {len(done)}   to do now: {len(jobs)}')
    print(f'workers       : {a.jobs}\n')
    if not jobs:
        print('nothing to do')
        return

    new = not os.path.exists(out)
    fh = open(out, 'a', newline='')
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
        fh.flush()

    t0 = time.time()
    n = 0
    with Pool(a.jobs) as pool:
        for row in pool.imap_unordered(one_run, jobs):
            n += 1
            if row is None:
                continue
            w.writerow(row)
            fh.flush()
            el = time.time() - t0
            left = el / n * (len(jobs) - n)
            print(f'[{n:4d}/{len(jobs)}] eta={row["eta"]:2d} '
                  f'xx{row["xx_pct"]:+3d}% yy{row["yy_pct"]:+3d}%  '
                  f'Ep/q={row["Ep_per_quark"]:9.4f}  r1={row["r1"]:8.4f}  '
                  f'({row["seconds"]:.0f}s, ~{left/60:.0f} min left)')
    fh.close()
    print(f'\ndone in {(time.time()-t0)/60:.1f} min -> {out}')


if __name__ == '__main__':
    main()
