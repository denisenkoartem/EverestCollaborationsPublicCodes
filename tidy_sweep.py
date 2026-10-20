#!/usr/bin/env python3
"""
tidy_sweep.py -- sort a sweep CSV into a stable order and build pivot tables.

    python3 tidy_sweep.py sweep.csv
    python3 tidy_sweep.py sweep.csv --etas 1,4,8,12
    python3 tidy_sweep.py sweep.csv --limit 20      only |xx|,|yy| <= 20

Writes:
    <name>_sorted.csv    same rows, ordered by xx_pct, yy_pct, eta
    <name>_pivot.csv     one row per (xx,yy), one column per eta
"""
import argparse, csv, os

ap = argparse.ArgumentParser()
ap.add_argument('csv')
ap.add_argument('--etas', default=None, help='comma list, e.g. 1,4,8,12')
ap.add_argument('--value', default='Ep_per_quark')
ap.add_argument('--limit', type=int, default=None,
                help='keep only |xx_pct| and |yy_pct| <= this')
a = ap.parse_args()

rows = list(csv.DictReader(open(a.csv)))
fields = list(rows[0].keys())
for r in rows:
    r['_x'], r['_y'], r['_e'] = int(r['xx_pct']), int(r['yy_pct']), int(r['eta'])

if a.limit is not None:
    rows = [r for r in rows if abs(r['_x']) <= a.limit and abs(r['_y']) <= a.limit]

rows.sort(key=lambda r: (r['_x'], r['_y'], r['_e']))
stem = os.path.splitext(os.path.basename(a.csv))[0]

with open(f'{stem}_sorted.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
    w.writeheader()
    w.writerows(rows)

etas = ([int(x) for x in a.etas.split(',')] if a.etas
        else sorted({r['_e'] for r in rows}))
xs = sorted({r['_x'] for r in rows})
ys = sorted({r['_y'] for r in rows})
g = {(r['_x'], r['_y'], r['_e']): float(r[a.value]) for r in rows}

with open(f'{stem}_pivot.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['xx_pct', 'yy_pct'] + [f'eta{e}' for e in etas])
    for x in xs:
        for y in ys:
            w.writerow([x, y] + [f'{g[(x, y, e)]:.6f}' if (x, y, e) in g else ''
                                 for e in etas])

print(f'{a.value},  rows = (xx%, yy%),  columns = eta')
print('  xx%   yy% ' + ''.join(f'{"eta"+str(e):>11}' for e in etas) + '   note')
for x in xs:
    for y in ys:
        line = f'  {x:+4d}  {y:+4d} ' + ''.join(
            f'{g[(x, y, e)]:11.4f}' if (x, y, e) in g else f'{"":>11}' for e in etas)
        note = ''
        if (x, y, 1) in g and (x, y, 4) in g and g[(x, y, 4)] >= g[(x, y, 1)]:
            note = '  <-- eta4 >= eta1'
        print(line + note)

print(f'\nwrote {stem}_sorted.csv and {stem}_pivot.csv')
