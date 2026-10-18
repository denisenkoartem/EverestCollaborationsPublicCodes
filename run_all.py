#!/usr/bin/env python3
"""
run_all.py -- run any or all of the 14 eta calculations.

  python3 run_all.py            run all 14
  python3 run_all.py 3          run eta = 3 only
  python3 run_all.py 1 2 5      run a subset
  python3 run_all.py --root     use Brent root-finding instead of FindMinimum

Writes eta_results_python.dat (raw values) and
data_paper_python.dat (same columns as data_paper.dat).
"""

import sys, time
import numpy as np
from tf_core import Model, P_VOL
from configs import CONFIGS

# final-loop values read out of each notebook's own saved output, for checking
REFERENCE = {
    # eta 1: notebook prints x2bar (inner) and x2 (middle); mapped to x2/x3 here
    1:  dict(x1=0.502821140, x2=0.288060384, x3=0.294051189, m2=12.822,
             T1=116.562550, V1=-41.215531, Volumeenergy=63.387596,
             Energywithoutmass=138.734615),
    2:  dict(x1=0.398164724, x2=0.192907181, x3=0.206159501, m2=17.922,
             T1=225.354997, V1=-69.808330, Volumeenergy=125.896028,
             Energywithoutmass=281.442696),
    3:  dict(x1=0.353869499, x2=0.194000887, x3=0.209978744, m2=27.643,
             T1=343.998557, V1=-87.729432, Volumeenergy=198.854579,
             Energywithoutmass=455.123704),
    4:  dict(x1=0.306772766, x2=0.152885920, x3=0.169327959, m2=28.234,
             T1=405.318470, V1=-114.464817, Volumeenergy=230.321291,
             Energywithoutmass=521.174943),
    5:  dict(x1=0.299818122, x2=0.162548702, x3=0.179602650, m2=39.112,
             T1=569.390298, V1=-123.412397, Volumeenergy=335.952084,
             Energywithoutmass=781.929985),
    6:  dict(x1=0.279074739, x2=0.131592682, x3=0.148373947, m2=37.812,
             T1=664.028327, V1=-146.177120, Volumeenergy=390.146557,
             Energywithoutmass=907.997765),
    7:  dict(x1=0.268744686, x2=0.145106543, x3=0.161743337, m2=48.902,
             T1=793.332339, V1=-154.291141, Volumeenergy=474.219527,
             Energywithoutmass=1113.260725),
    8:  dict(x1=0.245144084, x2=0.121306299, x3=0.137006858, m2=44.412,
             T1=800.052099, V1=-181.137416, Volumeenergy=470.119127,
             Energywithoutmass=1089.033809),
    9:  dict(x1=0.247715929, x2=0.133493874, x3=0.149518535, m2=57.642,
             T1=1015.598502, V1=-182.150100, Volumeenergy=613.918210,
             Energywithoutmass=1447.366612),
    10: dict(x1=0.236463594, x2=0.111190073, x3=0.126669670, m2=52.644,
             T1=1096.448159, V1=-204.896187, Volumeenergy=659.259741,
             Energywithoutmass=1550.811713),
    11: dict(x1=0.231968180, x2=0.124881060, x3=0.140296653, m2=65.751,
             T1=1238.399285, V1=-208.041509, Volumeenergy=753.068205,
             Energywithoutmass=1783.425982),
    12: dict(x1=0.214747471, x2=0.106167108, x3=0.120662790, m2=57.713,
             T1=1193.005591, V1=-236.782814, Volumeenergy=711.066731,
             Energywithoutmass=1667.289508),
    13: dict(x1=0.219641772, x2=0.118177949, x3=0.133042103, m2=73.321,
             T1=1460.415421, V1=-232.334114, Volumeenergy=892.884360,
             Energywithoutmass=2120.965668),
    14: dict(x1=0.211826987, x2=0.099578455, x3=0.113925858, m2=65.413,
             T1=1527.936214, V1=-255.927513, Volumeenergy=928.888598,
             Energywithoutmass=2200.897299),
}

CHECK = ['x1', 'x2', 'x3', 'm2', 'T1', 'V1', 'Volumeenergy', 'Energywithoutmass']


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    solver = 'root' if '--root' in sys.argv else 'minimize'
    which = [int(a) for a in args] if args else sorted(CONFIGS)

    results = {}
    for k in which:
        cfg = CONFIGS[k]
        t0 = time.time()
        res = Model(cfg, solver=solver).run()
        res['eta'] = k
        res['N'] = cfg.N
        results[k] = res
        ref = REFERENCE.get(k, {})
        worst = max((abs(res[q] - ref[q]) / abs(ref[q]) for q in CHECK if q in ref),
                    default=float('nan'))
        print(f'eta={k:2d}  N={cfg.N:3d}  Ewom={res["Energywithoutmass"]:12.6f}  '
              f'x1={res["x1"]:.9f}  worst rel diff vs notebook = {worst:.2e}  '
              f'[{time.time()-t0:.0f}s]')

    # ---------- detailed comparison ----------
    print('\n' + '=' * 78)
    print('COMPARISON WITH NOTEBOOK SAVED OUTPUT  (relative difference)')
    print('=' * 78)
    print(f'{"eta":>4} ' + ' '.join(f'{q[:9]:>10}' for q in CHECK))
    for k in which:
        ref = REFERENCE.get(k)
        if not ref:
            continue
        row = []
        for q in CHECK:
            row.append(f'{abs(results[k][q]-ref[q])/abs(ref[q]):10.1e}' if q in ref else f'{"-":>10}')
        print(f'{k:4d} ' + ' '.join(row))

    # ---------- write output files ----------
    with open('eta_results_python.dat', 'w') as f:
        f.write('# computed by the Python port (tf_core.py)\n')
        f.write('#eta   N      x1              x2              x3              '
                'R               m2              T1              V1              '
                'Volume          Ewom            Energy\n')
        for k in sorted(results):
            r = results[k]
            f.write(f'{k:<6d}{r["N"]:<7d}' + ''.join(
                f'{r[q]:<16.9f}' for q in ['x1', 'x2', 'x3', 'R', 'm2']) +
                ''.join(f'{r[q]:<16.6f}' for q in
                        ['T1', 'V1', 'Volumeenergy', 'Energywithoutmass', 'Energy']) + '\n')

    with open('data_paper_python.dat', 'w') as f:
        f.write('#     r1(*10^-3)\tr2(*10^-3)\tr3(*10^-3)\t     T\t\t  U'
                '              V\t\tEp\n')
        for k in sorted(results):
            r = results[k]
            N = r['N']
            f.write(f'{N:<5d} '
                    f'{r["R"]*r["x1"]*1e3:<14.9g} {r["R"]*r["x2"]*1e3:<14.9g} '
                    f'{r["R"]*r["x3"]*1e3:<14.9g} {r["T1"]/N:<15.9g} '
                    f'{r["V1"]/N:<15.9g} {r["Volumeenergy"]/N:<12.9g} '
                    f'{r["Energywithoutmass"]/N:<12.9g}\n')

    print('\nwrote eta_results_python.dat and data_paper_python.dat')


if __name__ == '__main__':
    main()
