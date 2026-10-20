# Thomas–Fermi multiquark calculation — Python implementation

A Python port of the 14 Mathematica notebooks (`eta = 1.nb` … `eta = 14.nb`), plus the
parameter sweeps and the plotting scripts.

Requires `numpy`, `scipy` and `matplotlib`.

---

## Core

| file | what it is |
|---|---|
| `tf_core.py` | the solver — the three-region Thomas–Fermi problem, the consistency conditions, the K coefficients, and the energies. One `Model` class; nothing physics-specific is hardcoded in it. |
| `configs.py` | the 14 parameter sets (n, g, couplings, Q constants, seed), extracted from the `.nb` files rather than retyped. |
| `run_all.py` | reproduces all 14 notebooks and writes `eta_results_python.dat` and `data_paper_python.dat`. |

```
python3 run_all.py            # all 14
python3 run_all.py 1 4 8      # a subset
python3 run_all.py --root     # Brent root-finding instead of FindMinimum
```

`run_all.py` prints, for every η, the relative difference against the values
stored in that notebook's own saved output.

## Sweeps

| file | what it does |
|---|---|
| `sweep.py` | varies α_s and B^(1/4) on a grid (default 0, ±10, ±20, ±30 %). |
| `mass_sweep.py` | varies the quark masses, in three modes: both together, light only, heavy only. |

```
python3 sweep.py --jobs 8                 # 14 η × 49 parameter sets
python3 sweep.py --span 20 --step 10      # a smaller grid
python3 mass_sweep.py --etas 1,4,8,12 --jobs 8
```

Both write one CSV row per run as it finishes and skip rows already present, so an
interrupted sweep resumes where it stopped.

## Plots

| file | figures |
|---|---|
| `plot_energy_terms.py` | kinetic, potential and volume contributions vs B^(1/4), vs α_s, and vs the quark masses |
| `plot_six_panels.py` | total energy and total mass per quark, one panel per sweep mode |
| `plot_sweep_analysis.py` | η comparison and the energy-vs-N figures |
| `plot_mass_sweep.py`, `plot_energies.py`, `plot_sweep.py` | earlier versions, kept for reference |
| `tidy_sweep.py` | sorts a sweep CSV and builds a pivot table (one row per parameter pair, one column per η) |

## Data

| file | contents |
|---|---|
| `eta_results_python.dat` | raw final-loop values for all 14 η: x₁, x₂, x₃, R, m₂, T, V, V_vol, E_p, E_total |
| `data_paper_python.dat` | the same in the column layout of `data_paper.dat` |
| `mass_sweep_results.csv` | the quark-mass sweep |
| `sweep_results*.csv` | the α_s / B sweeps |

---

## How the Mathematica constructs were mapped

| Mathematica | Python |
|---|---|
| `NDSolve` | `scipy.integrate.solve_ivp`, LSODA, rtol 1e-11, atol 1e-13 |
| `NIntegrate` | `scipy.integrate.quad` on the dense-output interpolant |
| `FindMinimum` on `Abs[...]` | `scipy.optimize.minimize`, Nelder–Mead (derivative-free, which suits the kink in `Abs`) |
| `For[n,1,3]` / `While[Energy<=H]` | `Model.run()`, transcribed literally including `HX`, `F`, `SN`, `nloop`, `m2B` |

Two details that matter and are easy to get wrong:

- `Abs[y^(3/2)]` means `|y|^(3/2)`, not `|y|` raised afterwards.
- the notebooks use exact π in `R` but `N[Pi,7] = 3.141593` in the volume term; using π in both
  shifts V_vol in the sixth digit.

A `solver` switch on `Model` selects how the region boundaries are located: `'minimize'`
reproduces the notebook's `FindMinimum`-on-`Abs` exactly, `'root'` uses Brent's method on the
un-`Abs`'d consistency condition (same mathematics, quadratic convergence). Comparing the two
separates genuine physics from solver tolerance.

## Validation

Every η was checked against the values in that notebook's own saved output.
E_p agrees to 1 part in 10⁷ or better for all 14. Where x₁ and m₂ drift into the
10⁻⁴ range it is the stopping point of `FindMinimum` differing between Mathematica and
scipy; the energy is stationary at the optimum, so it stays pinned even when the inputs
move in the fifth digit.

Two further checks run inside every evaluation:

- each region's normalisation is computed both as an integral and as `x·y'(x) − y(x)`;
  the two agree to 1e-8 or better.
- the rest mass recovered as `(E_total − E_p)/N` is 804.8 MeV at every η, exactly.
