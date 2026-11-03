# Thomas–Fermi multiquark calculation

Code for *System Energies of Hidden-Charm Pentaquark Families in the Thomas–Fermi
Quark Model*.

The repository holds both the original **Mathematica notebooks** and a **Python
port** of them, together with the parameter sweeps and the plotting scripts.

---

## Repository layout

| | |
|---|---|
| `mathematica/` | the original notebooks (see below) |
| everything else | the Python port, sweeps and plots |

The Python code requires `numpy`, `scipy` and `matplotlib`.

---

## Mathematica notebooks

The calculation runs in two stages, and the two stages use different Mathematica
solvers.

**Stage 1 — consistency conditions.** `mathematica/betas_gammas.nb` solves the
consistency conditions with **`NSolve`** to obtain the proportionality constants
(the beta and gamma parameters) and the Q constants for every eta. These are the
values reported in Table 1 of the paper.

**Stage 2 — the Thomas–Fermi equations.** Each `mathematica/eta_*.nb` takes the
Stage 1 constants as fixed input and solves the three-region Thomas–Fermi problem
with **`NDSolve`**: it locates the region boundaries from the quark-number
conditions, then minimises the total energy over the initial slope `m2` using
`FindMinimum`.

| file | eta | N = 5*eta |
|---|---|---|
| `eta___1.nb`  | 1  | 5  |
| `eta____2.nb` | 2  | 10 |
| `eta____3.nb` | 3  | 15 |
| `eta____4.nb` | 4  | 20 |
| `eta____5.nb` | 5  | 25 |
| `eta____6.nb` | 6  | 30 |
| `eta____7.nb` | 7  | 35 |
| `eta____8.nb` | 8  | 40 |
| `eta____9.nb` | 9  | 45 |
| `eta____10.nb`| 10 | 50 |
| `eta____11.nb`| 11 | 55 |
| `eta____12.nb`| 12 | 60 |
| `eta____13.nb`| 13 | 65 |
| `eta____14.nb`| 14 | 70 |

Filenames keep their original underscore counts (`eta___1.nb` has three, the rest
have four).

These notebooks implement the discontinuous case, in which the light-quark density
is discontinuous at the inner boundary.

---

## Python port

### Core

| file | what it is |
|---|---|
| `tf_core.py` | the solver — integrates the three regions, locates the region boundaries from the quark-number conditions, and computes the K coefficients and the energies. One `Model` class; nothing physics-specific is hardcoded in it. |
| `configs.py` | the 14 parameter sets (n, g, couplings, Q constants, seed). The Q constants and the beta/gamma couplings are the Stage 1 results, taken from the notebooks rather than recomputed; `tf_core.py` consumes them and does not re-solve the consistency conditions. |
| `run_all.py` | reproduces all 14 notebooks and writes `eta_results_python.dat` and `data_paper_python.dat`. |

```
python3 run_all.py            # all 14
python3 run_all.py 1 4 8      # a subset
python3 run_all.py --root     # Brent root-finding instead of FindMinimum
```

`run_all.py` prints, for every eta, the relative difference against the values
stored in that notebook's own saved output.

### Sweeps

| file | what it does |
|---|---|
| `sweep.py` | varies alpha_s and B^(1/4) on a grid (default 0, ±10, ±20, ±30 %). |
| `mass_sweep.py` | varies the quark masses, in three modes: both together, light only, heavy only. |

```
python3 sweep.py --jobs 8                 # 14 eta x 49 parameter sets
python3 sweep.py --span 20 --step 10      # a smaller grid
python3 mass_sweep.py --etas 1,4,8,12 --jobs 8
```

Both write one CSV row per run as it finishes and skip rows already present, so an
interrupted sweep resumes where it stopped.

### Plots

| file | figures |
|---|---|
| `plot_energy_terms.py` | kinetic, potential and volume contributions vs B^(1/4), vs alpha_s, and vs the quark masses |
| `plot_six_panels.py` | total energy and total mass per quark, one panel per sweep mode |
| `plot_sweep_analysis.py` | eta comparison and the energy-vs-N figures |
| `plot_mass_sweep.py`, `plot_energies.py`, `plot_sweep.py` | earlier versions, kept for reference |
| `tidy_sweep.py` | sorts a sweep CSV and builds a pivot table (one row per parameter pair, one column per eta) |

### Data

| file | contents |
|---|---|
| `eta_results_python.dat` | final-loop values for all 14 eta: x1, x2, x3, R, m2, T, V, V_vol, E_p, E_total |
| `data_paper_python.dat` | the same in the column layout of `data_paper.dat` |
| `mass_sweep_results.csv` | the quark-mass sweep |
| `sweep_results*.csv` | the alpha_s / B sweeps |

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
- the notebooks use exact pi in `R` but `N[Pi,7] = 3.141593` in the volume term;
  using pi in both shifts V_vol in the sixth digit.

A `solver` switch on `Model` selects how the region boundaries are located:
`'minimize'` reproduces the notebook's `FindMinimum`-on-`Abs` exactly, `'root'`
uses Brent's method on the un-`Abs`'d condition (same mathematics, quadratic
convergence). Comparing the two separates genuine physics from solver tolerance.

---

## Validation

Every eta was checked against the values stored in that notebook's own saved
output. Running `run_all.py` reproduces the comparison table.

- **E_p** agrees to **2.1e-07 or better** across all 14 configurations.
- `x1` and `m2` drift further in a few cases, up to 7.1e-05 and 1.4e-04 at
  eta = 11. This is the stopping point of `FindMinimum` differing between
  Mathematica and scipy; the energy is stationary at the optimum, so E_p stays
  pinned even when the inputs move in the fifth digit.

Two further checks run inside every evaluation:

- Each region's normalisation is computed both as an integral and as
  `x*y'(x) - y(x)`. The middle and outer regions agree to 1.4e-12 or better. The
  inner region differs by up to **1.1e-07**, and that difference is not solver
  error: it is exactly the boundary term `x0*m2` from starting the integration at
  `x0 = 1e-9` rather than at zero. Predicted and measured values agree to within
  the three figures quoted, and the deviation is linear in `x0`.
- The rest mass recovered as `(E_total - E_p)/N` is 804.800000 MeV at every eta,
  exactly.
