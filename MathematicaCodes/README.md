# Thomas–Fermi Code Map (For Python Files)

What each file in the pentaquark package does, and which values in the calculation are genuinely chosen rather than determined. Numbers throughout are computed from the code itself, at η = 1 — Case I, where the inner region is the anti-charm one. Region labels follow the code, which differ from the manuscript's at this η; the [seed section](#part-2--the-seed) sets out the mapping.

**Contents**

- [Part 1 — Pipeline](#part-1--pipeline)
- [Part 2 — The seed](#part-2--the-seed)

---

## Part 1 — Pipeline

Files flow through four stages: **Imported** (never run directly) → **Run these** (each writes its own output file) → **Produced** (results, so figures redraw without re-solving) → **Figures** (read the files above).

### Stage 1 · Imported — never run directly

#### `tf_core.py`
*Imported module · the solver*

Holds the `Config` dataclass and the `Model` class. For a given m₂ it locates the three region edges, integrates each region, and computes T, U, V_vol and E_p. The m₂ minimisation lives here too, transcribed from the notebook loop.

- **Reads:** nothing — pure computation
- **Writes:** nothing — returns values
- **Note:** imported by every other file

#### `configs.py`
*Imported module · the 14 parameter sets*

One `Config` per η = 1…14: quark numbers and degeneracies, the Q constants, the flavour overlap factors, and the seed m₂. Extracted from the `.nb` files rather than retyped.

- **Reads:** nothing
- **Writes:** nothing
- **Note:** imported alongside tf_core

### Stage 2 · Run these — each writes its own output file

#### `run_all.py`
*Script you run · reproduce all 14 η*

Runs every configuration and prints, for each η, the relative difference against the values saved in that notebook — so it doubles as the check against the Mathematica.

```bash
python3 run_all.py
python3 run_all.py 1 4 8
python3 run_all.py --root
```

- **Reads:** configs.py
- **Writes:** eta_results_python.dat, data_paper_python.dat

#### `sweep.py`
*Script you run · vary α_s and B*

Scans the strong coupling and the bag constant on a grid — by default 0, ±10, ±20, ±30 %, so 49 parameter sets across all 14 η. Writes each row as it finishes and skips rows already present, so an interrupted run resumes where it stopped.

```bash
python3 sweep.py --jobs 8
python3 sweep.py --span 20 --step 10
```

- **Reads:** configs.py
- **Writes:** sweep_results_pm30.csv

#### `mass_sweep.py`
*Script you run · vary the quark masses*

Three modes: both masses scaled together, light only, or heavy only. The light mass enters through a₁ = 1/m_light, so it moves R, c, d and the volume term; the heavy mass enters only through the ratio α₂.

```bash
python3 mass_sweep.py --etas 1,4,8,12 --jobs 8
```

- **Reads:** configs.py
- **Writes:** mass_sweep_results.csv

#### `tidy_sweep.py`
*Script you run · tidy a sweep CSV*

Sorts a sweep file and builds a pivot table — one row per parameter pair, one column per η. Worth running because the raw rows come out in completion order, not parameter order.

```bash
python3 tidy_sweep.py <file.csv>
```

- **Reads:** a sweep CSV
- **Writes:** a sorted CSV and a pivot table

#### `eta2.py`
*Script you run · the first port*

The standalone η = 2 conversion, written before the solver was generalised. Kept as an independent reference: it shares no code with tf_core.py, so agreement between them is a real check.

```bash
python3 eta2.py
```

- **Reads:** nothing
- **Writes:** prints to screen

### Stage 3 · Produced — results, so figures redraw without re-solving

#### `eta_results_python.dat`
*File produced · the 14 η results*

Final values for every configuration: x₁, x₂, x₃, R, m₂, T, U, V_vol, E_p and the total energy. Written by `run_all.py`.

#### `data_paper_python.dat`
*File produced · paper column layout*

The same results arranged in the column order of the original `data_paper.dat`, for direct comparison against the table in the paper. Written by `run_all.py`.

#### `sweep_results*.csv`
*File produced · α_s / B scans*

One row per (α_s, B, η) combination, with the energies and region edges for each. Written by `sweep.py`.

#### `mass_sweep_results.csv`
*File produced · quark-mass scan*

One row per (mode, percentage, η), carrying the scaled masses alongside the resulting energies. Written by `mass_sweep.py`.

### Stage 4 · Figures — read the files above

#### `plot_energy_terms.py`
*Plotting script · the three terms separately*

Kinetic, potential and volume contributions plotted against B^1/4, against α_s, and against the quark masses — the figures with the individual terms replacing the total.

```bash
python3 plot_energy_terms.py
```

- **Reads:** sweep and mass-sweep CSVs
- **Writes:** energy_terms_vs_B, _vs_alphas, _vs_mass

#### `plot_six_panels.py`
*Plotting script · energy and mass per quark*

Total energy and total mass per quark, one panel per sweep mode.

```bash
python3 plot_six_panels.py
```

- **Reads:** mass_sweep_results.csv
- **Writes:** a six-panel figure

#### `plot_sweep_analysis.py`
*Plotting script · η comparison*

The η comparison and the energy-versus-N figures, including the η = 1, 4, 8, 12 curves.

```bash
python3 plot_sweep_analysis.py
```

- **Reads:** a sweep CSV
- **Writes:** sweep comparison figures

#### `plot_sweep.py`, `plot_energies.py`, `plot_mass_sweep.py`
*Plotting script · earlier versions*

Superseded by the three above, kept so older figures can be regenerated exactly as they were first produced.

- **Reads:** the CSVs
- **Writes:** earlier figure versions
- **Note:** kept for reference

---

## Part 2 — The seed

### Why exactly one number has to be chosen

Each region obeys a **second-order** differential equation:

> **Region equation:**  y″(x) = Q · |y(x)|^(3/2) / √x

A second-order equation needs **two** pieces of information at the starting point before it has a unique solution. Count what the physics actually supplies:

| Condition | Status | Why |
|---|---|---|
| **y₂(0) = 0** | Forced by the physics | Since y/x is the potential, a non-zero value gives y/x ~ y₀/x — a 1/r Coulomb tail, which is the field of a point colour source at the centre. No configuration here has one, so this condition is not a choice. |
| **y₂′(0) = m₂** | Nothing fixes it | No physical requirement says what the slope at the origin should be. Every value gives a different, perfectly valid solution of the equation. |

Two needed, one supplied — so exactly one degree of freedom survives, and that leftover is m₂. This is not an artefact of how the code was written. Any method of solving this problem meets the same gap, which is why a *physical* criterion from outside the differential equation is needed to close it. That criterion is energy minimisation: the real configuration is the one with the lowest energy.

### What m₂ means

Because y₂ vanishes at the origin, the definition of the derivative gives the limit directly:

> **With y₂(0) = 0:**  y₂/x → y₂′(0) = m₂  as x → 0
>
> **And since y/x is the potential:**  ρ(0) ∝ m₂^(3/2)

So m₂ is the potential at the centre, and it sets the central density. It is a physical quantity, not a numerical knob — which is why the rest of the solution follows from it.

### Which region is which — read this before the numbers

Every number from here on is computed at **η = 1**, which is Case I. The flavour sitting in each region differs between the two cases, so the symbols below do not line up with the manuscript in the way you might expect.

The two naming schemes differ in kind. The **paper names each boundary by flavour**: x₂ is the charm boundary, x̄₂ the anti-charm boundary and x₁ the light one, always. The **code names them by position**: `x2` is the inner edge, `x3` the middle and `x1` the outer, always. Since the flavour ordering flips between the two cases, the schemes line up differently in each:

| Region | Case I — η = 1 | Case II — η = 2…14 |
|---|---|---|
| inner | anti-charm — paper **x̄₂**, code `x2` | charm — paper x₂, code `x2` |
| middle | charm — paper **x₂**, code `x3` | anti-charm — paper x̄₂, code `x3` |
| outer | light — paper x₁, code `x1` | light — paper x₁, code `x1` |

So where this document writes `x2` at η = 1 it means the **anti-charm** edge, which the paper calls x̄₂ — not the paper's x₂, which is this document's `x3`. Note also that the paper uses only three boundary symbols, x₂, x̄₂ and x₁; there is no x₃ in the manuscript, and the code's `x3` is always the middle edge whichever flavour occupies it. Nothing about the physics or the argument changes; only the labels do.

Written out as ranges: in Case I the inner region spans 0 < x < x̄₂, the middle x̄₂ < x < x₂ and the outer x₂ < x < x₁. In Case II the inner spans 0 < x < x₂, the middle x₂ < x < x̄₂ and the outer x̄₂ < x < x₁.

### Everything downstream depends on m₂

At η = 1 (Case I). Raising m₂ makes the centre denser, and the whole system contracts in response. All three region edges shrink together, and the three energy terms move in different directions.

*Figure 1 — Kinetic T, potential U and volume V_vol (MeV per quark) against m₂.* The competition. Contraction cuts the volume term steeply while driving the kinetic term up; the potential term deepens slowly.

*Figure 2 — Total energy per quark E_p (MeV) against m₂.* Their sum, E_p = T + U + V_vol, has a single minimum. The search converges to m₂ = 12.822; the true vertex of the curve lies at 12.817, 2.4 × 10⁻⁵ MeV lower.

(The full data behind both figures is in the [appendix](#appendix--full-m-scan-data).)

| m₂ | x₂ | x₃ | x₁ | T | U | V_vol | E_p |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8.000 | 0.3679 | 0.3757 | 0.6520 | 69.86 | -32.08 | 138.17 | 175.959 |
| 10.000 | 0.3276 | 0.3345 | 0.5760 | 89.13 | -36.13 | 95.31 | 148.303 |
| 11.000 | 0.3118 | 0.3183 | 0.5467 | 98.83 | -38.01 | 81.45 | 142.270 |
| 12.000 | 0.2981 | 0.3043 | 0.5213 | 108.55 | -39.80 | 70.63 | 139.377 |
| **12.822** | **0.2881** | **0.2941** | **0.5028** | **116.56** | **-41.22** | **63.39** | **138.735** |
| 14.000 | 0.2753 | 0.2810 | 0.4794 | 128.07 | -43.17 | 54.95 | 139.854 |
| 16.000 | 0.2570 | 0.2623 | 0.4461 | 147.69 | -46.30 | 44.27 | 145.662 |
| 18.000 | 0.2419 | 0.2469 | 0.4188 | 167.37 | -49.24 | 36.63 | 154.764 |
| 20.000 | 0.2292 | 0.2339 | 0.3959 | 187.12 | -52.02 | 30.94 | 166.040 |
Every column moves with m₂ — the edges contract by a third across this range, while T nearly triples and V_vol falls to a fifth.

One relation in that table is exact. The volume energy tracks the volume of the outer region, so across the whole range:

> **Constant to four figures at every m₂:**  V_vol / x₁³ = 498.6

That is the mechanism behind the minimum. V_vol falls as the cube of the outer edge, so contraction is cheap at first and the energy drops; but the kinetic term keeps climbing, and past m₂ ≈ 12.8 it wins. The balance between those two is what picks the physical configuration.

### But the seed you start from does not change the answer

m₂ determines the solution — yet the *seed*, the value the search starts from, does not. The search walks downhill to the same minimum from anywhere on the slope. Starting η = 1 from six seeds spanning a factor of more than three:

| Seed m₂⁽⁰⁾ | Converged m₂ | E_p (MeV) | x₁ (outer edge) |
|---:|---:|---:|---:|
| 5, 7, 9, 12 | 12.822 | 138.73462386 | 0.502821170 |
| 15, 18 | 12.788 | 138.73534216 | 0.503546930 |

Seeds below the minimum settle on 12.822, seeds above it on 12.788. The two differ because the search stops as soon as the energy stops falling, and approaching from opposite sides it halts on grid points 0.034 apart — to a search whose final step is 0.001, both are the bottom.

What that costs is the point, and the E_p curve above already shows why it is so little. The bowl is almost flat near its base, so a 0.034 error in m₂ moves the energy by

> **Spread between the two converged values:**  ΔE_p = 7.18 × 10⁻⁴ MeV = 5.2 × 10⁻⁶ × E_p

about five parts in a million, against the tens of MeV that separate one η from the next. The seed has to land the search in the right basin; beyond that it does not reach the result. The same value 45 is used for η = 10 through 14, whose converged answers run from 52.6 to 73.3 — at η = 13 it starts 28.3 away and still arrives.

### The other three starting values are not seeds

Once m₂ is fixed, the solution y(x) is fully determined. The region edges are then not choices but questions about that curve — each has one answer, which the code has to search for.

1. **m₂** — *chosen.* From the config, or from the current position of the outer search.
2. **x₂ — inner edge** — *determined.* The unique x where `∫√x |y₂|^(3/2) dx = n₂/3η`. Search starts at `0.01`, bracketed to `[10⁻⁶, 2.0]`.
3. **x₃ — middle edge** — *determined.* The unique root of the middle-region condition. Starts from `x₂`.
4. **x₁ — outer edge** — *determined.* The unique root of the outer-region condition. Starts from `x₃`.
5. **T, U, V_vol → E_p** — *output.* The energies follow from the three solved regions.

↻ *The outer loop then nudges m₂ and repeats the whole chain — three passes with steps 0.1, 0.01 and 0.001, each continuing while the energy keeps falling.*

### The starting guess cannot reach the answer

At η = 1 (Case I), so `x2` here is the anti-charm edge. Varying where its search begins, over a factor of a thousand, at m₂ = 12.822. Every digit is identical, because the answer was already there before the search started.

| x₂ start | x₂ found | x₁ found | E_p (MeV) |
|---:|---:|---:|---:|
| 0.001 | 0.288060387398 | 0.502821169908 | 138.73462386 |
| 0.005 | 0.288060387398 | 0.502821169908 | 138.73462386 |
| 0.010 | 0.288060387398 | 0.502821169908 | 138.73462386 |
| 0.050 | 0.288060387398 | 0.502821169908 | 138.73462386 |
| 0.200 | 0.288060387398 | 0.502821169908 | 138.73462386 |
| 0.500 | 0.288060387398 | 0.502821169908 | 138.73462386 |
| 1.000 | 0.288060387398 | 0.502821169908 | 138.73462386 |

### Two unrelated algorithms land in the same place

At η = 1 (Case I), code names throughout. Nelder–Mead minimisation of |condition| against Brent root-finding — different methods, agreeing to the last digits. The roots are unique, so there is nothing for a guess to select.

| Method | x₂ | x₃ | x₁ | E_p |
|---|---:|---:|---:|---:|
| `minimize` | 0.288060387398 | 0.294051162776 | 0.502821169908 | 138.73462386 |
| `root` | 0.288060387398 | 0.294051162775 | 0.502821169908 | 138.73462386 |

### In short

| | m₂ | x₁, x₂, x₃ |
|---|---|---|
| Why it exists | the equation needs a second condition physics does not supply | answers to the quark-counting conditions |
| Does the choice change the result? | **Yes** — a different solution | No — the same root every time |
| What settles it | minimising the energy | the condition itself |
| Role of the starting value | a real seed, per η in `configs.py` | scaffolding for the search |

---

## Appendix — full m₂ scan data

Computed from `tf_core.py` at η = 1. This is the data plotted in Figures 1 and 2 (energies in MeV).

| m₂ | x₂ | x₃ | x₁ | T | U | V_vol | E_p |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8.0 | 0.36794 | 0.375743 | 0.65196 | 69.862 | -32.078 | 138.174 | 175.9589 |
| 8.5 | 0.356497 | 0.364036 | 0.630277 | 74.664 | -33.136 | 124.841 | 166.3687 |
| 9.0 | 0.346049 | 0.353348 | 0.610574 | 79.477 | -34.164 | 113.496 | 158.8087 |
| 9.5 | 0.336459 | 0.343539 | 0.592571 | 84.3 | -35.163 | 103.749 | 152.8867 |
| 10.0 | 0.327617 | 0.334496 | 0.576037 | 89.133 | -36.135 | 95.305 | 148.3028 |
| 10.5 | 0.31943 | 0.326124 | 0.560785 | 93.975 | -37.084 | 87.933 | 144.8249 |
| 11.0 | 0.311823 | 0.318345 | 0.546659 | 98.826 | -38.01 | 81.454 | 142.2701 |
| 11.5 | 0.304729 | 0.311092 | 0.533527 | 103.684 | -38.915 | 75.724 | 140.4931 |
| 12.0 | 0.298094 | 0.304308 | 0.521279 | 108.549 | -39.8 | 70.627 | 139.377 |
| 12.5 | 0.29187 | 0.297945 | 0.50982 | 113.421 | -40.667 | 66.072 | 138.8266 |
| 12.822 | 0.28806 | 0.294051 | 0.502821 | 116.563 | -41.216 | 63.388 | 138.7346 |
| 13.0 | 0.286017 | 0.291962 | 0.499071 | 118.3 | -41.516 | 61.98 | 138.7639 |
| 13.5 | 0.280499 | 0.286322 | 0.488961 | 123.185 | -42.349 | 58.289 | 139.1244 |
| 14.0 | 0.275287 | 0.280994 | 0.47943 | 128.075 | -43.167 | 54.946 | 139.8544 |
| 14.5 | 0.270352 | 0.27595 | 0.470424 | 132.971 | -43.97 | 51.908 | 140.9085 |
| 15.0 | 0.265671 | 0.271166 | 0.461899 | 137.872 | -44.76 | 49.136 | 142.2487 |
| 15.5 | 0.261223 | 0.266621 | 0.453812 | 142.778 | -45.536 | 46.601 | 143.8423 |
| 16.0 | 0.256989 | 0.262294 | 0.446127 | 147.688 | -46.299 | 44.273 | 145.6617 |
| 16.5 | 0.252954 | 0.25817 | 0.438814 | 152.603 | -47.051 | 42.131 | 147.6829 |
| 17.0 | 0.249101 | 0.254233 | 0.431842 | 157.522 | -47.792 | 40.155 | 149.8852 |
| 17.5 | 0.245417 | 0.250469 | 0.425186 | 162.445 | -48.521 | 38.327 | 152.2507 |
| 18.0 | 0.241891 | 0.246867 | 0.418824 | 167.372 | -49.24 | 36.632 | 154.7637 |
| 18.5 | 0.238512 | 0.243414 | 0.412734 | 172.303 | -49.95 | 35.057 | 157.4105 |
| 19.0 | 0.235269 | 0.240101 | 0.406897 | 177.237 | -50.649 | 33.591 | 160.179 |
| 19.5 | 0.232154 | 0.236918 | 0.401298 | 182.175 | -51.34 | 32.223 | 163.0585 |
| 20.0 | 0.229159 | 0.233858 | 0.395919 | 187.116 | -52.021 | 30.944 | 166.0397 |
