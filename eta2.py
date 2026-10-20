#!/usr/bin/env python3
"""
Python port of eta____2.nb  --  Thomas-Fermi multiquark energy calculation.

Configuration: n1=2, g1=3, n2=1, g2=2, n3=1, g3=2  ->  N = 10 quarks, eta = 2

Structure mirrors the Mathematica notebook exactly:

  * three radial regions with the Thomas-Fermi equation
        y'' (x) = Q * |y(x)|^(3/2) / sqrt(x)
    solved outward, each region's boundary matched to the previous one
  * the region boundaries x2, x3, x1 are located by MINIMISING |consistency|
    (the notebook uses FindMinimum on an Abs[...] objective; see SOLVER below)
  * an outer search over the shooting parameter m2 = y2'(0+), stepping in
    decreasing decades (n = 1, 2, 3) until the total energy stops decreasing

Requires: numpy, scipy
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize, brentq

# ----------------------------------------------------------------------------
# SOLVER: how the region boundaries are located.
#
#   'minimize'  reproduces the notebook: FindMinimum on |consistency condition|.
#               Faithful, but Abs has a kink at the solution so convergence is
#               limited (this is the source of the FindMinimum::lstol warnings).
#
#   'root'      solves consistency == 0 with Brent's method instead. Same
#               mathematics, much better convergence. Use this to check whether
#               a discrepancy is physics or solver tolerance.
# ----------------------------------------------------------------------------
SOLVER = 'minimize'

# ============================================================================
# Constants  (notebook lines 1-40)
# ============================================================================
n1, g1 = 2, 3
n2, g2 = 1, 2
n3, g3 = 1, 2

gamma12bar = 0.4974508021184714   # light  <-> anti-charm
beta12     = 0.2485750240167111   # light  <-> charm
beta22bar  = 0.882615032977455    # charm  <-> anti-charm

xx    = 0.3114
alfas = xx
zz    = 306.0                     # light quark mass (MeV)
aa    = 0.197
yy    = 96.84
B0    = yy                        # B^(1/4) in MeV;  B = B0**4 is the bag constant
a1    = 1.0 / zz
P     = 3.141593                  # N[Pi, 7] -- the notebook uses this in the volume term
alfa1 = 1.0
alfa2 = aa

eta   = (n1 * g1 + n3 * g3 + n2 * g2) / 5.0
theta = (-6 * eta) / (5 * (5 * eta - 1))

Q1 = -0.5
Q3 = -1.4247837854483154
Q2 = -3.6375498477951878

g = alfas
a = a1                            # notebook sets a = a1 inside the j-loop
B = B0 ** 4

R = (3 * a) / (8 * alfas) * (3 * np.pi * eta / 2) ** (2 / 3)
c = (12 / (5 * np.pi)) * (4 / 3) * g**2 * (4 / 3) * (alfas / a) \
    * (3 * np.pi * eta / 2) ** (1 / 3)
d = (4 / 3) * g**2 * (4 / 3) * (alfas / a) * (3 * np.pi * eta / 2) ** (1 / 3) \
    * (18 / (25 * np.pi)) * (eta / (eta - 1 / 5))

X0   = 1.0e-9                     # inner start point (notebook uses 1.*^-9)
RTOL, ATOL = 1e-11, 1e-13

# ============================================================================
# The Thomas-Fermi ODE and its solution over one region
# ============================================================================
def _rhs(x, Y, Q):
    # y'' = Q |y|^(3/2) / sqrt(x)      (Abs[y^(3/2)] in the notebook == |y|^(3/2))
    return [Y[1], Q * abs(Y[0]) ** 1.5 / np.sqrt(x)]


def integrate(Q, x_start, x_end, y0, dy0):
    """Solve the region ODE from x_start to x_end. Returns the solution object."""
    if x_end == x_start:
        return None, y0, dy0
    sol = solve_ivp(_rhs, (x_start, x_end), [y0, dy0], args=(Q,),
                    method='LSODA', dense_output=True, rtol=RTOL, atol=ATOL)
    if not sol.success:
        return sol, np.nan, np.nan
    return sol, sol.y[0, -1], sol.y[1, -1]


def region_integral(sol, x_start, x_end):
    """NIntegrate[ y(x)^2.5 / sqrt(x), {x, x_start, x_end} ]"""
    if sol is None or x_end == x_start:
        return 0.0
    f = lambda x: abs(sol.sol(x)[0]) ** 2.5 / np.sqrt(x)
    val, _ = quad(f, x_start, x_end, limit=400, epsabs=1e-12, epsrel=1e-11)
    return val


# ============================================================================
# Consistency conditions -- the notebook's f2, f3, f1
# ============================================================================
def cond_f2(xmax, m2):
    """Inner region. Notebook f2[xmax_, m2_]: |Integral - n2/(3 eta)|."""
    if xmax <= X0:
        return np.nan
    sol, _, _ = integrate(Q2, X0, xmax, 0.0, m2)
    if sol is None or not sol.success:
        return np.nan
    f = lambda x: np.sqrt(x) * abs(sol.sol(x)[0]) ** 1.5
    m, _ = quad(f, X0, xmax, limit=400, epsabs=1e-12, epsrel=1e-11)
    return m - n2 / (3.0 * eta)


def cond_f3(x3, x2, u2, s2):
    """Middle region. Notebook f3[x3_, x2_, u2_, s2_]."""
    sol, y, dy = integrate(Q3, x2, x3, beta22bar * u2, beta22bar * s2)
    if np.isnan(y):
        return np.nan
    return (x3 * dy - y
            - n3 * Q3 / (3 * eta)
            - beta22bar * n2 * Q2 / (3 * eta)
            + beta22bar ** 1.5 * n2 * Q3 / (3 * eta))


def cond_f1(x1, x3, s3, u3):
    """Outer region. Notebook f1[x1_, x3_, s3_, u3_]."""
    sol, y, dy = integrate(Q1, x3, x1, gamma12bar * u3, gamma12bar * s3)
    if np.isnan(y):
        return np.nan
    return (x1 * dy - y
            - n1 * Q1 / (3 * eta)
            - gamma12bar * n3 * Q3 / (3 * eta)
            + gamma12bar ** 1.5 * n3 * Q1 / (3 * eta)
            - (gamma12bar ** 1.5 * beta22bar ** 1.5 - beta12 ** 1.5)
              * n2 * Q1 / (3 * eta)
            + gamma12bar * beta22bar ** 1.5 * n2 * Q3 / (3 * eta)
            - gamma12bar * beta22bar * n2 * Q2 / (3 * eta))


def locate(cond, start, args, lo=None, hi=None):
    """Find where cond(...) == 0, mirroring FindMinimum on Abs[cond]."""
    obj = lambda v: abs(cond(v[0], *args)) if np.isfinite(cond(v[0], *args)) else 1e6

    if SOLVER == 'root' and lo is not None and hi is not None:
        try:
            flo, fhi = cond(lo, *args), cond(hi, *args)
            if np.isfinite(flo) and np.isfinite(fhi) and flo * fhi < 0:
                return brentq(lambda v: cond(v, *args), lo, hi, xtol=1e-14, rtol=1e-14)
        except Exception:
            pass

    res = minimize(obj, [start], method='Nelder-Mead',
                   options=dict(xatol=1e-12, fatol=1e-14, maxiter=2000))
    return float(res.x[0])


# ============================================================================
# One full evaluation at a given m2  -- the body of the notebook's While loop
# ============================================================================
def evaluate(m2):
    out = {'m2': m2}

    # ---- inner region: locate x2, then solve on [X0, x2] ----
    x2 = locate(cond_f2, 0.01, (m2,), lo=1e-6, hi=2.0)
    sol2, u2, s2 = integrate(Q2, X0, x2, 0.0, m2)
    if u2 < 0.0:
        u2 = 0.0

    # ---- middle region: locate x3 ----
    x3 = locate(cond_f3, x2, (x2, u2, s2), lo=x2, hi=x2 + 2.0)
    sol3, u3, s3 = integrate(Q3, x2, x3, beta22bar * u2, beta22bar * s2)
    if u3 < 0.0:
        u3 = 0.0

    # ---- outer region: locate x1 ----
    x1 = locate(cond_f1, x3, (x3, s3, u3), lo=x3, hi=x3 + 2.0)
    sol1, u1, s1 = integrate(Q1, x3, x1, gamma12bar * u3, gamma12bar * s3)
    if u1 < 0.0:
        u1 = 0.0

    # ---- normalisation checks (computational vs by-hand) ----
    normmx2 = x2 * s2 - u2
    normx2 = (n2 * Q2) / (3 * eta)
    normmx3 = x3 * s3 - u3
    normx3 = ((n3 * Q3) / (3 * eta) + beta22bar * (n2 * Q2) / (3 * eta)
              - beta22bar ** 1.5 * (n2 * Q3) / (3 * eta))
    normmx1 = x1 * s1 - u1
    normx1 = ((n1 * Q1) / (3 * eta)
              + gamma12bar * (n3 * Q3) / (3 * eta)
              - gamma12bar ** 1.5 * (n3 * Q1) / (3 * eta)
              + (gamma12bar ** 1.5 * beta22bar ** 1.5 - beta12 ** 1.5)
                * (n2 * Q1) / (3 * eta)
              - gamma12bar * beta22bar ** 1.5 * (n2 * Q3) / (3 * eta)
              + gamma12bar * beta22bar * (n2 * Q2) / (3 * eta))

    # ---- region integrals h4, h5, h6 ----
    h4 = region_integral(sol2, X0, x2)
    h5 = region_integral(sol3, x2, x3)
    h6 = region_integral(sol1, x3, x1)

    # ---- K coefficients (by-hand forms used in the energies) ----
    K22 = (4 * n2) / (7 * eta * Q2) * s2 - (4 / (7 * Q2)) * np.sqrt(x2) * u2 ** 2.5

    K32 = (-(4 / (7 * Q2)) * beta22bar ** 1.5 * np.sqrt(x2) * u2 ** 2.5
           + (n2 / (3 * eta)) * (beta22bar ** 1.5 / Q2 - beta22bar / Q3
                                 + 5 * beta22bar ** 1.5 / (7 * Q2)) * s2
           + (n2 / (3 * eta * Q3)) * s3)

    K21 = (-(4 / (7 * Q2)) * beta12 ** 1.5 * np.sqrt(x2) * u2 ** 2.5
           + (n2 / (3 * eta)) * ((5 / (7 * Q2)) * beta12 ** 1.5
                                 + beta12 ** 1.5 / Q2
                                 - gamma12bar ** 1.5 * beta22bar / Q3) * s2
           + (n2 / (3 * eta)) * (gamma12bar ** 1.5 / Q3 - gamma12bar / Q1) * s3
           + (n2 / (3 * eta * Q1)) * s1)

    K31 = ((4 / 7) * beta22bar ** 1.5
           * (gamma12bar ** 1.5 * beta22bar / Q3 - beta12 ** 1.5 / Q2)
           * np.sqrt(x2) * u2 ** 2.5
           - (4 / (7 * Q3)) * gamma12bar ** 1.5 * np.sqrt(x3) * u3 ** 2.5
           + (n2 / (3 * eta)) * beta22bar
             * ((5 / (7 * Q2)) * beta22bar ** 0.5 * beta12 ** 1.5
                - gamma12bar ** 1.5 * beta22bar ** 1.5 / Q3
                - (5 * Q2 / (7 * Q3 ** 2)) * gamma12bar ** 1.5 * beta22bar
                + (Q2 / Q3 ** 2) * gamma12bar ** 1.5 * beta22bar
                + (1 / Q2) * beta12 ** 1.5 * beta22bar ** 0.5
                - beta12 ** 1.5 / Q3) * s2
           + (1 / (3 * eta * Q3))
             * (n3 * gamma12bar ** 1.5
                - (n3 * Q3 / Q1) * gamma12bar
                + n2 * beta12 ** 1.5
                - (n2 * Q2 / Q3) * gamma12bar ** 1.5 * beta22bar
                + (5 / 7) * gamma12bar ** 1.5
                  * (n3 + (n2 * Q2 / Q3) * beta22bar - n2 * beta22bar ** 1.5)) * s3
           + (n3 / (3 * eta * Q1)) * s1)

    K33 = ((4 / 7) * beta22bar ** 2 * (beta22bar ** 0.5 / Q3 - beta22bar / Q2)
           * np.sqrt(x2) * u2 ** 2.5
           - (4 / (7 * Q3)) * np.sqrt(x3) * u3 ** 2.5
           + (n2 / (3 * eta)) * beta22bar ** 2
             * (beta22bar * (12 / (7 * Q2)) + (2 / 7) * (Q2 / Q3 ** 2)
                - 2 * beta22bar ** 0.5 / Q3) * s2
           + ((12 * n3) / (21 * eta * Q3)
              + beta22bar * (n2 / (3 * eta * Q3))
                * (beta22bar ** 0.5 - (2 / 7) * (Q2 / Q3)
                   - (5 / 7) * beta22bar ** 0.5)) * s3)

    K11 = ((4 / 7) * (beta22bar ** 2.5 * gamma12bar ** 3 / Q3 - beta12 ** 3 / Q2)
           * np.sqrt(x2) * u2 ** 2.5
           + (4 / 7) * (gamma12bar ** 2.5 / Q1 - gamma12bar ** 3 / Q3)
             * np.sqrt(x3) * u3 ** 2.5
           - (4 / (7 * Q1)) * np.sqrt(x1) * u1 ** 2.5
           + ((n2 / (3 * eta)) * (beta12 ** 3 / Q2
                                  - beta12 ** 1.5 * gamma12bar ** 1.5 * beta22bar / Q3
                                  + (Q2 / Q3 ** 2) * gamma12bar ** 3 * beta22bar ** 2
                                  - beta12 ** 1.5 * gamma12bar ** 1.5 * beta22bar / Q3)
              + (5 * n2 / (21 * eta)) * (beta12 ** 3 / Q2
                                         - (Q2 / Q3 ** 2) * beta22bar ** 2 * gamma12bar ** 3)) * s2
           + ((n2 / (3 * eta)) * (beta12 ** 1.5 * gamma12bar ** 1.5 / Q3
                                  - beta12 ** 1.5 * gamma12bar / Q1
                                  - beta22bar ** 1.5 * gamma12bar ** 3 / Q3
                                  + beta22bar ** 1.5 * gamma12bar ** 2.5 / Q1
                                  + beta12 ** 1.5 * gamma12bar ** 1.5 / Q3
                                  - (Q2 / Q3 ** 2) * gamma12bar ** 3 * beta22bar
                                  - beta12 ** 1.5 * gamma12bar / Q1
                                  + gamma12bar ** 2.5 * beta22bar ** 1.5 / Q1
                                  + (Q2 / Q1 ** 2) * gamma12bar ** 2 * beta22bar
                                  - (Q3 / Q1 ** 2) * beta22bar ** 1.5 * gamma12bar ** 2)
              + (n3 / (3 * eta)) * (gamma12bar ** 3 / Q3
                                    - 2 * gamma12bar ** 2.5 / Q1
                                    + (Q3 / Q1 ** 2) * gamma12bar ** 2)
              + (5 / (21 * eta)) * ((n3 / Q3) * gamma12bar ** 3
                                    + (n2 * Q2 / Q3 ** 2) * gamma12bar ** 3 * beta22bar
                                    - (n2 / Q3) * beta22bar ** 1.5 * gamma12bar ** 3
                                    - (n3 * Q3 / Q1 ** 2) * gamma12bar ** 2
                                    - (n2 * Q2 / Q1 ** 2) * gamma12bar ** 2 * beta22bar
                                    + (n2 * Q3 / Q1 ** 2) * gamma12bar ** 2 * beta22bar ** 1.5)) * s3
           + ((n1 / (3 * eta * Q1))
              - (n3 * Q3 / (3 * eta * Q1 ** 2)) * gamma12bar
              + (n3 / (3 * eta * Q1)) * gamma12bar ** 1.5
              + (n2 / (3 * eta * Q1)) * (beta12 ** 1.5
                                         - gamma12bar ** 1.5 * beta22bar ** 1.5
                                         - (Q2 / Q1) * beta22bar * gamma12bar
                                         + (Q3 / Q1) * beta22bar ** 1.5 * gamma12bar)
              + (5 / (21 * eta)) * ((n1 / Q1)
                                    + (n3 * Q3 / Q1 ** 2) * gamma12bar
                                    - (n3 / Q1) * gamma12bar ** 1.5
                                    + (n2 / Q1) * gamma12bar ** 1.5 * beta22bar ** 1.5
                                    - (n2 / Q1) * beta12 ** 1.5
                                    - (n2 * Q3 / Q1 ** 2) * beta22bar ** 1.5 * gamma12bar
                                    + (n2 * Q2 / Q1 ** 2) * gamma12bar * beta22bar)) * s1)

    # ---- energies ----
    T1 = c * ((4 / 7) * (alfa1 * g1 * (beta12 ** 2.5
                                       - gamma12bar ** 2.5 * beta22bar ** 2.5)
                         + alfa2 * g2) * np.sqrt(x2) * u2 ** 2.5
              + (4 / 7) * alfa2 * g3 * np.sqrt(x3) * u3 ** 2.5
              + (4 / 7) * alfa1 * g1 * np.sqrt(x1) * u1 ** 2.5
              + (5 * n2 / (21 * eta))
                * (beta22bar ** 2 * (Q2 / Q3)
                   * (alfa1 * g1 * gamma12bar ** 2.5 + alfa2 * g3)
                   - alfa1 * g1 * beta12 ** 2.5
                   - alfa2 * g3 * beta22bar ** 2.5
                   - alfa2 * g2) * s2
              + (5 / (21 * eta))
                * ((gamma12bar ** 2 * alfa1 * g1 / Q1)
                   * (n3 * Q3 + beta22bar * n2 * Q2 - beta22bar ** 1.5 * n2 * Q3)
                   - (alfa1 * g1 * gamma12bar ** 2.5 + alfa2 * g3)
                     * (n3 + beta22bar * (n2 * Q2 / Q3) - beta22bar ** 1.5 * n2)) * s3
              - (5 * alfa1 * g1 / (21 * eta))
                * (n1 + gamma12bar * (n3 * Q3 / Q1)
                   + gamma12bar ** 1.5 * (beta22bar ** 1.5 * n2 - n3)
                   - beta12 ** 1.5 * n2
                   - gamma12bar * beta22bar ** 1.5 * (n2 * Q3 / Q1)
                   + gamma12bar * beta22bar * (n2 * Q2 / Q1)) * s1)

    V1 = -d * (((g1 * n1 - 1) / n1) * g1 * K11
               + ((g3 * n3 - 1) / n3) * g3 * K33
               + ((g2 * n2 - 1) / n2) * g2 * K22
               + 2 * g1 * g3 * K31
               + (4 / 3) * g2 * g1 * K21
               + (4 / 3) * g2 * g3 * K32)

    Volumeenergy = 4 * P / 3 * R ** 3 * x1 ** 3 * B
    Energywithoutmass = T1 + V1 + Volumeenergy
    Energy6bare = Energywithoutmass + n1 * g1 * 306 + n3 * g3 * 1553 + n2 * g2 * 1553

    out.update(x1=x1, x2=x2, x3=x3, u1=u1, u2=u2, u3=u3, s1=s1, s2=s2, s3=s3,
               normmx1=normmx1, normx1=normx1, normmx2=normmx2, normx2=normx2,
               normmx3=normmx3, normx3=normx3, h4=h4, h5=h5, h6=h6,
               K11=K11, K21=K21, K22=K22, K31=K31, K32=K32, K33=K33,
               T1=T1, V1=V1, Volumeenergy=Volumeenergy,
               Energywithoutmass=Energywithoutmass, Energy=Energy6bare)
    return out


# ============================================================================
# Outer m2 search  -- the notebook's For[n,...] / While[Energy <= H,...]
# ============================================================================
def run(m2_start=14.0, verbose=True):
    m2 = m2_start
    m2B = m2
    Energy = 1e9
    H = 1e9 + 1
    SN = 0.0
    nloop = 1
    F = H
    Energy11 = None
    HX = None
    last = None

    for n in (1, 2, 3):
        while Energy <= H:
            F = H
            H = Energy
            HX = last['Energy'] if last else None

            m2 = m2 + SN * 10.0 ** (-n)
            last = evaluate(m2)
            Energy = last['Energy']

            if verbose:
                print(f"  n={n}  m2={m2:.9f}  Ewom={last['Energywithoutmass']:.6f}  "
                      f"x1={last['x1']:.9f}  x2={last['x2']:.9f}  x3={last['x3']:.9f}")

            if n == 1:
                if nloop == 1:
                    nloop += 1
                    SN = 1.0e-1
                    Energy11 = Energy
                elif nloop == 2:
                    nloop += 1
                    H = 1e9 + 1
                    m2 = m2B
                    SN = 1.0 if (Energy - Energy11) < 0 else -1.0

        if n < 3:
            m2 = m2 - SN * 10.0 ** (-n)
            Energy = F
            H = F

    return last


if __name__ == '__main__':
    print(f'eta = {eta}   N = {int(5*eta)}   R = {R:.12f}   B = B0^4 = {B:.6e}')
    print(f'solver for region boundaries: {SOLVER}\n')

    res = run()

    print('\n' + '=' * 62)
    print('FINAL VALUES (last loop)')
    print('=' * 62)
    for k in ['x1', 'x2', 'x3', 'm2']:
        print(f'  {k:<20s} {res[k]:.9f}')
    print()
    for k, lab in [('T1', 'Kinetic Energy'), ('V1', 'Potential Energy'),
                   ('Volumeenergy', 'Volume Energy'),
                   ('Energywithoutmass', 'Energy without mass'),
                   ('Energy', 'Energy')]:
        print(f'  {lab:<20s} {res[k]:.6f}')
    print()
    print('  normalisation checks (computational vs by hand):')
    for a_, b_ in [('normmx1', 'normx1'), ('normmx2', 'normx2'), ('normmx3', 'normx3')]:
        print(f'    {a_:<9s} {res[a_]:+.9f}   {b_:<8s} {res[b_]:+.9f}   '
              f'rel diff {abs(res[a_]-res[b_])/abs(res[b_]):.2e}')

    # ---- reference values read out of the notebook's own saved output ----
    ref = dict(x1=0.398164724, x2=0.192907181, x3=0.206159501, m2=17.922,
               T1=225.354997, V1=-69.808330, Volumeenergy=125.896028,
               Energywithoutmass=281.442696, Energy=8329.4427)
    print('\n' + '=' * 62)
    print('COMPARISON WITH eta____2.nb SAVED OUTPUT')
    print('=' * 62)
    print(f'  {"quantity":<20s} {"python":>16s} {"notebook":>16s} {"rel diff":>11s}')
    for k, v in ref.items():
        p = res[k]
        print(f'  {k:<20s} {p:16.6f} {v:16.6f} {abs(p-v)/abs(v):11.2e}')
