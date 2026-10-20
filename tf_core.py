#!/usr/bin/env python3
"""
tf_core.py -- Thomas-Fermi multiquark energy calculation.

Shared engine for all 14 eta notebooks. Every notebook implements the same
three-region calculation; they differ only in the constants (see configs.py).

Region naming follows notebooks 2-14 (inner / middle / outer):

    inner   0   -> x2    ODE constant Q2,  occupancy n2, g2
    middle  x2  -> x3    ODE constant Q3,  occupancy n3, g3
    outer   x3  -> x1    ODE constant Q1,  occupancy n1, g1

Notebook 1 uses the names (x2bar, x2, x1) / (Q2bar, Q2, Q1) /
(n2bar g2bar, n2 g2, n1 g1) / (beta12bar, gamma12) for exactly the same
quantities. Its config in configs.py is written in the naming above, so the
same engine runs it -- the K coefficients, T1 and V1 were verified term by
term to be identical under that renaming.

In every region:      y''(x) = Q * |y(x)|^(3/2) / sqrt(x)

Requires: numpy, scipy
"""

from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize, brentq

X0 = 1.0e-9            # inner start point (notebook: 1.*^-9)
RTOL, ATOL = 1e-11, 1e-13
P_VOL = 3.141593       # N[Pi, 7] -- the notebook uses this in the volume term only


@dataclass
class Config:
    name: str
    n1: int; g1: int
    n2: int; g2: int
    n3: int; g3: int
    gamma12bar: float
    beta12: float
    beta22bar: float
    Q1: float
    Q2: float
    Q3: float
    m2_start: float
    xx: float = 0.3114
    zz: float = 306.0          # light quark mass, also sets a1 = 1/zz
    aa: float = 0.197          # alfa2 = m_light / m_heavy
    yy: float = 96.84
    m_light: float = 306.0     # rest mass used in the additive term
    m_heavy: float = 1553.0

    @property
    def eta(self):
        return (self.n1 * self.g1 + self.n3 * self.g3 + self.n2 * self.g2) / 5.0

    @property
    def N(self):
        return int(round(5 * self.eta))


class Model:
    def __init__(self, cfg: Config, solver='minimize'):
        self.cfg = cfg
        self.solver = solver
        c = cfg
        self.alfas = c.xx
        self.a = 1.0 / c.zz                      # notebook sets a = a1 in the loop
        self.B = c.yy ** 4
        self.alfa1 = 1.0
        self.alfa2 = c.aa
        self.eta = c.eta
        g = self.alfas
        self.R = (3 * self.a) / (8 * self.alfas) \
                 * (3 * np.pi * self.eta / 2) ** (2 / 3)
        self.c = (12 / (5 * np.pi)) * (4 / 3) * g ** 2 * (4 / 3) \
                 * (self.alfas / self.a) * (3 * np.pi * self.eta / 2) ** (1 / 3)
        self.d = (4 / 3) * g ** 2 * (4 / 3) * (self.alfas / self.a) \
                 * (3 * np.pi * self.eta / 2) ** (1 / 3) \
                 * (18 / (25 * np.pi)) * (self.eta / (self.eta - 1 / 5))

    # ---------------- ODE over one region ----------------
    @staticmethod
    def _rhs(x, Y, Q):
        return [Y[1], Q * abs(Y[0]) ** 1.5 / np.sqrt(x)]

    def integrate(self, Q, xa, xb, y0, dy0):
        if xb == xa:
            return None, y0, dy0
        sol = solve_ivp(self._rhs, (xa, xb), [y0, dy0], args=(Q,), method='LSODA',
                        dense_output=True, rtol=RTOL, atol=ATOL)
        if not sol.success:
            return sol, np.nan, np.nan
        return sol, sol.y[0, -1], sol.y[1, -1]

    @staticmethod
    def _quad(f, a, b):
        if b == a:
            return 0.0
        v, _ = quad(f, a, b, limit=400, epsabs=1e-12, epsrel=1e-11)
        return v

    # ---------------- consistency conditions ----------------
    def cond_inner(self, xmax, m2):
        c = self.cfg
        if xmax <= X0:
            return np.nan
        sol, _, _ = self.integrate(c.Q2, X0, xmax, 0.0, m2)
        if sol is None or not sol.success:
            return np.nan
        m = self._quad(lambda x: np.sqrt(x) * abs(sol.sol(x)[0]) ** 1.5, X0, xmax)
        return m - c.n2 / (3.0 * self.eta)

    def cond_middle(self, x3, x2, u2, s2):
        c = self.cfg; e = self.eta; b = c.beta22bar
        sol, y, dy = self.integrate(c.Q3, x2, x3, b * u2, b * s2)
        if np.isnan(y):
            return np.nan
        return (x3 * dy - y - c.n3 * c.Q3 / (3 * e) - b * c.n2 * c.Q2 / (3 * e)
                + b ** 1.5 * c.n2 * c.Q3 / (3 * e))

    def cond_outer(self, x1, x3, s3, u3):
        c = self.cfg; e = self.eta; gm = c.gamma12bar; b = c.beta22bar
        sol, y, dy = self.integrate(c.Q1, x3, x1, gm * u3, gm * s3)
        if np.isnan(y):
            return np.nan
        return (x1 * dy - y - c.n1 * c.Q1 / (3 * e)
                - gm * c.n3 * c.Q3 / (3 * e)
                + gm ** 1.5 * c.n3 * c.Q1 / (3 * e)
                - (gm ** 1.5 * b ** 1.5 - c.beta12 ** 1.5) * c.n2 * c.Q1 / (3 * e)
                + gm * b ** 1.5 * c.n2 * c.Q3 / (3 * e)
                - gm * b * c.n2 * c.Q2 / (3 * e))

    def locate(self, cond, start, args, lo=None, hi=None):
        """FindMinimum on Abs[cond]; or Brent root-find if solver='root'."""
        if self.solver == 'root' and lo is not None:
            try:
                flo, fhi = cond(lo, *args), cond(hi, *args)
                if np.isfinite(flo) and np.isfinite(fhi) and flo * fhi < 0:
                    return brentq(lambda v: cond(v, *args), lo, hi,
                                  xtol=1e-14, rtol=1e-15)
            except Exception:
                pass

        def obj(v):
            r = cond(v[0], *args)
            return abs(r) if np.isfinite(r) else 1e6

        res = minimize(obj, [start], method='Nelder-Mead',
                       options=dict(xatol=1e-12, fatol=1e-14, maxiter=2000))
        return float(res.x[0])

    # ---------------- one evaluation at fixed m2 ----------------
    def evaluate(self, m2):
        c = self.cfg
        e = self.eta
        n1, g1, n2, g2, n3, g3 = c.n1, c.g1, c.n2, c.g2, c.n3, c.g3
        Q1, Q2, Q3 = c.Q1, c.Q2, c.Q3
        gm, b12, b22 = c.gamma12bar, c.beta12, c.beta22bar

        x2 = self.locate(self.cond_inner, 0.01, (m2,), 1e-6, 2.0)
        sol2, u2, s2 = self.integrate(Q2, X0, x2, 0.0, m2)
        if u2 < 0.0:
            u2 = 0.0

        x3 = self.locate(self.cond_middle, x2, (x2, u2, s2), x2, x2 + 2.0)
        sol3, u3, s3 = self.integrate(Q3, x2, x3, b22 * u2, b22 * s2)
        if u3 < 0.0:
            u3 = 0.0

        x1 = self.locate(self.cond_outer, x3, (x3, s3, u3), x3, x3 + 2.0)
        sol1, u1, s1 = self.integrate(Q1, x3, x1, gm * u3, gm * s3)
        if u1 < 0.0:
            u1 = 0.0

        # normalisation checks
        normmx2 = x2 * s2 - u2
        normx2 = (n2 * Q2) / (3 * e)
        normmx3 = x3 * s3 - u3
        normx3 = ((n3 * Q3) / (3 * e) + b22 * (n2 * Q2) / (3 * e)
                  - b22 ** 1.5 * (n2 * Q3) / (3 * e))
        normmx1 = x1 * s1 - u1
        normx1 = ((n1 * Q1) / (3 * e) + gm * (n3 * Q3) / (3 * e)
                  - gm ** 1.5 * (n3 * Q1) / (3 * e)
                  + (gm ** 1.5 * b22 ** 1.5 - b12 ** 1.5) * (n2 * Q1) / (3 * e)
                  - gm * b22 ** 1.5 * (n2 * Q3) / (3 * e)
                  + gm * b22 * (n2 * Q2) / (3 * e))

        f = lambda s: (lambda x: abs(s.sol(x)[0]) ** 2.5 / np.sqrt(x))
        h4 = self._quad(f(sol2), X0, x2) if sol2 is not None else 0.0
        h5 = self._quad(f(sol3), x2, x3) if sol3 is not None else 0.0
        h6 = self._quad(f(sol1), x3, x1) if sol1 is not None else 0.0

        # ---- K coefficients (by-hand forms, the ones V1 uses) ----
        K22 = (4 * n2) / (7 * e * Q2) * s2 - (4 / (7 * Q2)) * np.sqrt(x2) * u2 ** 2.5

        K32 = (-(4 / (7 * Q2)) * b22 ** 1.5 * np.sqrt(x2) * u2 ** 2.5
               + (n2 / (3 * e)) * (b22 ** 1.5 / Q2 - b22 / Q3
                                   + 5 * b22 ** 1.5 / (7 * Q2)) * s2
               + (n2 / (3 * e * Q3)) * s3)

        K21 = (-(4 / (7 * Q2)) * b12 ** 1.5 * np.sqrt(x2) * u2 ** 2.5
               + (n2 / (3 * e)) * ((5 / (7 * Q2)) * b12 ** 1.5 + b12 ** 1.5 / Q2
                                   - gm ** 1.5 * b22 / Q3) * s2
               + (n2 / (3 * e)) * (gm ** 1.5 / Q3 - gm / Q1) * s3
               + (n2 / (3 * e * Q1)) * s1)

        K31 = ((4 / 7) * b22 ** 1.5 * (gm ** 1.5 * b22 / Q3 - b12 ** 1.5 / Q2)
               * np.sqrt(x2) * u2 ** 2.5
               - (4 / (7 * Q3)) * gm ** 1.5 * np.sqrt(x3) * u3 ** 2.5
               + (n2 / (3 * e)) * b22
                 * ((5 / (7 * Q2)) * b22 ** 0.5 * b12 ** 1.5
                    - gm ** 1.5 * b22 ** 1.5 / Q3
                    - (5 * Q2 / (7 * Q3 ** 2)) * gm ** 1.5 * b22
                    + (Q2 / Q3 ** 2) * gm ** 1.5 * b22
                    + (1 / Q2) * b12 ** 1.5 * b22 ** 0.5
                    - b12 ** 1.5 / Q3) * s2
               + (1 / (3 * e * Q3))
                 * (n3 * gm ** 1.5 - (n3 * Q3 / Q1) * gm + n2 * b12 ** 1.5
                    - (n2 * Q2 / Q3) * gm ** 1.5 * b22
                    + (5 / 7) * gm ** 1.5
                      * (n3 + (n2 * Q2 / Q3) * b22 - n2 * b22 ** 1.5)) * s3
               + (n3 / (3 * e * Q1)) * s1)

        K33 = ((4 / 7) * b22 ** 2 * (b22 ** 0.5 / Q3 - b22 / Q2)
               * np.sqrt(x2) * u2 ** 2.5
               - (4 / (7 * Q3)) * np.sqrt(x3) * u3 ** 2.5
               + (n2 / (3 * e)) * b22 ** 2
                 * (b22 * (12 / (7 * Q2)) + (2 / 7) * (Q2 / Q3 ** 2)
                    - 2 * b22 ** 0.5 / Q3) * s2
               + ((12 * n3) / (21 * e * Q3)
                  + b22 * (n2 / (3 * e * Q3))
                    * (b22 ** 0.5 - (2 / 7) * (Q2 / Q3)
                       - (5 / 7) * b22 ** 0.5)) * s3)

        K11 = ((4 / 7) * (b22 ** 2.5 * gm ** 3 / Q3 - b12 ** 3 / Q2)
               * np.sqrt(x2) * u2 ** 2.5
               + (4 / 7) * (gm ** 2.5 / Q1 - gm ** 3 / Q3)
                 * np.sqrt(x3) * u3 ** 2.5
               - (4 / (7 * Q1)) * np.sqrt(x1) * u1 ** 2.5
               + ((n2 / (3 * e)) * (b12 ** 3 / Q2
                                    - b12 ** 1.5 * gm ** 1.5 * b22 / Q3
                                    + (Q2 / Q3 ** 2) * gm ** 3 * b22 ** 2
                                    - b12 ** 1.5 * gm ** 1.5 * b22 / Q3)
                  + (5 * n2 / (21 * e)) * (b12 ** 3 / Q2
                                           - (Q2 / Q3 ** 2) * b22 ** 2 * gm ** 3)) * s2
               + ((n2 / (3 * e)) * (b12 ** 1.5 * gm ** 1.5 / Q3
                                    - b12 ** 1.5 * gm / Q1
                                    - b22 ** 1.5 * gm ** 3 / Q3
                                    + b22 ** 1.5 * gm ** 2.5 / Q1
                                    + b12 ** 1.5 * gm ** 1.5 / Q3
                                    - (Q2 / Q3 ** 2) * gm ** 3 * b22
                                    - b12 ** 1.5 * gm / Q1
                                    + gm ** 2.5 * b22 ** 1.5 / Q1
                                    + (Q2 / Q1 ** 2) * gm ** 2 * b22
                                    - (Q3 / Q1 ** 2) * b22 ** 1.5 * gm ** 2)
                  + (n3 / (3 * e)) * (gm ** 3 / Q3 - 2 * gm ** 2.5 / Q1
                                      + (Q3 / Q1 ** 2) * gm ** 2)
                  + (5 / (21 * e)) * ((n3 / Q3) * gm ** 3
                                      + (n2 * Q2 / Q3 ** 2) * gm ** 3 * b22
                                      - (n2 / Q3) * b22 ** 1.5 * gm ** 3
                                      - (n3 * Q3 / Q1 ** 2) * gm ** 2
                                      - (n2 * Q2 / Q1 ** 2) * gm ** 2 * b22
                                      + (n2 * Q3 / Q1 ** 2) * gm ** 2 * b22 ** 1.5)) * s3
               + ((n1 / (3 * e * Q1))
                  - (n3 * Q3 / (3 * e * Q1 ** 2)) * gm
                  + (n3 / (3 * e * Q1)) * gm ** 1.5
                  + (n2 / (3 * e * Q1)) * (b12 ** 1.5 - gm ** 1.5 * b22 ** 1.5
                                           - (Q2 / Q1) * b22 * gm
                                           + (Q3 / Q1) * b22 ** 1.5 * gm)
                  + (5 / (21 * e)) * ((n1 / Q1) + (n3 * Q3 / Q1 ** 2) * gm
                                      - (n3 / Q1) * gm ** 1.5
                                      + (n2 / Q1) * gm ** 1.5 * b22 ** 1.5
                                      - (n2 / Q1) * b12 ** 1.5
                                      - (n2 * Q3 / Q1 ** 2) * b22 ** 1.5 * gm
                                      + (n2 * Q2 / Q1 ** 2) * gm * b22)) * s1)

        # ---- energies ----
        a1_, a2_ = self.alfa1, self.alfa2
        T1 = self.c * ((4 / 7) * (a1_ * g1 * (b12 ** 2.5 - gm ** 2.5 * b22 ** 2.5)
                                  + a2_ * g2) * np.sqrt(x2) * u2 ** 2.5
                       + (4 / 7) * a2_ * g3 * np.sqrt(x3) * u3 ** 2.5
                       + (4 / 7) * a1_ * g1 * np.sqrt(x1) * u1 ** 2.5
                       + (5 * n2 / (21 * e))
                         * (b22 ** 2 * (Q2 / Q3) * (a1_ * g1 * gm ** 2.5 + a2_ * g3)
                            - a1_ * g1 * b12 ** 2.5 - a2_ * g3 * b22 ** 2.5
                            - a2_ * g2) * s2
                       + (5 / (21 * e))
                         * ((gm ** 2 * a1_ * g1 / Q1)
                            * (n3 * Q3 + b22 * n2 * Q2 - b22 ** 1.5 * n2 * Q3)
                            - (a1_ * g1 * gm ** 2.5 + a2_ * g3)
                              * (n3 + b22 * (n2 * Q2 / Q3) - b22 ** 1.5 * n2)) * s3
                       - (5 * a1_ * g1 / (21 * e))
                         * (n1 + gm * (n3 * Q3 / Q1)
                            + gm ** 1.5 * (b22 ** 1.5 * n2 - n3)
                            - b12 ** 1.5 * n2
                            - gm * b22 ** 1.5 * (n2 * Q3 / Q1)
                            + gm * b22 * (n2 * Q2 / Q1)) * s1)

        V1 = -self.d * (((g1 * n1 - 1) / n1) * g1 * K11
                        + ((g3 * n3 - 1) / n3) * g3 * K33
                        + ((g2 * n2 - 1) / n2) * g2 * K22
                        + 2 * g1 * g3 * K31
                        + (4 / 3) * g2 * g1 * K21
                        + (4 / 3) * g2 * g3 * K32)

        Vol = 4 * P_VOL / 3 * self.R ** 3 * x1 ** 3 * self.B
        Ewom = T1 + V1 + Vol
        Etot = (Ewom + n1 * g1 * c.m_light
                + (n3 * g3 + n2 * g2) * c.m_heavy)

        return dict(m2=m2, x1=x1, x2=x2, x3=x3, u1=u1, u2=u2, u3=u3,
                    s1=s1, s2=s2, s3=s3, R=self.R,
                    normmx1=normmx1, normx1=normx1, normmx2=normmx2, normx2=normx2,
                    normmx3=normmx3, normx3=normx3, h4=h4, h5=h5, h6=h6,
                    K11=K11, K21=K21, K22=K22, K31=K31, K32=K32, K33=K33,
                    T1=T1, V1=V1, Volumeenergy=Vol,
                    Energywithoutmass=Ewom, Energy=Etot)

    # ---------------- outer m2 search ----------------
    def run(self, verbose=False):
        """Mirrors For[n,1,3] / While[Energy <= H] in the notebook."""
        m2 = self.cfg.m2_start
        m2B = m2
        Energy, H, F = 1e9, 1e9 + 1, 1e9 + 1
        SN, nloop = 0.0, 1
        Energy11 = None
        last = None

        for n in (1, 2, 3):
            while Energy <= H:
                F, H = H, Energy
                m2 = m2 + SN * 10.0 ** (-n)
                last = self.evaluate(m2)
                Energy = last['Energy']
                if verbose:
                    print(f"    n={n} m2={m2:.6f} Ewom={last['Energywithoutmass']:.6f}")
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
