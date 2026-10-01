#!/usr/bin/env python3
"""NLO DIS (1+1 Born) with 1-jettiness slicing, as a test of the SCET
ingredients before the 2+1 case.

Observable: tau_1^b of Kang, Lee, Stewart (arXiv:1303.6952), q_B = xP,
q_J = q + xP, i.e. DIS thrust in the Breit frame,
    tau = (2/Q^2) sum_i min(q_B.p_i, q_J.p_i).
At fixed (x, Q^2, y), photon exchange, massless quarks, mu_R = mu_F = Q:

    sigma_NLO = sigma_0 * Sigma_c(tau_cut)          [SCET cumulant, KLS (173)+(174)]
              + int_{tau > tau_cut} dsigma_real      [O(alpha_s) 2+1 tree]
              + O(tau_cut log tau_cut),

compared with the exact NLO from the MS-bar coefficient functions (and with
hoppet's structure functions as a check of those). Everything is in units of
2 pi alpha^2 / Q^4, i.e. R = (1+(1-y)^2) F2/x - y^2 FL/x.

For two final-state partons with x_p = x/xi and z = z_p (Breit frame)
    tau = min(z, a(1-z)) + min(1-z, a z),   a = (1 - x_p)/x_p,
which is piecewise linear in z with kinks at a/(1+a) and 1/(1+a).
"""
import argparse, math
import numpy as np
from scipy import integrate
import lhapdf

CF, TR, CA = 4.0 / 3.0, 0.5, 3.0
PI2 = math.pi ** 2
EQ2 = {1: 1 / 9, 2: 4 / 9, 3: 1 / 9, 4: 4 / 9, 5: 1 / 9}   # d u s c b
QUAD = dict(epsabs=0.0, epsrel=1e-11, limit=400)


class Pdfs:
    def __init__(self, name, Q):
        self.p = lhapdf.mkPDF(name, 0)
        self.Q = Q
        self.as_ = self.p.alphasQ(Q)

    only = None           # 'q' or 'g' to switch off the other channel (diagnostics)

    def q(self, xi):      # sum_q e_q^2 (q + qbar)(xi), number densities
        if xi >= 1.0 or self.only == 'g':
            return 0.0
        return sum(e * (self.p.xfxQ(i, xi, self.Q) + self.p.xfxQ(-i, xi, self.Q)) for i, e in EQ2.items()) / xi

    def g(self, xi):      # sum_q e_q^2 * g(xi): one gluon term per flavour (q and qbar counted in the kernels)
        if xi >= 1.0 or self.only == 'q':
            return 0.0
        return sum(EQ2.values()) * self.p.xfxQ(21, xi, self.Q) / xi


def conv(f, x, reg, plus=None, delta=0.0):
    """int_x^1 dz/z f(x/z) [reg(z) + plus(z)_+ ] + delta f(x); plus(z) is the
    function under [.]_+ on [0,1] (singular at z = 1)."""
    res = integrate.quad(lambda z: f(x / z) / z * reg(z), x, 1.0, **QUAD)[0] if reg else 0.0
    if plus:
        h1 = f(x)
        res += integrate.quad(lambda z: plus(z) * (f(x / z) / z - h1), x, 1.0, **QUAD)[0]
        res -= h1 * integrate.quad(plus, 0.0, x, **QUAD)[0]
    return res + delta * f(x)


def exact_nlo(P, x, y):
    """O(alpha_s) part of R from the MS-bar coefficient functions (alpha_s/2pi normalisation)."""
    a = P.as_ / (2 * math.pi)
    Yp = 1 + (1 - y) ** 2
    # F2/x: quark C_q (per q and qbar), gluon C_g per q and per qbar -> 2 C_g per flavour
    c2q = CF * conv(P.q, x, lambda z: -(1 + z) * math.log(1 - z) - (1 + z * z) / (1 - z) * math.log(z) + 3 + 2 * z,
                    plus=lambda z: (2 * math.log(1 - z) - 1.5) / (1 - z), delta=-(4.5 + PI2 / 3))
    c2g = 2 * TR * conv(P.g, x, lambda z: (z * z + (1 - z) ** 2) * math.log((1 - z) / z) - 8 * z * z + 8 * z - 1)
    cLq = CF * conv(P.q, x, lambda z: 2 * z)
    cLg = 2 * TR * conv(P.g, x, lambda z: 4 * z * (1 - z))
    return a * (Yp * (c2q + c2g) - y * y * (cLq + cLg)), a * (c2q + c2g), a * (cLq + cLg)


def cumulant(P, x, y, tau, taua=False):
    """sigma_0 Sigma_c^(1)(tau) for tau_1^b at mu = Q (KLS (173) + (174)), O(alpha_s) part, in units of R.
    taua=True: tau_1^a (KLS (173) alone, i.e. without the ln z terms of (174))."""
    a = P.as_ / (2 * math.pi)
    lzb = 0.0 if taua else 1.0
    Yp = 1 + (1 - y) ** 2
    L = math.log(tau)
    # quark: delta term -(CF/2)(9 + 2pi^2/3 + 6L + 4L^2) [alpha_s/4pi -> alpha_s/2pi: factor 1/2]
    dq = -0.5 * CF * (9 + 2 * PI2 / 3 + 6 * L + 4 * L * L)
    # CF [ L1(1-z)(1+z^2) + (1 - z - (1+z^2)/(1-z) ln z) + ln(tau z) Pqq(z) ],  Pqq = [(1+z^2)/(1-z)]_+
    # ln(tau z) Pqq = L [(1+z^2)/(1-z)]_+ + (1+z^2) ln z/(1-z)  (regular), and L1(1-z)(1+z^2): [ln(1-z)/(1-z)]_+ (1+z^2)
    # Write [g]_+ h with h(1) != 0 via conv(plus=g*h): [g]_+ h = [g h]_+ + h(1) delta * int... handled below.
    def reg(z):
        return 1 - z - (1 + z * z) / (1 - z) * math.log(z) + lzb * (1 + z * z) * math.log(z) / (1 - z)
    # (1+z^2)[ln(1-z)/(1-z)]_+ = [(1+z^2) ln(1-z)/(1-z)]_+ + delta(1-z) int_0^1 ((1+z^2) - 2) ln(1-z)/(1-z) dz,
    #   int_0^1 (z^2 - 1) ln(1-z)/(1-z) dz = -int_0^1 (1+z) ln(1-z) dz = +7/4
    cq = CF * (conv(P.q, x, reg, plus=lambda z: ((1 + z * z) * math.log(1 - z) + L * (1 + z * z)) / (1 - z),
                    delta=+7 / 4))
    # note: [(1+z^2)/(1-z)]_+ L is already a pure plus distribution (Pqq without the 3/2 delta? KLS: Pqq = [(1+z^2)/(1-z)]_+)
    quark = Yp * (dq * P.q(x) + cq)
    # gluon: TF [ ln(tau (1-z)) Pqg(z) + 2 z(1-z) ] per q and per qbar
    glu = Yp * 2 * TR * conv(P.g, x, lambda z: (math.log(tau * (1 - z)) - (1 - lzb) * math.log(z)) * (z * z + (1 - z) ** 2)
                             + 2 * z * (1 - z))
    return a * (quark + glu)


def tau1b(xp, z):
    a = (1 - xp) / xp
    return min(z, a * (1 - z)) + min(1 - z, a * z)


def z_intervals(xp, tc):
    """Sub-intervals of [0,1] with tau_1^b(xp, z) > tc (tau is piecewise linear with kinks z1, z2)."""
    a = (1 - xp) / xp
    knots = sorted({0.0, 1.0, a / (1 + a), 1 / (1 + a)})
    out = []
    for lo, hi in zip(knots[:-1], knots[1:]):
        tl, th = tau1b(xp, lo), tau1b(xp, hi)
        if tl > tc and th > tc:
            out.append((lo, hi))
        elif tl > tc or th > tc:   # one crossing on a linear piece
            zc = lo + (tc - tl) * (hi - lo) / (th - tl)
            out.append((lo, zc) if tl > tc else (zc, hi))
    merged = []
    for iv in out:
        if merged and abs(merged[-1][1] - iv[0]) < 1e-15:
            merged[-1] = (merged[-1][0], iv[1])
        else:
            merged.append(iv)
    return merged


def real_above(P, x, y, tc):
    """O(alpha_s) real emission with tau_1^b > tc, in units of R (alpha_s/2pi normalisation).
    Quark channel (q and qbar):  2xF1-type CF[(xp^2+z^2)/((1-xp)(1-z)) + 2 xp z + 2] ... written as F2 and FL:
      F2: CF[(xp^2+z^2)/((1-xp)(1-z)) + 2 + 6 xp z],   FL: CF 4 xp z
    Gluon channel (per flavour, z over [0,1] covers q and qbar):
      F2: TR[(xp^2+(1-xp)^2)(z^2+(1-z)^2)/(z(1-z)) + 8 xp(1-xp)],   FL: TR 8 xp(1-xp)
    (the gluon F2 constant 8 xp(1-xp), i.e. no polynomial term in the transverse part, follows from requiring
    the tau -> 0 limit with the standard C_g and the beam-function constant 2z(1-z); the first version had 16,
    which gave an offset of exactly 8 TR xp(1-xp) per flavour; to be checked pointwise against DISENT's MATTHR)
    (FL: int dz gives CF 2 xp and 8 TR xp(1-xp) per flavour = Altarelli-Martinelli).
    Returns (R part, F2/x part, FL/x part)."""
    a = P.as_ / (2 * math.pi)
    Yp = 1 + (1 - y) ** 2
    xmax = 1 / (1 + tc)       # above this, tau_1^b <= a < tc everywhere

    k2q = lambda xp, z: CF * ((xp * xp + z * z) / ((1 - xp) * (1 - z)) + 2 + 6 * xp * z)
    kLq = lambda xp, z: CF * 4 * xp * z
    k2g = lambda xp, z: TR * ((xp * xp + (1 - xp) ** 2) * (z * z + (1 - z) ** 2) / (z * (1 - z)) + 8 * xp * (1 - xp))
    kLg = lambda xp, z: TR * 8 * xp * (1 - xp)

    def inner(xp, which):
        kq, kg = (k2q, k2g) if which == 2 else (kLq, kLg)
        sq = sg = 0.0
        for lo, hi in z_intervals(xp, tc):
            sq += integrate.quad(lambda z: kq(xp, z), lo, hi, **QUAD)[0]
            sg += integrate.quad(lambda z: kg(xp, z), lo, hi, **QUAD)[0]
        return (P.q(x / xp) * sq + P.g(x / xp) * sg) / xp

    if xmax <= x:
        return 0.0, 0.0, 0.0
    # integrate in u = -ln(1 - xp) to resolve the xp -> 1 region
    ulo, uhi = -math.log(1 - x), -math.log(1 - xmax)
    res = []
    for which in (2, 0):
        f = lambda u: inner(1 - math.exp(-u), which) * math.exp(-u)
        res.append(a * integrate.quad(f, ulo, uhi, epsabs=0.0, epsrel=1e-9, limit=400)[0])
    f2, fl = res
    return Yp * f2 - y * y * fl, f2, fl


def hoppet_nlo(name, x, Q):
    """O(alpha_s) F2/x and FL/x for photon exchange from hoppet (same machinery as disorder)."""
    import hoppet as hp
    p = lhapdf.mkPDF(name, 0)
    ids = [-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6]
    hp.Start(0.05, 2)
    hp.SetCoupling(p.alphasQ(91.1876), 91.1876, 2)
    hp.Assign(lambda xx, qq: [p.xfxQ(i, xx, qq) for i in ids])
    hp.StartStrFct(2)
    hp.InitStrFct(2, True, 1.0, 1.0)
    sf = hp.StrFctNLO(x, Q, Q, Q)
    f2, f1 = sf[hp.iF2EM], sf[hp.iF1EM]
    r = p.alphasQ(Q) / hp.AlphaS(Q)      # normalise to the LHAPDF alpha_s used here
    return r * f2 / x, r * (f2 - 2 * x * f1) / x


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pdf', default='NNPDF30_nlo_as_0118')
    ap.add_argument('--x', type=float, default=0.01)
    ap.add_argument('--Q', type=float, default=20.0)
    ap.add_argument('--y', type=float, default=0.5)
    ap.add_argument('--taus', default='1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,1e-6')
    ap.add_argument('--hoppet', action='store_true', help='also compare the exact NLO with hoppet')
    ap.add_argument('--cumulant-only', action='store_true', help='print only the tau_1^a and tau_1^b cumulants')
    ap.add_argument('--only', choices=['q', 'g'], help='one channel only (diagnostics; the Born then uses the quarks or nothing)')
    a = ap.parse_args()
    P = Pdfs(a.pdf, a.Q)
    P.only = a.only
    born = (1 + (1 - a.y) ** 2) * P.q(a.x)
    if a.cumulant_only:
        for t in map(float, a.taus.split(',')):
            print(f'tau {t:9.1e}  cumulant tau1a {cumulant(P, a.x, a.y, t, taua=True):.10e}  tau1b {cumulant(P, a.x, a.y, t):.10e}')
        return
    ex, f2, fl = exact_nlo(P, a.x, a.y)
    print(f'x = {a.x}, Q = {a.Q}, y = {a.y}, {a.pdf}, alpha_s(Q) = {P.as_:.6f}')
    print(f'Born R0 = {born:.8e};  exact O(as): {ex:.8e} ({ex / (born or ex):+.6f} of Born);  F2/x(1) = {f2:.6e}, FL/x(1) = {fl:.6e}')
    if a.hoppet:
        hf2, hfl = hoppet_nlo(a.pdf, a.x, a.Q)
        print(f'hoppet: F2/x(1) = {hf2:.6e} (ratio {hf2 / f2:.6f}), FL/x(1) = {hfl:.6e} (ratio {hfl / fl:.6f})')
    print(f'{"tau_cut":>9s} {"below":>14s} {"above":>14s} {"sum":>14s} {"(sum-exact)/exact":>18s} {"/Born":>12s} {"FL_above/FL-1":>14s}')
    for t in map(float, a.taus.split(',')):
        b = cumulant(P, a.x, a.y, t)
        r, rf2, rfl = real_above(P, a.x, a.y, t)
        s = b + r
        print(f'{t:9.1e} {b:14.6e} {r:14.6e} {s:14.6e} {(s - ex) / ex:18.3e} {(s - ex) / (born or ex):12.3e} {rfl / fl - 1:14.3e}')


if __name__ == '__main__':
    main()
