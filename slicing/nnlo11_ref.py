#!/usr/bin/env python3
"""Exact inclusive coefficients for the NNLO 1+1 study (nnlo11.f90): photon
exchange, mu_R = mu_F = Q, from hoppet's structure functions with the LHAPDF
set assigned on the grid (no evolution). Prints C_k = [Y+ F2^(k) - y^2 FL^(k)] /
[Y+ F2^(0)] / (alpha_s/2pi)^k for k = 1, 2, with hoppet's own alpha_s(Q) in the
normalisation, for two grid spacings."""
import argparse, math, lhapdf, hoppet as hp

def coeffs(name, x, Q, y, dy):
    p = lhapdf.mkPDF(name, 0)
    ids = [-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6]
    hp.Start(dy, 3)
    hp.SetCoupling(p.alphasQ(91.1876), 91.1876, 3)
    hp.Assign(lambda xx, qq: [p.xfxQ(i, xx, qq) for i in ids])
    hp.StartStrFct(3)
    hp.InitStrFct(3, True, 1.0, 1.0)
    a = hp.AlphaS(Q) / (2 * math.pi)
    Yp = 1 + (1 - y) ** 2
    out = []
    for f in (hp.StrFctLO, hp.StrFctNLO, hp.StrFctNNLO):
        sf = f(x, Q, Q, Q)
        f2, f1 = sf[hp.iF2EM], sf[hp.iF1EM]
        fl = f2 - 2 * x * f1
        out.append(Yp * f2 - y * y * fl)
    s0 = out[0]
    return out[1] / s0 / a, out[2] / s0 / a ** 2, a, p.alphasQ(Q) / (2 * math.pi)

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--pdf', default='NNPDF30_nlo_as_0118')
ap.add_argument('--x', type=float, default=0.01)
ap.add_argument('--Q2', type=float, default=400.0)
ap.add_argument('--s', type=float, default=4 * 27.5 * 920.0)
a = ap.parse_args()
Q = math.sqrt(a.Q2); y = a.Q2 / (a.x * a.s)
for dy in (0.05, 0.025):
    c1, c2, ahp, alh = coeffs(a.pdf, a.x, Q, y, dy)
    print(f'x {a.x} Q2 {a.Q2} y {y:.6f} dy {dy}: C1 = {c1:.8f}  C2 = {c2:.6f}  (alpha_s/2pi hoppet {ahp:.8f}, LHAPDF {alh:.8f})')
