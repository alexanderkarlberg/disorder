#!/usr/bin/env python3
"""Combine nnlo11.dat files (runs with the same number of events, different
seeds) and print the O(alpha_s) and O(alpha_s^2) coefficients of the 1+1
observables, per Born (sigma_0) and in units of (alpha_s/2pi)^k, three ways:

  DISENT+P2B : C_k^incl O(Born) + E_k / (sigma_0 a^k)
  slicing+P2B: the same with the tau_2-sliced O(alpha_s^2) 2+1 part (E2s)
  tau_1 slicing: LP_k(tau_cut) O(Born) + A_k(tau_cut) / (sigma_0 a^k)

C_k^incl from nnlo11_ref.py (hoppet), LP_k from dis_tau1_lp (MCFM's SCET
pieces), both passed on the command line. Differences to DISENT+P2B use the
per-event paired sums written by nnlo11 (correlated errors)."""
import sys, glob, argparse
import numpy as np

NAMES = ['total', '>=1 jet', 'pt 5-10', 'pt 10-14', 'pt 14-15', 'pt 15-16', 'pt 16-17', 'pt 17-20', 'pt 20-40',
         'eta -1..-0.6', 'eta -0.6..-0.4', 'eta -0.4..-0.2', 'eta -0.2..0.2', 'eta 0.2..1', 'eta 1..2.5', '>=2 jets']


def load(f):
    w = [float(v.replace('D', 'E')) for v in open(f).read().split()]
    nevt, nb, nt1, nt2 = (int(v) for v in w[:4]); i = 4
    def take(n):
        nonlocal i
        a = np.array(w[i:i + n]); i += n; return a
    tc1 = take(nt1); tc2 = take(nt2); rel = take(nt2).astype(int)
    as2pi = take(1)[0]; fbB = take(nb)
    def arr(shape):            # Fortran column-major (2, ..., nb)
        n = int(np.prod(shape)); return take(n).reshape(shape[::-1]).T
    B0, E1, E2d = arr((2, nb)), arr((2, nb)), arr((2, nb))
    E2s, D2s = arr((2, nt2, nb)), arr((2, nt2, nb))
    A1, A2, Dp = arr((2, nt1, nb)), arr((2, nt1, nb)), arr((2, nt1, nb))
    return dict(nevt=nevt, nb=nb, tc1=tc1, tc2=tc2, rel=rel, as2pi=as2pi, fbB=fbB, B0=B0, E1=E1, E2d=E2d,
                E2s=E2s, D2s=D2s, A1=A1, A2=A2, Dp=Dp)


def comb(runs, key):
    """mean over runs of the integrals, error from the per-event variances"""
    I = np.array([r[key][0] for r in runs])
    V = np.array([r[key][1] - r[key][0] ** 2 / r['nevt'] for r in runs])
    return I.mean(0), np.sqrt(np.clip(V, 0, None).sum(0)) / len(runs)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+')
    ap.add_argument('--C1', type=float, required=True, help='inclusive O(alpha_s) coefficient per Born (hoppet)')
    ap.add_argument('--C2', type=float, required=True, help='inclusive O(alpha_s^2) coefficient per Born (hoppet)')
    ap.add_argument('--lp', required=True, help='dis_tau1_lp output (tau_cut LP1 LP2 lines)')
    ap.add_argument('--bins', default='all', help='comma-separated bin numbers (1-based) or all')
    a = ap.parse_args()
    files = sorted(sum((glob.glob(f) for f in a.files), []))
    runs = [load(f) for f in files]
    r0 = runs[0]
    lp = {}
    for line in open(a.lp):
        p = line.split()
        if len(p) == 3:
            try:
                lp[float(p[0])] = (float(p[1]), float(p[2]))
            except ValueError:
                pass
    def lpv(t, k):
        for tt, v in lp.items():
            if abs(tt / t - 1) < 1e-6:
                return v[k]
        raise KeyError(t)
    s = {k: comb(runs, k) for k in ('B0', 'E1', 'E2d', 'E2s', 'D2s', 'A1', 'A2', 'Dp')}
    B0 = s['B0'][0]; sig0 = B0[0]; a1 = r0['as2pi']; a2 = a1 ** 2
    O = B0 / sig0                      # O(Born) per bin
    bins = range(r0['nb']) if a.bins == 'all' else [int(b) - 1 for b in a.bins.split(',')]
    print(f'{len(runs)} runs x {r0["nevt"]} events; sigma_0 = {sig0:.6e}; alpha_s/2pi = {a1:.8f}')
    print(f'Born 1+1 bins: {r0["fbB"].astype(int).tolist()}')
    for b in bins:
        print(f'\n== bin {b + 1}: {NAMES[b]}   O(Born) = {O[b]:.0f}')
        c1d = a.C1 * O[b] + s['E1'][0][b] / (sig0 * a1); e1d = s['E1'][1][b] / (sig0 * a1)
        c2d = a.C2 * O[b] + s['E2d'][0][b] / (sig0 * a2); e2d = s['E2d'][1][b] / (sig0 * a2)
        print(f'  O(as)   DISENT+P2B {c1d: .5f} +- {e1d:.5f}')
        for it, t in enumerate(r0['tc1']):
            c = lpv(t, 0) * O[b] + s['A1'][0][it, b] / (sig0 * a1); e = s['A1'][1][it, b] / (sig0 * a1)
            print(f'          tau_1 slicing {t:7.0e}: {c: .5f} +- {e:.5f}   (- P2B {c - c1d: .5f})')
        print(f'  O(as^2) DISENT+P2B {c2d: .4f} +- {e2d:.4f}')
        for it, t in enumerate(r0['tc2']):
            c = a.C2 * O[b] + s['E2s'][0][it, b] / (sig0 * a2)
            dd = s['D2s'][0][it, b] / (sig0 * a2); ed = s['D2s'][1][it, b] / (sig0 * a2)
            lab = f'tau_2 > {t:.0e} tau_1' if r0['rel'][it] else f'tau_2 > {t:.0e}      '
            print(f'          slicing+P2B {lab}: {c: .4f}   - DISENT+P2B {dd: .4f} +- {ed:.4f}')
        for it, t in enumerate(r0['tc1']):
            c = lpv(t, 1) * O[b] + s['A2'][0][it, b] / (sig0 * a2)
            dd = (lpv(t, 1) - a.C2) * O[b] + s['Dp'][0][it, b] / (sig0 * a2); ed = s['Dp'][1][it, b] / (sig0 * a2)
            print(f'          tau_1 slicing {t:7.0e}: {c: .4f}   - DISENT+P2B {dd: .4f} +- {ed:.4f}')


if __name__ == '__main__':
    main()
