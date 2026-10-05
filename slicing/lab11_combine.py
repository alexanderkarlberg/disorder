#!/usr/bin/env python3
"""NLO and NNLO 1+1 by pure tau_1 slicing (slicing/nnlo11.f90 -integrated,
nnlo11i.dat per seed) against disorder -p2b (-nlocoef, -nnlocoef, histograms
of analysis/lab11_analysis.f), lab-frame jet observables integrated over Q^2
and y. Units: sigma per bin in pb; the O(alpha_s^k) coefficient includes
alpha_s^k (LHAPDF's alpha_s in nnlo11, disorder's own coupling in disorder).

Seeds are combined with equal weights (error = scatter/sqrt(N)).
Usage: lab11_combine.py --slicing 's*/nnlo11i.dat' --c1 'c1_*.dat' --c2 'c2_*.dat' [--json f]"""
import glob, argparse, json
import numpy as np

NAMES = ['total', '>=1 jet', 'pt 5-8', 'pt 8-11', 'pt 11-15', 'pt 15-20', 'pt 20-30', 'pt 30-50', 'pt 50-100',
         'y -1..-0.5', 'y -0.5..0', 'y 0..0.5', 'y 0.5..1', 'y 1..1.5', 'y 1.5..2.5', '>=2 jets']


def load_slicing(fn):
    w = [float(v.replace('D', 'E')) for v in open(fn).read().split()]
    nevt, nb, nt = int(w[0]), int(w[1]), int(w[2]); i = 3
    def take(n):
        nonlocal i
        a = np.array(w[i:i + n]); i += n; return a
    tc = take(nt); take(8 + 7)
    def arr(shape):
        n = int(np.prod(shape)); return take(n).reshape(shape[::-1]).T
    B0 = arr((2, nb)); A1 = arr((2, nt, nb)); A2 = arr((2, nt, nb)); L1 = arr((2, nt, nb)); L2 = arr((2, nt, nb))
    return tc, B0[0], A1[0] + L1[0], A2[0] + L2[0]          # sums over events = sigma


def load_disorder(fn):
    vals, cur = [], None
    for l in open(fn):
        if l.startswith('# ') and 'index' in l:
            cur = l.split()[1]
        elif cur and l.strip() and not l.startswith('#'):
            a, b, v, e = map(float, l.split())
            vals.append(v * (b - a))
    return np.array(vals)


def comb(a):
    a = np.array(a); n = len(a)
    return a.mean(0), a.std(0, ddof=1) / np.sqrt(n), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slicing'); ap.add_argument('--c1'); ap.add_argument('--c2'); ap.add_argument('--json')
    a = ap.parse_args()
    S = [load_slicing(f) for f in sorted(glob.glob(a.slicing))]
    tc = S[0][0]
    lo = comb([s[1] for s in S]); n1 = comb([s[2] for s in S]); n2 = comb([s[3] for s in S])
    d1 = comb([load_disorder(f) for f in sorted(glob.glob(a.c1))]) if a.c1 else None
    d2 = comb([load_disorder(f) for f in sorted(glob.glob(a.c2))]) if a.c2 else None
    print('slicing %d seeds, disorder nlocoef %s, nnlocoef %s seeds' % (lo[2], d1 and d1[2], d2 and d2[2]))
    out = {'tau_cut': list(tc), 'bins': {}}
    for b, name in enumerate(NAMES):
        print('\n%-12s LO %.4g +- %.2g' % (name, lo[0][b], lo[1][b]))
        for k, (sl, ref) in enumerate(((n1, d1), (n2, d2)), 1):
            s = '  O(as^%d) ' % k
            if ref:
                s += 'disorder %10.4g +- %7.2g |' % (ref[0][b], ref[1][b])
            print(s)
            for it, t in enumerate(tc):
                v, e = sl[0][it, b], sl[1][it, b]
                pull = (v - ref[0][b]) / np.hypot(e, ref[1][b]) if ref else float('nan')
                print('     tau_cut %7.0e  slicing %10.4g +- %7.2g   pull %+5.1f' % (t, v, e, pull))
        out['bins'][name] = {'lo': [lo[0][b], lo[1][b]],
                             'nlo': [list(n1[0][:, b]), list(n1[1][:, b])], 'nnlo': [list(n2[0][:, b]), list(n2[1][:, b])],
                             'disorder_nlo': [d1[0][b], d1[1][b]] if d1 else None,
                             'disorder_nnlo': [d2[0][b], d2[1][b]] if d2 else None}
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1)


if __name__ == '__main__':
    main()
