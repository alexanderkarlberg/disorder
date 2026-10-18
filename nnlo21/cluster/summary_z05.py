#!/usr/bin/env python3
"""Z05 (N3LO 1+1 P2B in the ZEUS set-up of 1803.09973, LABSET=zeus05, SCALES=3): per part the seed statistics of
every cell, for the combination at MPP, and the N3LO P2B part (NNLO 2+1 above and below the cut) per row.

Usage: summary_z05.py <Z05 directory> [out.json.gz]
The directory holds one subdirectory per part (b2, vi, kp, r, rc1, rc2, b0, b1, lo) with the run logs *.out.
Writes, per part: n, and per (scale, row, tau_cut column) the mean, sd and 1.48 MAD over the seeds of the SCELL
cells (sigma per bin, pb; scales in NNLOJET's order (1,1) (2,2) (.5,.5) (2,1) (.5,1) (1,2) (1,.5); rows 1-31 of
nlo31's zeus05 set: >= 1 jet, eta_inc 8, E_T,inc 8, Q2 7, x 7; tau_cut 2e-2 ... 1e-5), and the fraction of seeds
beyond 20 MAD per part (total row, 1e-4)."""
import sys, os, glob, json, gzip
import numpy as np

SC = [(1.0, 1.0), (2.0, 2.0), (0.5, 0.5), (2.0, 1.0), (0.5, 1.0), (1.0, 2.0), (1.0, 0.5)]
TCS = [2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 3e-5, 1e-5]
PARTS = ['b2', 'vi', 'kp', 'r', 'rc1', 'rc2', 'b0', 'b1', 'lo']


def scells(path):
    d = {}
    with open(path, errors='replace') as f:
        for l in f:
            if l.startswith(' SCELL '):
                w = l.split()
                d[(float(w[1]), float(w[2]), int(w[3]))] = [float(v) for v in w[4:]]
    if len(d) != 7*31:
        return None
    return np.array([[d[(r, f, i)] for i in range(1, 32)] for r, f in SC])


def main():
    D = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(D, 'Z05-summary.json.gz')
    res = {}
    for p in PARTS:
        a = [c for c in (scells(f) for f in sorted(glob.glob(os.path.join(D, p, '*.out')))) if c is not None]
        if len(a) < 2:
            continue
        a = np.array(a)                                  # (n, 7, 31, 10)
        med = np.median(a, 0)
        mad = 1.48*np.median(np.abs(a - med), 0)
        x = a[:, 0, 0, 7]
        far = np.abs(x - med[0, 0, 7]) > 20*max(mad[0, 0, 7], 1e-300)
        res[p] = dict(n=len(a), mean=a.mean(0).tolist(), sd=a.std(0, ddof=1).tolist(), mad=mad.tolist(),
                      far20=float(far.mean()))
        print('%-4s n = %6d   seeds beyond 20 MAD (total row, 1e-4): %.4f' % (p, len(a), far.mean()))
    with gzip.open(out, 'wt') as f:
        json.dump(dict(scales=SC, tau_cut=TCS, parts=res), f)
    print('wrote', out)
    # the N3LO P2B part (central scale): b2 + vi + kp + r, plus rc1 + rc2 below tau_cut 1e-4; sigma per bin [pb]
    main4 = [p for p in ('b2', 'vi', 'kp', 'r') if p in res]
    if len(main4) < 4:
        return
    m = sum(np.array(res[p]['mean'])[0] for p in main4)
    v = sum((np.array(res[p]['sd'])[0])**2/res[p]['n'] for p in main4)
    for p in ('rc1', 'rc2'):
        if p in res:
            small = np.array(TCS) < 1e-4*(1 - 1e-9)
            m = m + np.array(res[p]['mean'])[0]*small
            v = v + (np.array(res[p]['sd'])[0])**2/res[p]['n']*small
    names = ['>=1 jet'] + ['eta %d' % i for i in range(1, 9)] + ['E_T %d' % i for i in range(1, 9)]
    print('\nN3LO P2B part (central scale), sigma per bin [pb], tau_cut', ' '.join('%9.0e' % t for t in TCS[4:]))
    for r in range(17):
        print('  %-8s' % names[r], ' '.join('%9.1f' % m[r, c] for c in range(4, 10)))
        print('  %-8s' % '  +-', ' '.join('%9.1f' % np.sqrt(v[r, c]) for c in range(4, 10)))


if __name__ == '__main__':
    main()
