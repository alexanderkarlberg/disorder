#!/usr/bin/env python3
"""Z05: per part the covariance between the 10 tau_cut columns of every row (central scale), for the tau_cut fit at
MPP (every seed fills all columns, so the columns are correlated; summary_z05.py has the sd per cell only).

Usage: cov_z05.py <Z05 directory> [out.json.gz]
Reads the same files as summary_z05.py (one subdirectory per part, run logs *.out). Writes, per part: n, and per row
(1-31 of nlo31's zeus05 set) the 10 x 10 covariance of the seeds' sigma per bin (pb; tau_cut 2e-2 ... 1e-5) and their
Spearman rank correlation (robust against the heavy tails)."""
import sys, os, glob, json, gzip
import numpy as np
from summary_z05 import scells, PARTS, TCS


def ranks(x):
    """ranks of the seeds per column (n, 10); ties broken by order, as they are rare"""
    r = np.empty_like(x)
    for c in range(x.shape[1]):
        r[np.argsort(x[:, c]), c] = np.arange(len(x))
    return r


def main():
    D = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(D, 'Z05-cov.json.gz')
    res = {}
    for p in PARTS:
        a = [c[0] for c in (scells(f) for f in sorted(glob.glob(os.path.join(D, p, '*.out')))) if c is not None]
        if len(a) < 2:
            continue
        a = np.array(a)                                  # (n, 31, 10), central scale
        cov = [np.cov(a[:, r, :].T).tolist() for r in range(31)]
        spr = [np.nan_to_num(np.corrcoef(ranks(a[:, r, :]).T)).tolist() for r in range(31)]
        res[p] = dict(n=len(a), cov=cov, spearman=spr)
        c = np.array(cov[7]); d = np.sqrt(np.diag(c))
        print('%-4s n = %6d   row 8 (eta 2-2.5): corr(1e-3, 1e-4) = %6.3f, corr(2e-2, 1e-4) = %6.3f'
              % (p, len(a), c[4, 7]/max(d[4]*d[7], 1e-300), c[0, 7]/max(d[0]*d[7], 1e-300)))
    with gzip.open(out, 'wt') as f:
        json.dump(dict(tau_cut=TCS, rows='1-31 of the zeus05 set', scale='central', parts=res), f)
    print('wrote', out)


if __name__ == '__main__':
    main()
