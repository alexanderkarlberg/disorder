#!/usr/bin/env python3
"""Outlier check before combining (instructions 1.B): for each set of job
directories (globs of run.log, only directories with 'done'), compare per bin
and tau_cut the mean, the median and the 1%-trimmed mean (at least one seed per side) of the seeds, in
units of the equal-weight error. Reads ZCELL (mode 2), LCELL (mode 3) or
CELL (mode 1) rows. Prints the worst bins and the most extreme seeds.
  check_parts.py 'glob [glob...]' [...]"""
import glob, os, sys
import numpy as np


def rows(f):
    r = []
    for l in open(f):
        w = l.split()
        if l.startswith(' ZCELL '):
            r.append([float(x) for x in w[4:]])
        elif l.startswith(' LCELL '):
            r.append([float(x) for x in w[2:]])
        elif l.startswith(' CELL '):
            r.append([float(x) for x in w[3:]])
    return np.array(r) if r else None


def trim(x, f=0.01):
    # at least one seed on each side
    c = max(1, int(round(f*len(x))))
    s = np.sort(x)
    return s[c:len(x) - c].mean()


for arg in sys.argv[1:]:
    fs = sorted({f for p in arg.split() for f in glob.glob(p) if os.path.exists(os.path.join(os.path.dirname(f), 'done'))})
    data = [(f, rows(f)) for f in fs]
    data = [(f, r) for f, r in data if r is not None]
    shapes = {r.shape for f, r in data}
    if len(shapes) != 1:
        print(arg, 'inconsistent shapes', shapes); continue
    a = np.array([r for f, r in data]); n = len(a)
    if n < 3:
        print('%s: %d seeds' % (arg, n)); continue
    m = a.mean(0); e = a.std(0, ddof=1)/np.sqrt(n); md = np.median(a, 0)
    t = np.apply_along_axis(trim, 0, a)
    with np.errstate(divide='ignore', invalid='ignore'):
        dm = np.nan_to_num(np.abs(m - md)/e); dt = np.nan_to_num(np.abs(m - t)/e)
    b, k = np.unravel_index(np.argmax(dt), dt.shape)
    z = (a[:, b, k] - md[b, k])/(1.4826*np.median(np.abs(a[:, b, k] - md[b, k])) + 1e-300)
    worst = np.argsort(-np.abs(z))[:3]
    flag = 'OK' if dt.max() < 2 and np.abs(z).max() < 50 else 'CHECK'
    print('%-5s %-45s %5d seeds: max|mean-trim|/err %.2f (row %d, tau %d), max|mean-median|/err %.2f; extreme seeds %s'
          % (flag, arg[:45], n, dt.max(), b, k, dm.max(),
             ', '.join('%s z=%.0f' % (os.path.basename(os.path.dirname(data[j][0])), z[j]) for j in worst)))
