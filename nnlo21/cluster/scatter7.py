#!/usr/bin/env python3
"""Per-seed scatter (sd, 1.48 MAD) of the total row (CELL 0.05 0.50) of r-x runs: scatter7.py DIR [DIR...]  (DIR = .../<set>/r-x)"""
import glob, os, sys, re
import numpy as np
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
for d in sys.argv[1:]:
    v, tc = [], None
    for s in sorted(glob.glob(d + '/s*/')):
        if not os.path.exists(s + 'done'): continue
        for l in open(s + 'run.log'):
            w = l.split()
            if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
            elif l.startswith(' CELL ') and w[1] == '0.05' and w[2] == '0.50': v.append([ff(x) for x in w[3:]])
    a = np.array(v); print(d, 'N=%d finite=%d' % (len(a), np.isfinite(a).all(1).sum()))
    for j, t in enumerate(tc):
        if t in (1e-4, 1e-3): x = a[:, j]; md = np.median(x); print('  tau_cut %g  mean %.4g  sd %.4g  1.48MAD %.4g  median %.4g' % (t, x.mean(), x.std(ddof=1), 1.4826*np.median(np.abs(x - md)), md))
