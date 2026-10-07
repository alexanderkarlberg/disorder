#!/usr/bin/env python3
"""Size the B production (NNLO coefficient b2 + vi + kp + r) from pilot runs:
per-seed sigma per bin and tau_cut, mean CPU per job; optimal allocation
N_i ∝ sigma_i/sqrt(t_i) for a target error on the total.
  size_b.py TARGET_PB 'b2glob' 'viglob' 'kpglob' 'rglob'"""
import sys, glob, re, os, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from combine_zeus import read

tgt = float(sys.argv[1])
S, T, keys, tc = {}, {}, None, None
for name, g in zip(('b2', 'vi', 'kp', 'r'), sys.argv[2:6]):
    fs = [f for f in glob.glob(g) if os.path.exists(os.path.join(os.path.dirname(f), 'done'))]
    runs = [read(f) for f in fs]
    tc = runs[0][1]; keys = list(runs[0][2])
    a = np.array([[r[2][k] for k in keys] for r in runs])
    S[name] = a.std(0, ddof=1)
    cp = []
    for f in fs:
        t = open(os.path.join(os.path.dirname(f), 'time.log')).read()
        cp.append(float(re.search(r'User time \(seconds\): ([\d.]+)', t).group(1)))
    T[name] = np.mean(cp)
    print('%-3s %3d seeds, CPU/job %.0f s; per-seed sigma of the total at tau_cut %s: %s' % (
        name, len(fs), T[name], tc, ' '.join('%.1f' % v for v in S[name][0])))
for k in (7, 8):
    i = 0
    w = sum(S[p][i, k]*math.sqrt(T[p]) for p in S)
    N = {p: S[p][i, k]*w/math.sqrt(T[p])/tgt**2 for p in S}
    cost = sum(N[p]*T[p] for p in S)/3600
    print('tau_cut %.0e, total to +-%.2f pb: %s; %.0f core-h' % (tc[k], tgt, ' '.join('%s %.0f' % (p, N[p]) for p in N), cost))
