#!/usr/bin/env python3
"""CH8 pilot: CPU per job per part (CH8q, CH8g, J8 PH9) and q + g vs J8 means (total row, tau_cut 1e-3..1e-5), first 5 seeds."""
import math, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g8_pilot import load, cpu, P
parts = ['b2', 'vi', 'kp', 'r', 'rcorr', 'd2', 'b1', 'lo']
njob = dict(b2=200, vi=240, kp=200, r=1600, rcorr=200, d2=100, b1=100, lo=100)
tot = 0
print('%-6s %8s %8s %8s %8s %8s' % ('part', 'J8 h', 'q h', 'g h', 'q/J8', 'g/J8'))
for p in parts:
    j, q, g = cpu('PH9', p).mean(), cpu('CH8q', p), cpu('CH8g', p)
    print('%-6s %8.3f %8.3f %8.3f %8.2f %8.2f  (N q %d g %d)' % (p, j, q.mean(), g.mean(), q.mean() / j, g.mean() / j, len(q), len(g)))
    tot += njob[p] * (q.mean() + g.mean())
print('estimated total core-h (q + g production, pilot CPU): %.0f' % tot)
for c in ('total 0-0', 'ptavg 8-15'):
    print('--- cell', c, ': per part, q+g (5 pilot seeds) vs J8 (all seeds); units pb/GeV2 (sliced), error = sd/sqrtN in quadrature')
    for p in parts:
        tc, nm, J, _ = load('PH9', p); Q = load('CH8q', p)[2]; G = load('CH8g', p)[2]
        i = nm.index(c)
        s = []
        for t in (1e-3, 1e-4, 1e-5):
            jx = list(np.isclose(tc, t)).index(True)
            e = lambda A: (A[:, i, jx].mean(), A[:, i, jx].std(ddof=1) / math.sqrt(len(A)))
            (a, ea), (b, eb), (m, em) = e(Q), e(G), e(J)
            s.append('%g: q %.3g g %.3g sum %.4g+-%.2g J8 %.4g+-%.2g pull %+.1f' % (t, a, b, a + b, math.hypot(ea, eb), m, em, (a + b - m) / math.hypot(ea, eb, em)))
        print(p); [print('   ', x) for x in s]
