#!/usr/bin/env python3
"""G8 pilot analysis: G8 (flat) vs G8x4 vs J8 (psmc) per part and cell, per-seed stats; slope-error budget (item 4 of UPDATE 15 Oct).
  g8_pilot.py -> prints tables (stdout)"""
import glob, math, os, re, sys
import numpy as np
sys.argv = ['x', 'PH9']
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs/ours8j'
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def rows(f):
    tc, nm, r = None, [], []
    for l in open(f):
        w = l.split()
        if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
        elif l.startswith(' ZCELL '): nm.append('%s %g-%g' % (w[1], float(w[2]), float(w[3]))); r.append([ff(x) for x in w[4:]])
    return tc, nm, np.array(r)
def load(setn, part, maxn=None):
    out = []; seeds = []; tc = nm = None
    for d in sorted(glob.glob('%s/%s/%s/s*/' % (P, setn, part))):
        if os.path.exists(d + 'done'):
            tc, nm, a = rows(d + 'run.log'); out.append(a); seeds.append(int(re.search(r's(\d+)/$', d).group(1)))
    return tc, nm, np.array(out), seeds
def cpu(setn, part):
    v = []
    for d in sorted(glob.glob('%s/%s/%s/s*/' % (P, setn, part))):
        if os.path.exists(d + 'done'):
            u = s = None
            for l in open(d + 'time.log'):
                if 'User time' in l: u = float(l.split()[-1])
                if 'System time' in l: s = float(l.split()[-1])
            if u is not None: v.append((u + s) / 3600)
    return np.array(v)
def mad(x): return 1.4826 * np.median(abs(x - np.median(x)))
if __name__ == '__main__':
    cells = ['total 0-0', 'ptavg 8-15']
    for part in ('vi', 'kp', 'r'):
        tc, nm, J, _ = load('PH9', part); tc = np.array(tc)
        G = load('G8', part)[2]; G4 = load('G8x4', part)[2]
        print('=== %s: J8 psmc N=%d, G8 flat N=%d, G8x4 N=%d; CPU h/job: G8 %.2f, G8x4 %.2f, J8 %.2f' % (part, len(J), len(G), len(G4), cpu('G8', part).mean() if len(G) else 0, cpu('G8x4', part).mean() if len(G4) else float('nan'), cpu('PH9', part).mean()))
        for c in cells:
            i = nm.index(c); print(' -- cell', c)
            print('  %-8s %-26s %-26s %-26s %-9s %-9s' % ('tau_cut', 'J8 mean+-err [median]', 'G8 mean+-err [median]', 'G8x4 mean+-err [median]', 'pull G-J', 'pull G4-J'))
            if len(G4) == 0: G4 = np.full_like(J[:1], np.nan)
            for j, t in enumerate(tc):
                def st(A):
                    x = A[:, i, j]; return (x.mean(), x.std(ddof=1) / math.sqrt(len(x)), np.median(x)) if len(x) > 1 and np.isfinite(x).all() else (float('nan'),) * 3
                sj, sg, s4 = st(J), st(G), st(G4)
                def pl(a, b): return (a[0] - b[0]) / math.hypot(a[1], b[1])
                print('  %-8.0e %-26s %-26s %-26s %+7.1f   %+7.1f' % (t, '%.4g+-%.2g [%.3g]' % sj, '%.4g+-%.2g [%.3g]' % sg, '%.4g+-%.2g [%.3g]' % s4, pl(sg, sj), pl(s4, sj)))
