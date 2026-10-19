#!/usr/bin/env python3
"""Z05 wave 1 (UPDATE 18 Oct item 4): per part the finished jobs, CPU per job (user+sys from the /usr/bin/time -v report
at the end of each <part>_<seed>.out, core-h), the projected total for the full seed ranges of z05.py, and the means at
tau_cut 1e-4 (central scale: LCELL rows) of rows 1 (>= 1 jet), 8 (eta 2-2.5), 10 (E_T 6-8) against the MPP pilot.
  z05_wave1.py [Z05 dir]"""
import glob, math, os, re, sys
import numpy as np
Z = sys.argv[1] if len(sys.argv) > 1 else '/ptmp/mpp/akarlber/disorder-nnlo21/Z05'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
RANGE = {'r': (60001, 78700), 'vi': (80001, 82729), 'kp': (83001, 85479), 'b2': (86001, 86729), 'rc1': (87001, 88000),
         'rc2': (88001, 89000), 'b0': (89001, 89200), 'b1': (89201, 89600), 'lo': (89601, 90000)}
MPP = {'b2': [(73284, 41), (11133, 7), (37003, 39)], 'vi': [(-13213, 198), (-2848, 13), (-5443, 224)],
       'kp': [(6922, 39), (1241, 6), (3642, 47)], 'r': [(-63501, 2308), (-9640, 259), (-31578, 2605)]}
ROWS = [1, 8, 10]; J4 = 7   # tau_cut column of 1e-4
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def read(f):
    L, cpu = {}, None; u = s = None
    for l in open(f, errors='replace'):
        if l.startswith(' LCELL '):
            w = l.split(); L[int(w[1])] = [ff(x) for x in w[2:]]
        elif 'User time (seconds)' in l: u = float(l.split()[-1])
        elif 'System time (seconds)' in l: s = float(l.split()[-1])
    if u is not None: cpu = (u + (s or 0))/3600
    return L, cpu
tot = 0.0
print('%-4s %6s %8s %8s %10s %10s' % ('part', 'done', 'core-h', 'max', 'jobs', 'proj core-h'))
res = {}
for p, (a, b) in RANGE.items():
    fs = sorted(glob.glob('%s/%s/%s_*.out' % (Z, p, p)))
    D = [read(f) for f in fs]; c = [x[1] for x in D if x[1] is not None]
    n = b - a + 1; m = np.mean(c) if c else float('nan'); tot += m*n if c else 0
    print('%-4s %6d %8.3f %8.3f %10d %10.0f' % (p, len(fs), m, max(c) if c else float('nan'), n, m*n))
    res[p] = np.array([[x[0][r][J4] for r in ROWS] for x in D if all(r in x[0] for r in ROWS)])
print('projected total %.0f core-h (cap 50k)' % tot)
print('\nmeans at tau_cut 1e-4, central scale, sigma per bin [pb]: ours (n) vs MPP pilot, pull')
for p, ref in MPP.items():
    v = res[p]
    if len(v) < 2: print(p, 'too few'); continue
    for k, r in enumerate(ROWS):
        m = v[:, k].mean(); e = v[:, k].std(ddof=1)/math.sqrt(len(v)); M, E = ref[k]
        print('%-3s row %2d: %10.1f +- %7.1f (n %4d)  MPP %8.0f +- %5.0f  pull %+5.1f' % (p, r, m, e, len(v), M, E, (m - M)/math.hypot(e, E)))
