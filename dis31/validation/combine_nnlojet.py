#!/usr/bin/env python3
"""Combine independent NNLOJET jobs (one directory per seed) for one channel.

For each job: the cross section (*.cross.s<seed>.dat, first two columns, fb)
and the q2 histogram (*.q2.s<seed>.dat: xlow xcen xhigh value error, value =
dsigma/dQ^2 in fb/GeV^2), summed over the sub-channels the job produced
(e.g. RRa + RRb). Jobs with equal statistics are combined with equal
weights, error = scatter/sqrt(N) (robust for heavy-tailed weights);
inverse-variance weighting is printed as a cross check. Output in pb, the
Q^2 bins as sigma per bin (as nlo31).
Usage: combine_nnlojet.py <dir with job directories> [glob of job dirs]"""
import glob, math, os, sys
import numpy as np

base = sys.argv[1]
pat = sys.argv[2] if len(sys.argv) > 2 else 's*'
tot, err, hists = [], [], []
for d in sorted(glob.glob(os.path.join(base, pat))):
    cr = sorted(glob.glob(os.path.join(d, '*.cross.s*.dat')))
    qh = sorted(glob.glob(os.path.join(d, '*.q2.s*.dat')))
    if not cr:
        continue
    t = e2 = 0
    for f in cr:
        row = [l.split() for l in open(f) if not l.startswith('#') and l.strip()][0]
        t += float(row[0]); e2 += float(row[1])**2
    h = None
    for f in qh:
        rows = np.array([[float(x) for x in l.split()[:5]] for l in open(f) if not l.startswith('#') and l.strip()])
        v = rows[:, 3]*(rows[:, 2] - rows[:, 0])/1000      # pb per bin
        h = v if h is None else h + v
    tot.append(t/1000); err.append(math.sqrt(e2)/1000); hists.append(h)
n = len(tot)
if n == 0:
    sys.exit('no finished jobs')
tot, err = np.array(tot), np.array(err)
mean = tot.mean(); sem = tot.std(ddof=1)/math.sqrt(n) if n > 1 else err[0]
w = 1/err**2
print('%d jobs: equal weights %.5f +- %.5f pb;  inverse variance %.5f +- %.5f pb'
      % (n, mean, sem, (tot*w).sum()/w.sum(), math.sqrt(1/w.sum())))
print('  median of the per-job errors %.4f pb (scatter/job %.4f pb)' % (np.median(err), tot.std(ddof=1) if n > 1 else 0))
if all(h is not None for h in hists):
    H = np.array(hists)
    print('  Q2 bins (pb):', ' '.join('%.4f(%.4f)' % (m, s) for m, s in
                                       zip(H.mean(axis=0), H.std(axis=0, ddof=1)/math.sqrt(n) if n > 1 else 0*H[0])))
