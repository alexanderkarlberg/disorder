#!/usr/bin/env python3
"""Z05 cost study (20 Oct): per part the heavy-tail statistics of the seeds, central scale, every row and tau_cut
column: the largest seed's share of sum (x - mean)^2; the error ratio err(n)/err(n/4) (err = sd/sqrt(n), n/4: the
median over 400 random subsets of n/4 seeds; 0.5 for well-behaved seeds, about 1 or more when a few giant weights
carry the variance); mean and sd of the column differences x(tau_cut) - x(1e-4). For a part p next to a part p + 's'
or p + 't' (the same seeds with ZSMEAR, e.g. rx and rxs): mean and sd of the per-seed difference (smoothed - sharp).

Usage: tails_z05.py <out.json.gz> <directory> ...
Every subdirectory (part) with >= 16 finished run logs <part>_<seed>.out is read (the LCELL rows; in mode 3 they
equal the central-scale SCELL rows). Writes per part: n and per (row 1-31 of the zeus05 set, tau_cut column) the
statistics above (sigma per bin, pb)."""
import sys, os, re, glob, json, gzip
import numpy as np
from summary_z05 import TCS

REF = 7                                                  # tau_cut 1e-4


def lcells(path):
    r = []
    with open(path, errors='replace') as f:
        for l in f:
            if l.startswith(' LCELL '):
                r.append([float(v) for v in l.split()[2:]])
    return np.array(r) if len(r) == 31 else None


def load(d):
    """{seed: (31, 10)} of the finished logs in d"""
    out = {}
    for f in glob.glob(os.path.join(d, '*.out')):
        m = re.search(r'_(\d+)\.out$', f)
        c = lcells(f) if m else None
        if c is not None:
            out[int(m.group(1))] = c
    return out


def stats(a, rng):
    n = len(a); d = a - a.mean(0); ss = (d**2).sum(0)
    share = (d**2).max(0)/np.where(ss > 0, ss, np.inf)
    m = n//4
    eq = np.median([a[rng.choice(n, m, replace=False)].std(0, ddof=1)/np.sqrt(m) for _ in range(400)], 0)
    err = a.std(0, ddof=1)/np.sqrt(n)
    dc = a - a[:, :, REF:REF + 1]
    return dict(n=n, mean=a.mean(0).tolist(), sd=a.std(0, ddof=1).tolist(), share=share.tolist(),
                ratio4=(err/np.where(eq > 0, eq, np.inf)).tolist(), dmean=dc.mean(0).tolist(), dsd=dc.std(0, ddof=1).tolist())


def main():
    out = sys.argv[1]; rng = np.random.default_rng(1)
    seeds, res = {}, {}
    for D in sys.argv[2:]:
        for d in sorted(glob.glob(os.path.join(D, '*/'))):
            p = os.path.basename(d.rstrip('/'))
            s = load(d)
            if len(s) < 16:
                continue
            seeds[p] = s
            res[p] = stats(np.array([s[k] for k in sorted(s)]), rng)
            sh = np.array(res[p]['share']); r4 = np.array(res[p]['ratio4'])
            # rows 5, 6 (eta 0.5-1, 1-1.5), 8 (eta 2-2.5), 17 (E_T 71-127), 20 (Q2 100-250); tau_cut 1e-3
            print('%-5s n = %6d  tau_cut 1e-3, rows 5 6 8 17 20: largest-seed share' % (p, len(s)),
                  ' '.join('%4.2f' % sh[i - 1, 4] for i in (5, 6, 8, 17, 20)), '  err ratio n/4 -> n',
                  ' '.join('%5.2f' % r4[i - 1, 4] for i in (5, 6, 8, 17, 20)))
    pairs = {}
    for p, q in [(p, p + x) for p in seeds for x in 'st' if p + x in seeds]:
        k = sorted(set(seeds[p]) & set(seeds[q]))
        if len(k) < 16:
            continue
        d = np.array([seeds[q][i] - seeds[p][i] for i in k])
        pairs[q] = dict(n=len(k), mean=d.mean(0).tolist(), sd=d.std(0, ddof=1).tolist())
        e = d.std(0, ddof=1)/np.sqrt(len(k)); e0 = np.array(res[p]['sd'])/np.sqrt(res[p]['n'])
        print('%-5s pairs %6d  tau_cut 1e-3, rows 5 6 8 17 20: (smoothed - sharp)/err(sharp)' % (q, len(k)),
              ' '.join('%6.2f +- %5.2f' % (d.mean(0)[i - 1, 4]/e0[i - 1, 4], e[i - 1, 4]/e0[i - 1, 4]) for i in (5, 6, 8, 17, 20)))
    with gzip.open(out, 'wt') as f:
        json.dump(dict(tau_cut=TCS, ref=TCS[REF], rows='1-31 of the zeus05 set', scale='central', parts=res, pairs=pairs), f)
    print('wrote', out)


if __name__ == '__main__':
    main()
