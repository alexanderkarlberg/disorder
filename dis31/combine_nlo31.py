#!/usr/bin/env python3
"""Combine nlo31 outputs: per part (lo, vi, kp, r) the seeds by inverse
variance (the VEGAS errors), then the NLO correction vi + kp + r and LO + NLO.
Usage: combine_nlo31.py file1.out file2.out ...
The Q^2 bins are combined with the same weights as the totals of each seed."""
import re, sys, collections, math

res = collections.defaultdict(list)
for fn in sys.argv[1:]:
    part = tot = err = None
    bins = []
    for line in open(fn):
        m = re.match(r'\s*RESULT (\S+) sigma\(>=3 jets\) \[pb\] =\s*(\S+)\s+\+-\s+(\S+)', line)
        if m:
            part, tot, err = m.group(1), float(m.group(2)), float(m.group(3))
        m = re.match(r'\s*Q2bin\s+(\S+)\s+(\S+)\s+(\S+)', line)
        if m:
            bins.append(float(m.group(3)))
    if part is None:
        print('no result in', fn); continue
    res[part].append((tot, err, bins, fn))

comb = {}
for part, lst in sorted(res.items()):
    w = [1/e**2 for _, e, _, _ in lst]
    tot = sum(t*wi for (t, _, _, _), wi in zip(lst, w))/sum(w)
    err = math.sqrt(1/sum(w))
    nb = len(lst[0][2])
    bins = [sum(b[2][i]*wi for b, wi in zip(lst, w))/sum(w) for i in range(nb)]
    # seed scatter as a cross check of the VEGAS errors
    sc = 0
    if len(lst) > 1:
        sc = math.sqrt(sum((t - tot)**2 for t, _, _, _ in lst)/(len(lst) - 1)/len(lst))
    comb[part] = (tot, err, bins)
    print('%-3s %2d seeds: %14.6e +- %.3e   (seed scatter %.3e)' % (part, len(lst), tot, err, sc))
    print('     Q2 bins:', ' '.join('%.4e' % b for b in bins))
if all(p in comb for p in ('vi', 'kp', 'r')):
    t = sum(comb[p][0] for p in ('vi', 'kp', 'r'))
    e = math.sqrt(sum(comb[p][1]**2 for p in ('vi', 'kp', 'r')))
    print('NLO correction (vi + kp + r): %.6e +- %.3e pb' % (t, e))
    if 'lo' in comb:
        print('LO: %.6e +- %.3e pb;  LO + NLO: %.6e +- %.3e pb;  K = %.4f' % (
            comb['lo'][0], comb['lo'][1], comb['lo'][0] + t, math.hypot(e, comb['lo'][1]),
            (comb['lo'][0] + t)/comb['lo'][0]))
