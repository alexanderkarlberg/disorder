#!/usr/bin/env python3
"""Combine nlo31 outputs: per part (lo, vi, kp, r) the seeds with equal
weights (runs with equal statistics), error = seed scatter/sqrt(N); then the
NLO correction vi + kp + r and LO + NLO. The real part has heavy tails and
VEGAS errors correlated with the values, which biases inverse-variance
weighting (3 Oct 2026: -29.37 against -29.75 +- 0.40 pb for 24 seeds), so
that is only printed as a cross check (and used when there is one seed).
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
    n = len(lst)
    w = [1/e**2 for _, e, _, _ in lst]
    tiv = sum(t*wi for (t, _, _, _), wi in zip(lst, w))/sum(w)
    eiv = math.sqrt(1/sum(w))
    nb = len(lst[0][2])
    if n > 1:
        tot = sum(t for t, _, _, _ in lst)/n
        err = math.sqrt(sum((t - tot)**2 for t, _, _, _ in lst)/(n - 1)/n)
        bins = [sum(b[2][i] for b in lst)/n for i in range(nb)]
        berr = [math.sqrt(sum((b[2][i] - bins[i])**2 for b in lst)/(n - 1)/n) for i in range(nb)]
    else:
        tot, err, bins, berr = tiv, eiv, lst[0][2], [0]*nb
    comb[part] = (tot, err, bins)
    print('%-3s %2d seeds: %14.6e +- %.3e   (inverse variance %.6e +- %.3e)' % (part, n, tot, err, tiv, eiv))
    print('     Q2 bins:', ' '.join('%.4e(%.1e)' % (b, e) for b, e in zip(bins, berr)))
if all(p in comb for p in ('vi', 'kp', 'r')):
    t = sum(comb[p][0] for p in ('vi', 'kp', 'r'))
    e = math.sqrt(sum(comb[p][1]**2 for p in ('vi', 'kp', 'r')))
    print('NLO correction (vi + kp + r): %.6e +- %.3e pb' % (t, e))
    if 'lo' in comb:
        print('LO: %.6e +- %.3e pb;  LO + NLO: %.6e +- %.3e pb;  K = %.4f' % (
            comb['lo'][0], comb['lo'][1], comb['lo'][0] + t, math.hypot(e, comb['lo'][1]),
            (comb['lo'][0] + t)/comb['lo'][0]))
