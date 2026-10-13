#!/usr/bin/env python3
"""UPDATE 13 Oct (d) decision statistics: total row (CELL 0.05 0.50) of r-x runs at tau_cut 1e-3, 1e-4.
  pilot200.py   -> Ex7 (9101-9300), Ex6 (2000 seeds), Ex6 first 200 seeds (lowest seed numbers).
Plain equal-weight means. Outlier = |x - median| > 20 * (1.48 MAD)."""
import glob, os, re
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def load(d):
    v, tc = [], None
    for s in sorted(glob.glob(d + '/s*/'), key=lambda s: int(re.search(r's(\d+)/$', s).group(1))):
        if not os.path.exists(s + 'done'): continue
        for l in open(s + 'run.log'):
            w = l.split()
            if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
            elif l.startswith(' CELL ') and w[1] == '0.05' and w[2] == '0.50': v.append([ff(x) for x in w[3:]])
    return np.array(v), tc
def stats(a, tc, name):
    for j, t in enumerate(tc):
        if t not in (1e-3, 1e-4): continue
        x = a[:, j]; n = len(x); md = np.median(x); sd = x.std(ddof=1); m = 1.4826 * np.median(np.abs(x - md))
        yield dict(name=name, tau=t, n=n, mean=x.mean(), err=sd / np.sqrt(n), sd=sd, mad=m, med=md, frac=np.mean(np.abs(x - md) > 20 * m))
a7, tc = load(P + '/ours7/Ex7-pilot/r-x'); a6, tc6 = load(P + '/ours6/Ex6/r-x')
assert tc == tc6
rows = list(stats(a7, tc, 'Ex7 200')) + list(stats(a6, tc, 'Ex6 2000')) + list(stats(a6[:200], tc, 'Ex6 first 200'))
print('# total row; pb; plain equal-weight means; Ex7 N=%d finite=%d; Ex6 N=%d finite=%d' % (len(a7), np.isfinite(a7).all(1).sum(), len(a6), np.isfinite(a6).all(1).sum()))
print('%-14s %7s %5s %12s %10s %10s %9s %9s %9s' % ('set', 'tau_cut', 'N', 'mean', 'sd/sqrtN', 'sd', '1.48MAD', 'median', 'frac>20M'))
for r in sorted(rows, key=lambda r: (-r['tau'], r['name'])):
    print('%-14s %7g %5d %12.5g %10.4g %10.4g %9.4g %9.4g %9.4f' % (r['name'], r['tau'], r['n'], r['mean'], r['err'], r['sd'], r['mad'], r['med'], r['frac']))
print('# decision: (a) sd7^2*0.65 < sd6^2*1.43 ; (b) frac7 <= frac6 ; (c) |mean7-mean6| vs sqrt(err7^2+err6^2)')
for t in (1e-3, 1e-4):
    e7, e6 = [next(r for r in rows if r['name'] == n and r['tau'] == t) for n in ('Ex7 200', 'Ex6 2000')]
    print('tau_cut %g: (a) %.4g < %.4g : %s ; (b) %.4f <= %.4f : %s ; (c) diff %.4g, sigma %.4g, pull %.2f' % (
        t, e7['sd']**2 * 0.65, e6['sd']**2 * 1.43, e7['sd']**2 * 0.65 < e6['sd']**2 * 1.43, e7['frac'], e6['frac'], e7['frac'] <= e6['frac'],
        e7['mean'] - e6['mean'], np.hypot(e7['err'], e6['err']), (e7['mean'] - e6['mean']) / np.hypot(e7['err'], e6['err'])))
