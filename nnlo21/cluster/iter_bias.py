#!/usr/bin/env python3
"""Inverse-variance iteration weighting bias (instructions UPDATE 7 Oct):
per seed, equal-weight average of VEGAS iterations 2..6 (the ' iteration'
lines = the target cell: total at tau_cut 1e-5 in mode 2, all tau_zQ bins at
1e-5 in mode 1) against the reported cell. Plain equal-weight means over seeds.
  iter_bias.py label='globs' ..."""
import glob, os, sys, re
import numpy as np


def ff(x):
    try:
        return float(x)
    except ValueError:
        return float(re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', x))


def one(f):
    it, rep = [], None
    for l in open(f):
        w = l.split()
        if len(w) >= 3 and w[0] == 'iteration':
            it.append(ff(w[2]))
        elif l.startswith(' ZCELL total'):
            rep = ff(w[-1])
        elif l.startswith(' CELL ') and abs(ff(w[1]) - 0.05) < 1e-9 and abs(ff(w[2]) - 0.5) < 1e-9:
            rep = ff(w[-1])
    return np.array(it), rep


for arg in sys.argv[1:]:
    lab, pats = arg.split('=', 1)
    fs = sorted({f for p in pats.split() for f in glob.glob(p) if os.path.exists(os.path.join(os.path.dirname(f), 'done'))})
    eq, rp, mx = [], [], []
    for f in fs:
        it, rep = one(f)
        if len(it) >= 6 and rep is not None:
            eq.append(it[1:6].mean()); rp.append(rep); mx.append(np.max(np.abs(it)))
    eq, rp, mx = np.array(eq), np.array(rp), np.array(mx); d = eq - rp; n = len(d)
    fin = np.isfinite(d)
    se = lambda x: x.std(ddof=1)/np.sqrt(len(x))
    lim = 20*np.median(mx[np.isfinite(mx)])
    ok = np.isfinite(d) & (mx <= lim)
    k = max(1, int(round(0.01*n))); ds = np.sort(d[np.isfinite(d)])
    print('%-10s DIAG: median(eq - rep) %.2f; 1%%-trimmed mean %.2f; mean over %d unflagged seeds (all |it| <= 20 x median max|it| = %.3g): %.2f +- %.2f (eq %.2f, rep %.2f)' % (
          lab, np.median(d[np.isfinite(d)]), ds[k:len(ds)-k].mean(), ok.sum(), lim, d[ok].mean(), se(d[ok]), eq[ok].mean(), rp[ok].mean()))
    print('%-10s %5d seeds: equal-weight it2-6 %.6g +- %.3g | reported %.6g +- %.3g | eq - rep %.4g +- %.3g'
          ' | medians %.6g / %.6g | finite-only diff %.4g +- %.3g (%d)' % (
          lab, n, eq.mean(), se(eq), rp.mean(), se(rp), d.mean(), se(d), np.median(eq), np.median(rp),
          d[fin].mean(), se(d[fin]), fin.sum()))
