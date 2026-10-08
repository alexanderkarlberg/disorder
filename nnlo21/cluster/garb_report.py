#!/usr/bin/env python3
"""UPDATE 8 Oct, step 3: per set, before (bin2/old) against after (bin3,
DIPGARB): flagged seeds (any |iteration| > 20x the median max|iteration|),
per-seed garbage-drop counts, and for the total row at selected tau_cut the
plain mean (result), median, 1% trimmed mean, robust (1.4826 MAD) and sample
per-seed sigma.
  garb_report.py label=before_dir:after_dir ..."""
import glob, os, re, sys
import numpy as np

TCI = [(3, '2e-3'), (4, '1e-3'), (6, '2e-4'), (7, '1e-4'), (8, '3e-5'), (9, '1e-5')]


def ff(x):
    try:
        return float(x)
    except ValueError:
        return float(re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', x))


def load(d):
    out = []
    for f in sorted(glob.glob(d + '/s*/run.log')):
        if not os.path.exists(os.path.dirname(f) + '/done'):
            continue
        it, tot, drop = [], None, None
        for l in open(f):
            w = l.split()
            if len(w) >= 3 and w[0] == 'iteration':
                it.append(abs(ff(w[2])))
            elif l.startswith(' ZCELL total'):
                tot = [ff(x) for x in w[4:]]
            elif l.startswith(' CELL ') and w[1] == '0.05' and w[2] == '0.50':
                tot = [ff(x) for x in w[3:]]
            elif 'events dropped (garbage dipole)' in l:
                drop = (int(w[-3]), int(w[-1]))
        if tot is not None and it:
            out.append((os.path.basename(os.path.dirname(f)), np.array(it), np.array(tot), drop))
    return out


def stats(D):
    mx = np.array([d[1].max() for d in D]); lim = 20*np.median(mx[np.isfinite(mx)])
    X = np.array([d[2] for d in D])
    flag = [d[0] for d, m in zip(D, mx) if not (m <= lim) or not np.all(np.isfinite(d[2]))]
    n = len(X); k = max(1, int(round(0.01*n)))
    rows = []
    for j, t in TCI:
        v = X[:, j]; fin = v[np.isfinite(v)]
        with np.errstate(all='ignore'):
            mean = v.mean(); err = v.std(ddof=1)/np.sqrt(n); ss = v.std(ddof=1)
        md = np.median(fin); mad = 1.4826*np.median(np.abs(fin - md)); tr = np.sort(fin)[k:len(fin) - k].mean()
        rows.append((t, mean, err, md, tr, mad, ss))
    return n, lim, flag, rows


for arg in sys.argv[1:]:
    lab, dirs = arg.split('=', 1)
    bd, ad = dirs.split(':')
    print('=' * 100); print(lab)
    for name, d in (('before', bd), ('after ', ad)):
        D = load(d)
        if len(D) < 3:
            print('  %s: %d seeds' % (name, len(D))); continue
        n, lim, flag, rows = stats(D)
        dr = [x[3] for x in D if x[3]]
        dtxt = ''
        if dr:
            c = np.array([a for a, b in dr]); tot = np.array([b for a, b in dr])
            dtxt = '; garbage-dipole drops per seed: mean %.1f, median %.0f, max %d (fraction %.1e)' % (
                c.mean(), np.median(c), c.max(), c.sum()/tot.sum())
        print('  %s: %d seeds, flagged %d (20x rule, limit %.3g)%s' % (name, n, len(flag), lim, dtxt))
        if flag:
            print('     flagged: ' + ' '.join(flag[:40]) + (' ...' if len(flag) > 40 else ''))
        for t, mean, err, md, tr, mad, ss in rows:
            ms = '%.6g ± %.3g' % (mean, err) if np.isfinite(mean) and abs(mean) < 1e8 else 'undefined'
            print('     tau %-5s plain %-22s median %-10.6g trim %-10.6g robust sigma %-9.4g sample sigma %.4g' % (
                t, ms, md, tr, mad, ss if np.isfinite(ss) and ss < 1e8 else float('inf')))
