#!/usr/bin/env python3
"""Outlier report (instructions UPDATE 7 Oct, item 4): per part and tau_cut
column (the total row of ZCELL/CELL output), the plain mean (= the result,
AK), the median and the 1% trimmed mean (diagnostics), and the seeds with any
VEGAS iteration beyond 20x the median |iteration| of that part.
  outlier_report.py label='globs' ..."""
import glob, os, re, sys
import numpy as np


def ff(x):
    try:
        return float(x)
    except ValueError:
        return float(re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', x))


def load(f):
    it, tot, tc = [], None, None
    for l in open(f):
        w = l.split()
        if len(w) >= 3 and w[0] == 'iteration':
            it.append(abs(ff(w[2])))
        elif l.startswith(' tau_cut'):
            tc = [ff(x) for x in w[1:]]
        elif l.startswith(' ZCELL total'):
            tot = [ff(x) for x in w[4:]]
        elif l.startswith(' CELL ') and abs(ff(w[1]) - 0.05) < 1e-9 and abs(ff(w[2]) - 0.5) < 1e-9:
            tot = [ff(x) for x in w[3:]]
    return np.array(it), (np.array(tot) if tot else None), tc


for arg in sys.argv[1:]:
    lab, pats = arg.split('=', 1)
    fs = sorted({f for p in pats.split() for f in glob.glob(p) if os.path.exists(os.path.join(os.path.dirname(f), 'done'))})
    D = [(f,) + load(f) for f in fs]
    D = [d for d in D if d[2] is not None]
    x = np.array([d[2] for d in D]); n = len(x); tc = D[0][3]
    with np.errstate(all='ignore'):
        mean = x.mean(0); err = x.std(0, ddof=1)/np.sqrt(n)
    med = np.median(x, 0); k = max(1, int(round(0.01*n)))
    trim = np.sort(x, 0)[k:n - k].mean(0)
    mi = np.median([np.max(d[1]) for d in D if len(d[1])])
    flag = [os.path.basename(os.path.dirname(d[0])) for d in D if len(d[1]) and (np.max(d[1]) > 20*mi or not np.all(np.isfinite(d[1])))]
    print('== %s: %d seeds; median max|iteration| %.4g; flagged (any |iteration| > 20x): %d' % (lab, n, mi, len(flag)))
    print('   tau_cut ' + ''.join('%12.0e' % t for t in tc))
    print('   mean    ' + ''.join('%12.5g' % v for v in mean) + '   (result)')
    print('    +-     ' + ''.join('%12.3g' % v for v in err))
    print('   median  ' + ''.join('%12.5g' % v for v in med))
    print('   trim1%   ' + ''.join('%12.5g' % v for v in trim))
    print('   flagged: ' + ' '.join(flag))
