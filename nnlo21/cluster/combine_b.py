#!/usr/bin/env python3
"""Combine B (ZEUS, mode 2, ZCELL rows) or F (fixed point, mode 1, CELL rows)
with the TECHDIFF correction as a separate part (instructions 1.B, fefa1bc).

  combine_b.py [--nnlojet DIR] [--json F] label='glob glob' ...

labels: b0, b1, lo (NLO = b1 + lo), b2, vi, kp, r, rcorr (NNLO = b2 + vi + kp
+ r + rcorr), any other label is printed on its own (e.g. conv, edge).
Only job directories with 'done'. Equal weights per part (error = seed
scatter/sqrt(N)); parts add, errors in quadrature. Also prints, per part, the
largest |mean - trimmed mean| in units of the error (outlier check)."""
import glob, json, math, os, sys
import numpy as np
import re as _re


def ff(x):
    """Fortran reals, including 3-digit exponents written without E (1.0-154)"""
    try:
        return float(x)
    except ValueError:
        return float(_re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from combine_zeus import nnlojet


def rows(f):
    tc, names, r = None, [], []
    for l in open(f):
        w = l.split()
        if l.startswith(' tau_cut'):
            tc = [float(x) for x in w[1:]]
        elif l.startswith(' ZCELL '):
            names.append((w[1], float(w[2]), float(w[3]))); r.append([ff(x) for x in w[4:]])
        elif l.startswith(' CELL '):
            names.append(('tauzQ', float(w[1]), float(w[2]))); r.append([ff(x) for x in w[3:]])
    return tc, names, (np.array(r) if r else None)


def iters(f):
    v = []
    for l in open(f):
        w = l.split()
        if len(w) >= 3 and w[0] == 'iteration':
            v.append(ff(w[2]))
    return np.array(v) if v else np.zeros(1)


def trim(x):
    c = max(1, int(round(0.01*len(x)))); s = np.sort(x)
    return s[c:len(x) - c].mean()


def main():
    a = sys.argv[1:]; nd = js = None; maxit = {}; trimmed = False
    while a and a[0].startswith('--'):
        o = a.pop(0)
        if o == '--nnlojet': nd = a.pop(0)
        elif o == '--json': js = a.pop(0)
        elif o == '--maxit':
            # DIAGNOSTIC ONLY (AK: plain means are the result): leave out seeds
            # with non-finite cells or a VEGAS iteration |I| > threshold
            l, v = a.pop(0).split('='); maxit[l] = float(v)
        elif o == '--trimmed':
            trimmed = True   # DIAGNOSTIC ONLY: 1% trimmed means per cell
    comb, tc, names = {}, None, None
    for arg in a:
        lab, pats = arg.split('=', 1)
        fs = sorted({f for p in pats.split() for f in glob.glob(p) if os.path.exists(os.path.join(os.path.dirname(f), 'done'))})
        data = [rows(f) for f in fs]
        data = [d for d in data if d[2] is not None]
        if lab in maxit:
            keep = [i for i, f in enumerate(fs) if np.all(np.isfinite(data[i][2])) and
                    np.max(np.abs(iters(f))) <= maxit[lab]]
            print('%-6s DIAGNOSTIC: %d of %d seeds left out (non-finite or |iteration| > %g)' % (lab, len(data) - len(keep), len(data), maxit[lab]))
            data = [data[i] for i in keep]
        tc, names = data[0][0], data[0][1]
        x = np.array([d[2] for d in data]); n = len(x)
        m = x.mean(0); e = x.std(0, ddof=1)/math.sqrt(n)
        t = np.apply_along_axis(trim, 0, x) if n >= 3 else m
        with np.errstate(divide='ignore', invalid='ignore'):
            dt = np.nan_to_num(np.abs(m - t)/e)
        if trimmed and n >= 3:
            c = max(1, int(round(0.01*n))); s = np.sort(x, axis=0)[c:n - c]
            m, e = s.mean(0), s.std(0, ddof=1)/math.sqrt(len(s))
        comb[lab] = (m, e, n)
        print('%-6s %5d seeds; max |mean - trimmed|/err = %.2f' % (lab, n, dt.max()))
    ref = nnlojet(nd)
    out = {'tau_cut': tc, 'bins': {}}
    def tot(ps):
        if not all(p in comb for p in ps): return None
        return sum(comb[p][0] for p in ps), np.sqrt(sum(comb[p][1]**2 for p in ps))
    sums = [('LO', ['b0']), ('NLO', ['b1', 'lo']), ('NNLO', ['b2', 'vi', 'kp', 'r', 'rcorr']),
            ('NNLO(no corr)', ['b2', 'vi', 'kp', 'r'])]
    sums += [(l, [l]) for l in comb if l not in ('b0', 'b1', 'lo', 'b2', 'vi', 'kp', 'r')]
    print('tau_cut        ' + ''.join('%10.0e' % v for v in tc))
    for i, nm in enumerate(names):
        key = '%s %g-%g' % nm
        print('\n' + key); out['bins'][key] = {}
        for lab, ps in sums:
            v = tot(ps)
            if v is None: continue
            print('  %-13s' % lab + ''.join('%10.4g' % x for x in v[0][i]))
            print('  %-13s' % '  +-' + ''.join('%10.2g' % x for x in v[1][i]))
            out['bins'][key][lab] = [list(v[0][i]), list(v[1][i])]
            o = lab.split('(')[0]
            r = ref.get(o, {}).get(nm)
            if r and lab == o:
                print('  %-13s %.5g +- %.2g   pulls ' % ('NNLOJET', r[0], r[1]) +
                      ' '.join('%+.1f' % ((v[0][i][k] - r[0])/math.hypot(v[1][i][k], r[1])) for k in range(len(tc))))
                out['bins'][key]['nnlojet_' + o] = r
    if js:
        json.dump(out, open(js, 'w'), indent=1)


main()
