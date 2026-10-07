#!/usr/bin/env python3
"""Item D: NNLO 1+1 (lab11 jet bins) by P2B on our tau_2-sliced NLO 2+1
against disorder -p2b -nnlocoef (DISENT). Equal-weight seed averages.

  combine_p2b11.py --b1 'glob [glob...]' --lo 'glob' --incl 'glob' --dis 'glob' [--lo0 'glob'] [--json f]

--b1/--lo: run.log of sliced21 b1 / nlo31 lo in mode 3 (LCELL rows);
--incl: disorder inclusive -nnlocoef histogram files; --dis: disorder -p2b
-nnlocoef histogram files; --lo0: disorder -lo files (LO jet rate for the
normalisation; per-mille statements are relative to LO >= 1 jet).
Prints per bin and tau_cut: ours, disorder, pull; chi2 over the 15 jet bins
(row 0, the total, is identically zero in P2B and is left out)."""
import argparse, glob, json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'slicing'))
from lab11_combine import load_disorder

PTE = [5, 8, 11, 15, 20, 30, 50, 100]
YE = [-1, -0.5, 0, 0.5, 1, 1.5, 2.5]
NAMES = ['total', '>=1 jet'] + ['pt %g-%g' % (a, b) for a, b in zip(PTE, PTE[1:])] + \
        ['y %g..%g' % (a, b) for a, b in zip(YE, YE[1:])] + ['>=2 jets']


def lcell(f):
    tc, r = None, []
    for l in open(f):
        if l.startswith(' tau_cut'):
            tc = [float(v) for v in l.split()[1:]]
        if l.startswith(' LCELL '):
            r.append([float(x) for x in l.split()[2:]])
    return (tc, np.array(r)) if len(r) == 16 else (None, None)


def files(pats):
    """space-separated globs; only job directories with a 'done' marker"""
    out = []
    for p in pats.split():
        out += [f for f in glob.glob(p) if os.path.exists(os.path.join(os.path.dirname(f), 'done'))]
    return sorted(set(out))


def me(a):
    a = np.array(a)
    return a.mean(0), a.std(0, ddof=1)/np.sqrt(len(a)), len(a)


def main():
    ap = argparse.ArgumentParser()
    for o in ('b1', 'lo', 'locorr', 'incl', 'dis', 'lo0', 'json'):
        ap.add_argument('--' + o)
    a = ap.parse_args()
    fb1 = [lcell(f) for f in files(a.b1)]
    flo = [lcell(f) for f in files(a.lo)]
    tc = [t for t, c in fb1 if t][0]
    B1 = me([c for t, c in fb1 if c is not None])
    if os.environ.get('TRIM_B1'):
        # diagnostic only (AK decides): 1% trimmed mean of b1 per cell (>= 1 seed per side)
        x = np.array([c for t, c in fb1 if c is not None]); n = len(x); k = max(1, int(round(0.01*n)))
        s = np.sort(x, axis=0)[k:n - k]
        print('b1 TRIMMED (diagnostic, %d of %d seeds per side removed per cell)' % (k, n))
        B1 = (s.mean(0), s.std(0, ddof=1)/np.sqrt(len(s)), len(s))
    LO3 = me([c for t, c in flo if c is not None])
    if a.locorr:
        # TECHDIFF correction of lo (default cut -> min(1e-10 W2, 1e-8 Q2)), added in quadrature
        LC = me([c for t, c in (lcell(f) for f in files(a.locorr)) if c is not None])
        print('lo TECHDIFF correction: %d seeds; largest |corr|/err %.1f' % (LC[2], np.max(np.abs(LC[0][1:])/np.maximum(LC[1][1:], 1e-300))))
        LO3 = (LO3[0] + LC[0], np.sqrt(LO3[1]**2 + LC[1]**2), LO3[2])
    I2 = me([load_disorder(f) for f in files(a.incl)])
    D2 = me([load_disorder(f) for f in files(a.dis)])
    L0 = me([load_disorder(f) for f in files(a.lo0)]) if a.lo0 else None
    jr = L0[0][1] if L0 else float('nan')
    print('b1 %d, lo %d seeds; disorder incl %d, p2b %d; LO >= 1 jet %.2f pb' % (B1[2], LO3[2], I2[2], D2[2], jr))
    ours = I2[0][:, None] + B1[0] + LO3[0]
    oe = np.sqrt(I2[1][:, None]**2 + B1[1]**2 + LO3[1]**2)
    out = {'tau_cut': tc, 'names': NAMES, 'ours': ours.tolist(), 'ours_err': oe.tolist(),
           'disorder': D2[0].tolist(), 'disorder_err': D2[1].tolist(), 'lo_jet_rate': jr,
           'b1': B1[0].tolist(), 'b1_err': B1[1].tolist(), 'lo': LO3[0].tolist(), 'lo_err': LO3[1].tolist()}
    print('%-10s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc))
    pulls = (ours - D2[0][:, None])/np.sqrt(oe**2 + D2[1][:, None]**2)
    for k, n in enumerate(NAMES):
        print('%-10s disorder %10.4f +- %.4f (%.2f permille of LO jet rate)' % (n, D2[0][k], D2[1][k], 1000*D2[1][k]/jr))
        print('%-10s' % '  ours' + ''.join('%11.4f' % v for v in ours[k]))
        print('%-10s' % '  +-' + ''.join('%11.4f' % v for v in oe[k]))
        print('%-10s' % '  pull' + ''.join('%11.1f' % v for v in pulls[k]))
    chi2 = (pulls[1:]**2).sum(0)
    print('%-10s' % 'chi2/15' + ''.join('%11.1f' % v for v in chi2))
    print('%-10s' % 'max err' + ''.join('%11.2f' % (1000*v/jr) for v in oe[1:].max(0)) + '   (permille of LO jet rate)')
    out['chi2'] = chi2.tolist()
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1)


main()
