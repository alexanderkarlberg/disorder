#!/usr/bin/env python3
"""Lab-frame 1+1 jet distributions (leading-jet p_T and rapidity, lab11 bins)
from P2B on our tau_2-sliced 2+1 against disorder -p2b, at NLO and (when the
runs exist) NNLO. Writes an SVG for the report page (CSS colour variables as
nnlo21/plot_tcut.py) and prints the numbers.

Ours:  NLO  = disorder inclusive -nlocoef + sliced21 b0 (mode 3, P2B)
       NNLO = disorder inclusive -nnlocoef + sliced21 b1 + nlo31 lo (mode 3)
Reference: disorder -p2b -nlocoef / -nnlocoef (DISENT cutoff 1e-10).
Differences in % of the LO cross section with >= 1 jet.
LO per bin: DISENT's 1+1 Born (nnlo11 -integrated, nnlo11i.dat).
Usage: p2b11_plot.py RUNS_DIR out.svg [tau_cut indices for NNLO, default 5,7,9]"""
import sys, glob, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'nnlo21'))
from lab11_combine import load_disorder, load_slicing
import plot_tcut as pt

PTE = [5, 8, 11, 15, 20, 30, 50, 100]
YE = [-1, -0.5, 0, 0.5, 1, 1.5, 2.5]
TCS = [2e-2, 1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 3e-5, 1e-5]


def lcell(f):
    r = [l.split() for l in open(f) if l.startswith(' LCELL ')]
    return np.array([[float(x) for x in w[2:]] for w in r]) if len(r) == 16 else None


def me(a):
    a = np.array(a)
    return a.mean(0), a.std(0, ddof=1) / np.sqrt(len(a)), len(a)


def main():
    R = sys.argv[1]; out = sys.argv[2]
    its = [int(v) for v in sys.argv[3].split(',')] if len(sys.argv) > 3 else [5, 7, 9]
    lo = me([load_slicing(f)[1] for f in glob.glob(R + '/lab11/s*/nnlo11i.dat')])
    # runs of different size are not mixed (equal-weight seed averages): the
    # high-statistics sets replace the first ones once >= 8 of them exist
    def pick(lo_stat, hi_stat):
        return hi_stat if len(hi_stat) >= 8 else lo_stat
    fb0 = pick(glob.glob(R + '/p2b1/b0_*.out'),
               [f for f in glob.glob(R + '/p2b2/s*/out.txt') if 'sliced21 part b0' in open(f).read() and lcell(f) is not None])
    b0 = me([c[:, 0] for c in map(lcell, fb0) if c is not None])
    i1 = me([load_disorder(f) for f in pick(glob.glob(R + '/p2b1/i1_disorder*.dat'), glob.glob(R + '/p2b2/s*/i1h_disorder*.dat'))])
    d1 = me([load_disorder(f) for f in pick(glob.glob(R + '/lab11c/s*/c1x_disorder*.dat'), glob.glob(R + '/p2b2/s*/c1h_disorder*.dat'))])
    nlo = (i1[0] + b0[0], np.hypot(i1[1], b0[1]))
    nn = None
    fb1 = [f for f in glob.glob(R + '/p2b2/s*/out.txt') if 'sliced21 part b1' in open(f).read() and lcell(f) is not None]
    flo = [f for f in glob.glob(R + '/p2b2/s*/out.txt') if 'nlo31 part lo' in open(f).read() and lcell(f) is not None]
    fi2 = glob.glob(R + '/p2b2/s*/i2_disorder*.dat') + glob.glob(R + '/p2b1/i2_disorder*.dat')
    fd2 = glob.glob(R + '/lab11b/s*/c2x_disorder*.dat')
    if len(fb1) > 1 and len(flo) > 1 and len(fi2) > 1:
        B1 = me([lcell(f) for f in fb1]); LO3 = me([lcell(f) for f in flo]); I2 = me([load_disorder(f) for f in fi2])
        D2 = me([load_disorder(f) for f in fd2])
        nn = (I2[0][:, None] + B1[0] + LO3[0], np.sqrt(I2[1][:, None]**2 + B1[1]**2 + LO3[1]**2), D2, len(fb1), len(flo))
    print('LO %d seeds, b0 %d, disorder incl NLO %d, disorder p2b NLO %d' % (lo[2], b0[2], i1[2], d1[2]))
    for k, n in enumerate(['total', '>=1 jet'] + ['pt %g-%g' % (a, b) for a, b in zip(PTE, PTE[1:])] +
                          ['y %g..%g' % (a, b) for a, b in zip(YE, YE[1:])] + ['>=2 jets']):
        o, oe = nlo[0][k], nlo[1][k]
        print('  NLO %-12s ours %10.4f +- %7.4f  disorder %10.4f +- %7.4f  pull %+5.1f' % (n, o, oe, d1[0][k], d1[1][k],
              (o - d1[0][k]) / np.hypot(oe, d1[1][k])))
    if nn:
        print('NNLO: b1 %d, lo %d seeds' % (nn[3], nn[4]))
    fig, axs = pt.plt.subplots(3 if nn else 2, 2, figsize=(7.0, 7.6 if nn else 5.4), sharex='col',
                               gridspec_kw={'height_ratios': [1.6, 1, 1] if nn else [1.6, 1]})
    ljet = lo[0][1]     # LO cross section with >= 1 jet: normalisation of the differences
    for col, (sl, edges, lab) in enumerate(((slice(2, 9), PTE, 'leading-jet p_T [GeV]'), (slice(9, 15), YE, 'leading-jet y (lab)'))):
        e = np.array(edges, float); w = np.diff(e); c = (e[1:] + e[:-1]) / 2
        ax = axs[0, col]
        ax.stairs(lo[0][sl] / w, e, color=pt.MUTED, lw=1.2, label='LO')
        ax.stairs(d1[0][sl] / w, e, color=pt.ACC[0], lw=1.2, label='NLO coefficient, disorder')
        ax.errorbar(c, nlo[0][sl] / w, nlo[1][sl] / w, fmt='o', ms=3, color=pt.ACC[1], label='NLO coefficient, P2B + ours')
        if nn:
            ax.stairs(nn[2][0][sl] / w, e, color=pt.ACC[2], lw=1.2, label='NNLO coefficient, disorder')
            ax.errorbar(c, nn[0][sl, its[-1]] / w, nn[1][sl, its[-1]] / w, fmt='s', ms=3, color=pt.ACC[3],
                        label='NNLO coefficient, P2B + τ₂ slicing (%.0e)' % TCS[its[-1]])
        ax.axhline(0, color=pt.MUTED, lw=0.6)
        if col == 0:
            ax.set_xscale('log'); ax.set_xlim(4.5, 110); ax.set_ylabel('dσ/dX [pb per unit]')
        ax.legend(frameon=False, fontsize=6)
        ax = axs[1, col]
        r = (nlo[0][sl] - d1[0][sl]) / ljet * 100; re = np.hypot(nlo[1][sl], d1[1][sl]) / ljet * 100
        ax.errorbar(c, r, re, xerr=w / 2, fmt='o', ms=3, color=pt.ACC[1])
        ax.axhline(0, color=pt.MUTED, lw=0.8)
        if col == 0:
            ax.set_ylabel('NLO: (ours − disorder)\n/ σ_LO(≥1 jet) [%]')
        if nn:
            ax = axs[2, col]
            for k, it in enumerate(its):
                r = (nn[0][sl, it] - nn[2][0][sl]) / ljet * 100
                re = np.hypot(nn[1][sl, it], nn[2][1][sl]) / ljet * 100
                ax.errorbar(c * (1 + 0.02 * (k - 1)) if col == 0 else c + 0.03 * (k - 1), r, re, fmt='o', ms=3,
                            color=pt.ACC[(k + 2) % 6], label='τ_cut %.0e' % TCS[it])
            ax.axhline(0, color=pt.MUTED, lw=0.8)
            if col == 0:
                ax.set_ylabel('NNLO: (ours − disorder)\n/ σ_LO(≥1 jet) [%]')
            ax.legend(frameon=False, fontsize=6)
        axs[-1, col].set_xlabel(lab)
    pt.save(fig, out)


if __name__ == '__main__':
    main()
