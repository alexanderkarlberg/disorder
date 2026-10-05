#!/usr/bin/env python3
"""Efficiency of P2B + tau_2 slicing against disorder for the lab11 bins:
CPU x error^2 ratio ours/disorder, independent b1 + lo (optimal split) and
balanced correlated sampling (sliced21 c1), thserv jobs only.
Usage: eff_plot.py RUNS out.svg [tau index, default 7 = 1e-4]"""
import sys, glob, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'nnlo21'))
from lab11_combine import NAMES, load_disorder
from p2b11_plot import lcell
import plot_tcut as pt


def jobs(globpat, test, loader):
    vals, times = [], []
    for d in glob.glob(globpat):
        if not os.path.exists(d + '/host') or not open(d + '/host').read().startswith('thserv'):
            continue
        o = d + '/out.txt'
        if not os.path.exists(o) or not test(open(o).read()):
            continue
        v = loader(d)
        if v is None:
            continue
        vals.append(v); times.append(os.path.getmtime(o) - os.path.getmtime(d + '/host'))
    return np.array(vals), np.mean(times)


def main():
    R, out = sys.argv[1], sys.argv[2]; it = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    cv, tcv = jobs(R + '/corrv/s*', lambda t: 'LCELL' in t, lambda d: lcell(d + '/out.txt'))
    b1, tb1 = jobs(R + '/p2b2/s*', lambda t: 'sliced21 part b1' in t, lambda d: lcell(d + '/out.txt'))
    lo, tlo = jobs(R + '/p2b2/s*', lambda t: 'nlo31 part lo' in t, lambda d: lcell(d + '/out.txt'))
    ds, tds = jobs(R + '/p2b2/s*', lambda t: 'TOTAL TIME' in t and 'nnlocoef' in t,
                   lambda d: (lambda f: load_disorder(f[0]) if f else None)(glob.glob(d + '/c2h_disorder*.dat')))
    ks = list(range(1, 16)); ind = []; cor = []
    for k in ks:
        vd = ds[:, k].var(ddof=1) * tds
        ind.append((np.sqrt(b1[:, k, it].var(ddof=1) * tb1) + np.sqrt(lo[:, k, it].var(ddof=1) * tlo))**2 / vd)
        cor.append(cv[:, k, it].var(ddof=1) * tcv / vd)
    fig, ax = pt.plt.subplots(figsize=(6.6, 3.6))
    x = np.arange(len(ks))
    ax.semilogy(x, ind, 'o', color=pt.MUTED, label='independent b1 + lo')
    ax.semilogy(x, cor, 's', color=pt.ACC[1], label='correlated (stratified, balanced VEGAS)')
    ax.axhline(1, color=pt.ACC[0], lw=1)
    ax.set_xticks(x); ax.set_xticklabels([NAMES[k].replace('>=', '≥') for k in ks], rotation=60, ha='right', fontsize=7)
    ax.set_ylabel('CPU × error²: ours / disorder')
    ax.legend(frameon=False, fontsize=7)
    pt.save(fig, out)


if __name__ == '__main__':
    main()
