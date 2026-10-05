#!/usr/bin/env python3
"""ZEUS-like dijet distributions (Q^2, ptavg_12, m12) at LO and NLO: tau_2
slicing (sliced21 b0; b1 + lo at the smallest tau_cut) against NNLOJET (LO, R + V).
Usage: plot_zeus_dist.py <runs dir> <nnlojet nlo dir> out.svg"""
import sys, glob, json, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_tcut as pt
from combine_zeus import read


def nnlojet(d, ch, obs):
    a = []
    for f in glob.glob('%s/s*_%s/DIS.zeus2j.%s.%s.s*.dat' % (d, ch, ch, obs)):
        rows = [[float(v) for v in l.split()] for l in open(f) if l.strip() and not l.startswith('#')]
        a.append([r[3] for r in rows]); edges = [r[0] for r in rows] + [rows[-1][2]]
    a = np.array(a) / 1000            # fb/unit -> pb/unit
    return np.array(edges), a.mean(0), a.std(0, ddof=1) / np.sqrt(len(a))


def ours(files, part_filter, it):
    vals = {}
    for f in files:
        r = read(f)
        if r is None or r[0] not in part_filter:
            continue
        for k, v in r[2].items():
            vals.setdefault((r[0], k), []).append(v[it])
    return vals


def main():
    R, N, out = sys.argv[1:4]
    lo = ours(glob.glob(R + '/zb0/b0_*.out'), ('b0',), -1)
    nl = ours(glob.glob(R + '/znlo2/s*/out.txt'), ('b1', 'lo'), -1)
    fig, axs = pt.plt.subplots(2, 3, figsize=(7.4, 4.6), gridspec_kw={'height_ratios': [1.6, 1]})
    for col, (obs, name, lab, logx) in enumerate((('q2', 'q2', 'Q² [GeV²]', True), ('ptavg', 'ptavg_12', 'p̄_T [GeV]', False),
                                                   ('m12', 'm12', 'm₁₂ [GeV]', False))):
        e, nlo_lo, nlo_loe = nnlojet(N, 'LO', name)
        _, r_, re_ = nnlojet(N, 'R', name); _, v_, ve_ = nnlojet(N, 'V', name)
        nj, nje = r_ + v_, np.hypot(re_, ve_)
        w = np.diff(e); c = (e[1:] + e[:-1]) / 2 if not logx else np.sqrt(e[1:] * e[:-1])
        keys = sorted([k for k in {k for (_, k) in lo} if k[0] == obs], key=lambda k: k[1])
        olo = np.array([np.mean(lo[('b0', k)]) for k in keys]) / w
        b1 = [np.array(nl[('b1', k)]) for k in keys]; l3 = [np.array(nl[('lo', k)]) for k in keys]
        onl = np.array([b.mean() + l.mean() for b, l in zip(b1, l3)]) / w
        onle = np.array([np.hypot(b.std(ddof=1) / np.sqrt(len(b)), l.std(ddof=1) / np.sqrt(len(l))) for b, l in zip(b1, l3)]) / w
        ax = axs[0, col]
        ax.stairs(nlo_lo, e, color=pt.MUTED, lw=1.2, label='LO, NNLOJET')
        ax.plot(c, olo, 'o', ms=3, color=pt.ACC[0], label='LO, ours')
        ax.stairs(nj, e, color=pt.ACC[2], lw=1.2, label='NLO coefficient, NNLOJET')
        ax.errorbar(c, onl, onle, fmt='s', ms=3, color=pt.ACC[1], label='NLO coefficient, ours (τ_cut 1e-5)')
        ax.axhline(0, color=pt.MUTED, lw=0.6)
        if logx:
            ax.set_xscale('log'); ax.set_yscale('symlog', linthresh=1e-3)
        ax.set_title(lab, fontsize=8)
        if col == 0:
            ax.set_ylabel('dσ/dX [pb per unit]'); ax.legend(frameon=False, fontsize=6)
        ax = axs[1, col]
        ax.errorbar(c, (onl - nj) / nlo_lo * 100, np.hypot(onle, nje) / nlo_lo * 100, xerr=w / 2, fmt='s', ms=3, color=pt.ACC[1])
        ax.errorbar(c, (olo - nlo_lo) / nlo_lo * 100, nlo_loe / nlo_lo * 100, fmt='o', ms=3, color=pt.ACC[0])
        ax.axhline(0, color=pt.MUTED, lw=0.8)
        if logx:
            ax.set_xscale('log')
        if col == 0:
            ax.set_ylabel('(ours − NNLOJET)/LO [%]')
        ax.set_xlabel(lab)
    pt.save(fig, out)


if __name__ == '__main__':
    main()
