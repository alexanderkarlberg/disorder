#!/usr/bin/env python3
"""SVG plots of the tau_cut tests of the DIS 2+1 slicing (for the report pages).

  plot_tcut.py nlo-disent <tau2_combine output> out.svg
      (below + above - DISENT)/DISENT against tau_cut, per tau_zQ bin
  plot_tcut.py nnlo <combine_tcut json> out.svg [sliced21 b0 output]
      the NLO (b1 + lo, against DISENT) and NNLO (b2 + vi + kp + r)
      coefficients against tau_cut, per tau_zQ bin
Colours are placeholders replaced by CSS variables (--ink, --muted, --a1..--a6),
so the page can theme them."""
import sys, re, json
import numpy as np
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt

INK, MUTED = '#010101', '#020202'
ACC = ['#0a0a01', '#0a0a02', '#0a0a03', '#0a0a04', '#0a0a05', '#0a0a06']
CSS = {INK: 'var(--ink)', MUTED: 'var(--muted)'}
CSS.update({c: 'var(--a%d)' % (i + 1) for i, c in enumerate(ACC)})
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK, 'axes.labelcolor': INK, 'xtick.color': INK,
                     'ytick.color': INK, 'text.color': INK, 'svg.fonttype': 'none',
                     'font.family': 'sans-serif'})


def save(fig, fn):
    fig.savefig(fn, transparent=True, bbox_inches='tight')
    s = open(fn).read()
    for k, v in CSS.items():
        s = s.replace(k, v).replace(k.upper(), v)
    s = re.sub(r'<\?xml[^>]*>\s*', '', s)
    s = re.sub(r'<!DOCTYPE[^>]*>\s*', '', s)
    s = re.sub(r'width="[\d.]+pt" height="[\d.]+pt"', 'width="100%"', s, count=1)
    open(fn, 'w').write(s)


def nlo_disent(src, out):
    txt = open(src).read()
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    blocks = re.split(r'\ntau_zQ in ', txt)[1:]
    for k, blk in enumerate(blocks):
        lab = re.match(r'(\[[\d.]+,[\d.]+\))', blk).group(1)
        rows = [l.split() for l in blk.split('\n') if re.match(r'\s+\d\.\de-0\d', l)]
        t = np.array([float(r[0]) for r in rows]); y = np.array([float(r[4]) for r in rows])
        e = np.array([float(r[5]) for r in rows])
        sh = 1 + 0.05 * (k - 2.5)
        ax.errorbar(t * sh, y, e, fmt='o', ms=3, lw=1, color=ACC[k % 6], label='τ_zQ ' + lab)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xscale('log'); ax.set_xlabel('τ_cut'); ax.set_ylabel('(below + above − DISENT)/DISENT')
    ax.legend(frameon=False, fontsize=7, ncol=2)
    save(fig, out)


def born_cells(fn):
    b = {}
    for l in open(fn):
        if l.startswith(' CELL '):
            w = l.split(); b['%.2f-%.2f' % (float(w[1]), float(w[2]))] = float(w[3])
    return b


def nnlo(src, out, bornfile=None):
    d = json.load(open(src))
    born = born_cells(bornfile) if bornfile else None
    tc = np.array(d['tau_cut'])
    bins = list(d['bins'])
    fig, axs = plt.subplots(2, 1, figsize=(6.2, 5.6), sharex=True)
    for k, b in enumerate(bins):
        v = d['bins'][b]
        sh = 1 + 0.05 * (k - 2.5)
        if v['nlo'] and v['disent_nlo']:
            r, e = np.array(v['nlo'][0]), np.array(v['nlo'][1])
            ref = v['disent_nlo'][0]
            axs[0].errorbar(tc * sh, r / ref - 1, e / abs(ref), fmt='o', ms=3, lw=1, color=ACC[k % 6], label='τ_zQ ' + b)
        if v['nnlo']:
            r, e = np.array(v['nnlo'][0]), np.array(v['nnlo'][1])
            norm = born[b] if born else 1
            axs[1].errorbar(tc * sh, r / norm, e / norm, fmt='o', ms=3, lw=1, color=ACC[k % 6], label='τ_zQ ' + b)
    axs[0].axhline(0, color=MUTED, lw=0.8)
    axs[0].set_ylabel('NLO: (b1 + lo)/DISENT − 1')
    axs[1].set_ylabel('NNLO correction / LO' if born else 'NNLO correction')
    axs[1].axhline(0, color=MUTED, lw=0.8)
    if born:
        axs[1].set_ylim(-1.5, 6)
    axs[1].set_xscale('log'); axs[1].set_xlabel('τ_cut')
    axs[0].legend(frameon=False, fontsize=7, ncol=2)
    save(fig, out)


if __name__ == '__main__':
    {'nlo-disent': nlo_disent, 'nnlo': nnlo}[sys.argv[1]](*sys.argv[2:])
