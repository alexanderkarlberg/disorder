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


def compare(old, new, bornfile, out, labels='uniform r (4 Oct),psmc + tau2 fix (5 Oct)'):
    lab = labels.split(',')
    born = born_cells(bornfile)['0.05-0.50']
    fig, axs = plt.subplots(2, 1, figsize=(6.2, 5.6), sharex=True, gridspec_kw={'height_ratios': [1, 1.6]})
    for k, (src, l) in enumerate(zip((old, new), lab)):
        d = json.load(open(src)); tc = np.array(d['tau_cut']); v = d['bins']['0.05-0.50']
        sh = 1 + 0.06*(k - 0.5)
        if v['nlo'] and v['disent_nlo']:
            r, e = np.array(v['nlo'][0]), np.array(v['nlo'][1]); ref = v['disent_nlo'][0]
            axs[0].errorbar(tc*sh, r/ref - 1, e/abs(ref), fmt='o', ms=3.5, lw=1, color=ACC[k], label=l)
        r, e = np.array(v['nnlo'][0]), np.array(v['nnlo'][1])
        axs[1].errorbar(tc*sh, r/born, e/born, fmt='o', ms=3.5, lw=1, color=ACC[k], label=l)
    axs[0].axhline(0, color=MUTED, lw=0.8); axs[1].axhline(0, color=MUTED, lw=0.8)
    axs[0].set_ylabel('NLO: (b1 + lo)/DISENT − 1'); axs[0].set_ylim(-0.7, 0.15)
    axs[1].set_ylabel('NNLO coefficient / LO'); axs[1].set_ylim(-0.45, 0.6)
    axs[1].set_xscale('log'); axs[1].set_xlabel('τ_cut')
    axs[1].legend(frameon=False, fontsize=8, loc='upper left')
    axs[1].annotate('uniform r at 1e-4, 3e-5, 1e-5: +1.3, +5.6, +18 (off scale)', (1.1e-5, -0.41), fontsize=7, color=MUTED)
    save(fig, out)


def dist(src, out, it='7'):
    it = int(it)
    d = json.load(open(src)); tc = d['tau_cut'][it]
    bins = sorted((b for b in d['bins'] if b != '0.05-0.50'), key=lambda b: float(b.split('-')[0]))
    lo = np.array([float(b.split('-')[0]) for b in bins]); hi = np.array([float(b.split('-')[1]) for b in bins])
    w = hi - lo; c = (lo + hi)/2
    dis = np.array([d['bins'][b]['disent_nlo'][0] for b in bins])/w
    dise = np.array([d['bins'][b]['disent_nlo'][1] for b in bins])/w
    our = np.array([d['bins'][b]['nlo'][0][it] for b in bins])/w
    oure = np.array([d['bins'][b]['nlo'][1][it] for b in bins])/w
    fig, axs = plt.subplots(2, 1, figsize=(6.2, 4.8), sharex=True, gridspec_kw={'height_ratios': [1.6, 1]})
    edges = np.append(lo, hi[-1])
    axs[0].stairs(dis, edges, color=ACC[0], lw=1.4, label='DISENT (NLO coefficient × 1/x)')
    axs[0].errorbar(c, our, oure, fmt='o', ms=3.5, lw=1, color=ACC[1], label='slicing: b1 + lo, τ_cut = %.0e' % tc)
    axs[0].set_ylabel('dσ/dx dQ² dτ_zQ, O(α_s²) [pb/GeV²]'); axs[0].legend(frameon=False, fontsize=8)
    axs[1].axhline(0, color=MUTED, lw=0.8)
    axs[1].errorbar(c, our/dis - 1, np.hypot(oure/dis, our*dise/dis**2), xerr=w/2, fmt='o', ms=3.5, lw=1, color=ACC[1])
    axs[1].set_ylabel('slicing/DISENT − 1'); axs[1].set_xlabel('τ_zQ'); axs[1].set_ylim(-0.08, 0.08)
    save(fig, out)


def zeusnlo(src, ref, out):
    d = json.load(open(src)); r = json.load(open(ref)); tc = np.array(d['tau_cut'])
    sel = [('total 0-0', 'total'), ('ptavg 8-15', 'p̄_T 8–15 (next to E_T > 8)'), ('m12 20-30', 'm₁₂ 20–30 (next to m₁₂ > 20)'),
           ('q2 125-250', 'Q² 125–250'), ('ptavg 15-22', 'p̄_T 15–22'), ('m12 45-65', 'm₁₂ 45–65')]
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for k, (b, lab) in enumerate(sel):
        v = d['bins'][b]; lo = r[b]['lo']; nj, nje = r[b]['nlo']
        y = (np.array(v['NLO'][0]) - nj)/lo*100; e = np.hypot(np.array(v['NLO'][1]), nje)/lo*100
        ax.errorbar(tc*(1 + 0.06*(k - 2.5)), y, e, fmt='o-', ms=3, lw=1, color=ACC[k % 6], label=lab)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xscale('log'); ax.set_xlabel('τ_cut'); ax.set_ylabel('(slicing − NNLOJET) NLO / LO  [%]')
    ax.set_ylim(-1, 12); ax.legend(frameon=False, fontsize=7, ncol=2)
    save(fig, out)


if __name__ == '__main__':
    {'nlo-disent': nlo_disent, 'nnlo': nnlo, 'compare': compare, 'dist': dist, 'zeusnlo': zeusnlo}[sys.argv[1]](*sys.argv[2:])
