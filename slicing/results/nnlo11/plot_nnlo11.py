#!/usr/bin/env python3
"""Plots of the NNLO DIS 1+1 three-way comparison (slicing/nnlo11.f90) from
nnlo11_results.csv (written by the export in the notebook entry of 2 Oct 2026).

All values are coefficients per Born in units of (alpha_s/2pi)^k, k = 1, 2, at
mu_R = mu_F = Q, photon exchange; bins are lab-frame jet observables (anti-k_t
R = 1, p_t > 5 GeV, -1 < eta < 2.5), the leading-jet bins refer to the hardest
jet inside that acceptance.

  python3 plot_nnlo11.py            # writes the PNG/PDF files next to this script
"""
import csv, os
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, 'nnlo11_results.csv'))))
for r in rows:
    for k in ('value', 'error', 'diff', 'differr'):
        r[k] = float(r[k])
    r['bin'] = int(r['bin']); r['order'] = int(r['order']); r['born_in_bin'] = int(r['born_in_bin'])

def get(point, cutoff, method, cut, order):
    out = {}
    for r in rows:
        if (r['point'], r['disent_cutoff'], r['method'], r['cut'], r['order']) == (point, cutoff, method, cut, order):
            out[r['bin']] = r
    return out

POINTS = [('x0.01_Q400', '1e-8', 'rho=0.001', r'$x=0.01$, $Q^2=400$ GeV$^2$ (240 runs, DISENT cutoff $10^{-8}$)'),
          ('x0.05_Q1000', '1e-10', 'rho=0.0003', r'$x=0.05$, $Q^2=1000$ GeV$^2$ (180 runs, DISENT cutoff $10^{-10}$)')]
STYLE = [('DISENT+P2B', '', 'k', 'o', 'DISENT + P2B'),
         ('slicing+P2B', 'tau2=1e-05', 'C0', 's', r'$\tau_2$ slicing + P2B, $\tau_2 > 10^{-5}$'),
         ('slicing+P2B', None, 'C2', 'D', None),
         ('pure tau1 slicing', 'tau1=0.0003', 'C1', '^', r'pure $\tau_1$ slicing, $\tau_{1,\rm cut} = 3\cdot10^{-4}$'),
         ('pure tau1 slicing', 'tau1=0.001', 'C3', 'v', r'pure $\tau_1$ slicing, $\tau_{1,\rm cut} = 10^{-3}$')]


def distributions(point, cutoff, rho, title):
    groups = [('leading-jet $p_t$ [GeV]', list(range(3, 10))), (r'leading-jet $\eta$ (lab)', list(range(10, 16)))]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharex='col', gridspec_kw=dict(height_ratios=[2, 1.3], hspace=0.05, wspace=0.18))
    for col, (xlabel, bins) in enumerate(groups):
        ref = get(point, cutoff, 'DISENT+P2B', '', 2)
        labs = []
        for b in bins:
            r = ref[b]
            labs.append(f'{float(r["lo"]):.3g}–{float(r["hi"]):.3g}' + (' ◆' if r['born_in_bin'] else ''))
        xs = np.arange(len(bins))
        for k, (meth, cut, c, mk, lab) in enumerate(STYLE):
            if cut is None:
                cut = rho; lab = rf'$\tau_2$ slicing + P2B, $\tau_2 > {float(rho.split("=")[1]):g}\,\tau_1$'
            d = get(point, cutoff, meth, cut, 2)
            if not d:
                continue
            off = (k - 2) * 0.12
            v = np.array([d[b]['value'] for b in bins]); e = np.array([d[b]['error'] for b in bins])
            dd = np.array([d[b]['diff'] for b in bins]); de = np.array([d[b]['differr'] for b in bins])
            axes[0, col].errorbar(xs + off, v, yerr=e, fmt=mk, color=c, ms=5, capsize=2, lw=1.2, label=lab)
            if meth != 'DISENT+P2B':
                axes[1, col].errorbar(xs + off, dd, yerr=de, fmt=mk, color=c, ms=5, capsize=2, lw=1.2)
        axes[0, col].axhline(0, color='0.6', lw=0.8)
        axes[1, col].axhline(0, color='k', lw=0.8)
        axes[1, col].set_xticks(xs); axes[1, col].set_xticklabels(labs, rotation=20, ha='right', fontsize=8.5)
        axes[1, col].set_xlabel(xlabel)
        axes[0, col].set_ylabel(r'$C_2$ per Born  [$(\alpha_s/2\pi)^2$]')
        axes[1, col].set_ylabel('method − DISENT+P2B')
        lim = max(10.0, 1.3 * np.max(np.abs([get(point, cutoff, 'slicing+P2B', 'tau2=1e-05', 2)[b]['diff'] for b in bins])) + 10)
        axes[1, col].set_ylim(-max(60, lim), max(60, lim))
        axes[0, col].grid(alpha=0.25); axes[1, col].grid(alpha=0.25)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', ncol=3, fontsize=8.5, bbox_to_anchor=(0.5, 0.955), frameon=False)
    fig.suptitle('NNLO DIS 1+1, O($\\alpha_s^2$) coefficient per bin: ' + title + '   (◆ = bin containing the Born jet)', fontsize=11, y=1.0)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(HERE, f'nnlo11_{point}_jets.{ext}'), dpi=130, bbox_inches='tight')
    plt.close(fig)


def pure_total():
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    sets = [('x0.01_Q400', '1e-8', 'C0', 'o', r'$x=0.01$, cutoff $10^{-8}$'), ('x0.01_Q400', '1e-10', 'C2', 's', r'$x=0.01$, cutoff $10^{-10}$'),
            ('x0.01_Q400', '1e-6', 'C3', 'v', r'$x=0.01$, cutoff $10^{-6}$'), ('x0.05_Q1000', '1e-10', 'C1', 'D', r'$x=0.05$, cutoff $10^{-10}$')]
    for k, (pt, cut, c, mk, lab) in enumerate(sets):
        d = defaultdict(dict)
        for r in rows:
            if r['point'] == pt and r['disent_cutoff'] == cut and r['method'] == 'pure tau1 slicing' and r['order'] == 2 and r['bin'] == 1:
                d[float(r['cut'].split('=')[1])] = (r['diff'], r['differr'])
        t = np.array(sorted(d)); v = np.array([d[x][0] for x in t]); e = np.array([d[x][1] for x in t])
        ax.errorbar(t * (1 + 0.06 * (k - 1.5)), v, yerr=e, fmt=mk + '-', color=c, ms=5, capsize=2, lw=1, label=lab)
    ax.set_xscale('log'); ax.set_yscale('symlog', linthresh=50)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xlabel(r'$\tau_{1,\rm cut}$'); ax.set_ylabel(r'(sum − exact) per Born  [$(\alpha_s/2\pi)^2$]')
    ax.set_title(r'Pure $\tau_1$ slicing at NNLO, inclusive: deviation from the exact $C_2$' + '\n' + r'(exact $C_2=-34.27$ at $x=0.01$, $-13.74$ at $x=0.05$; symlog axis, linear for $|y|<50$)', fontsize=10)
    ax.grid(alpha=0.25, which='both'); ax.legend(fontsize=8.5)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(HERE, f'nnlo11_pure_tau1_total.{ext}'), dpi=130, bbox_inches='tight')
    plt.close(fig)


def pull_matrix(point, cutoff, title, suffix=''):
    cuts = [('tau2=0.001', r'$10^{-3}$'), ('tau2=0.0001', r'$10^{-4}$'), ('tau2=1e-05', r'$10^{-5}$'), ('rho=0.01', r'$\rho\,10^{-2}$'),
            ('rho=0.003', r'$\rho\,3{\cdot}10^{-3}$'), ('rho=0.001', r'$\rho\,10^{-3}$'), ('rho=0.0003', r'$\rho\,3{\cdot}10^{-4}$'), ('rho=0.0001', r'$\rho\,10^{-4}$')]
    ref = get(point, cutoff, 'DISENT+P2B', '', 2)
    bins = list(range(2, 17))
    M = np.zeros((len(bins), len(cuts))); D = np.zeros_like(M)
    for j, (c, _) in enumerate(cuts):
        d = get(point, cutoff, 'slicing+P2B', c, 2)
        for i, b in enumerate(bins):
            M[i, j] = d[b]['diff'] / d[b]['differr'] if d[b]['differr'] > 0 else 0
            D[i, j] = d[b]['diff']
    fig, ax = plt.subplots(figsize=(9, 6.2))
    im = ax.imshow(np.clip(M, -5, 5), cmap='RdBu_r', vmin=-5, vmax=5, aspect='auto')
    for i in range(len(bins)):
        for j in range(len(cuts)):
            ax.text(j, i, f'{D[i, j]:+.1f}\n({M[i, j]:+.1f}σ)', ha='center', va='center', fontsize=6.5)
    ax.set_xticks(range(len(cuts))); ax.set_xticklabels([l for _, l in cuts], fontsize=8.5)
    ax.set_yticks(range(len(bins)))
    ax.set_yticklabels([ref[b]['label'] + ('  ◆' if ref[b]['born_in_bin'] else '') + f'   [{ref[b]["value"]:.1f}]' for b in bins], fontsize=8)
    ax.set_xlabel(r'$\tau_2$ cut: absolute $\tau_2 > c$, or relative $\tau_2 > \rho\,\tau_1$')
    ax.set_title(r'$\tau_2$ slicing + P2B minus DISENT + P2B (per Born; [ ] = DISENT + P2B $C_2$): ' + title, fontsize=9.5)
    fig.colorbar(im, ax=ax, label='pull (clipped at ±5)')
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(HERE, f'nnlo11_{point}_pulls{suffix}.{ext}'), dpi=130, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    for p, cut, rho, title in POINTS:
        distributions(p, cut, rho, title)
        pull_matrix(p, cut, title)
    # at x = 0.01 the relative cuts rho <= 3e-4 reach DISENT's technical cutoff 1e-8: the 1e-10 runs as well
    pull_matrix('x0.01_Q400', '1e-10', r'$x=0.01$, $Q^2=400$ GeV$^2$ (90 runs, DISENT cutoff $10^{-10}$)', suffix='_cut1e-10')
    pure_total()
    print('written to', HERE)
