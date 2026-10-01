#!/usr/bin/env python3
"""Combine tau2_nlo.dat files of independent runs (same number of events
each, different seeds): integrals are averaged, errors from the
per-event variances. Prints, per tau_zQ bin, DISENT's NLO coefficient
and for every tau_cut (below + above - ref)/ref with its error."""
import sys, glob
import numpy as np

def load(f):
    w = open(f).read().split()
    nevt, nt, nb = int(w[0]), int(w[1]), int(w[2]); i = 3
    def take(n):
        nonlocal i
        a = np.array([float(v.replace('D', 'E')) for v in w[i:i + n]]); i += n; return a
    taus = take(nt); blo = take(nb); bhi = take(nb)
    ref = take(nb), take(nb)
    ab = take(nt * nb).reshape(nb, nt).T, take(nt * nb).reshape(nb, nt).T
    be = take(nt * nb).reshape(nb, nt).T, take(nt * nb).reshape(nb, nt).T
    d = take(nt * nb).reshape(nb, nt).T, take(nt * nb).reshape(nb, nt).T
    return nevt, taus, blo, bhi, ref, ab, be, d

files = sorted(sum((glob.glob(a) for a in sys.argv[1:]), []))
runs = [load(f) for f in files]
taus, blo, bhi = runs[0][1:4]
def comb(k):
    # each run: integral I_r = s1, var_r = s2 - s1^2/n ; combined = mean over runs
    I = np.array([r[k][0] for r in runs]); n = np.array([r[0] for r in runs], float)
    V = np.array([r[k][1] - r[k][0] ** 2 / r[0] for r in runs])
    return I.mean(0), np.sqrt(np.clip(V, 0, None).sum(0)) / len(runs)
ref, eref = comb(4); ab, eab = comb(5); be, ebe = comb(6); d, ed = comb(7)
print(f'{len(runs)} runs, {sum(r[0] for r in runs):,} events')
for ib in range(len(blo)):
    print(f'\ntau_zQ in [{blo[ib]:.2f},{bhi[ib]:.2f}): NLO coefficient (DISENT) {ref[ib]:.6e} +- {eref[ib]:.2e}')
    print(f'  {"tau_cut":>8s} {"below":>13s} {"above":>13s} {"sum-ref":>13s} {"(sum-ref)/ref":>14s} {"+-":>9s}')
    for it, t in enumerate(taus):
        print(f'  {t:8.1e} {be[it, ib]:13.5e} {ab[it, ib]:13.5e} {d[it, ib]:13.5e} {d[it, ib] / ref[ib]:14.4e} {ed[it, ib] / abs(ref[ib]):9.2e}')
