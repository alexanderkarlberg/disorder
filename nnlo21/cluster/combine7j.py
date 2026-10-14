#!/usr/bin/env python3
"""J7 (14 Oct): ZEUS jets at the fixed window (x, Q2) = (0.01, 400) (ZFIX=2, bin6), photon (PH) and W-. ZCELL rows (pb).
  combine7j.py PH|W  -> results/J7_SET_{main,parts,slopes}.txt, J7_SET.json, results/j7/SET_<part>.zcell.gz (raw ZCELL rows of every done seed)
Plain equal-weight means over done seeds, error = seed scatter/sqrt(N), parts in quadrature. NNLO = b2+vi+kp+r+rcorr (+d2 for tau_cut < 1e-4),
NLO = b1+lo (the 100+100 seeds of mk_j7nlo.py). Slopes: (N(a)-N(b))/ln(a/b) in % of the bin NLO at tau_cut = 1e-5 (own NLO, no external reference),
pull = (N(a)-N(b))/sqrt|ea^2-eb^2|. Bins with NLO = 0 (outside the window's Q2 bin) are skipped."""
import glob, gzip, json, math, os, re, sys
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs/ours7j'; H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results'
S = sys.argv[1]
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def rows(f):
    tc, nm, r = None, [], []
    for l in open(f):
        w = l.split()
        if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
        elif l.startswith(' ZCELL '): nm.append('%s %g-%g' % (w[1], float(w[2]), float(w[3]))); r.append([ff(x) for x in w[4:]])
    return tc, nm, (np.array(r) if r else None)
parts = ['b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2']; X = {}; INFO = {}
os.makedirs(R + '/j7', exist_ok=True)
for p in parts:
    ds = sorted(glob.glob('%s/%s/%s/s*/' % (P, S, p)), key=lambda s: int(re.search(r's(\d+)/$', s).group(1)))
    dd, seeds = [], []
    for d in ds:
        if not os.path.exists(d + 'done'): continue
        x = rows(d + 'run.log')
        if x[2] is not None: dd.append(x); seeds.append(int(re.search(r's(\d+)/$', d).group(1)))
    arr = np.array([x[2] for x in dd]); tc, names = dd[0][0], dd[0][1]
    X[p] = arr; INFO[p] = (len(ds), len(dd), int((~np.isfinite(arr).all(axis=(1, 2))).sum()))
    with gzip.open('%s/j7/%s_%s.zcell.gz' % (R, S, p), 'wt') as f:
        f.write('# J7 %s part %s: raw ZCELL rows [pb] per seed (sigma per bin below the cut; b2/vi/kp/r/rcorr/d2 as in mk_j7.py, b1/lo as mk_j7nlo.py)\n# tau_cut %s\n# columns: seed observable lo hi values(tau_cut columns)\n' % (S, p, ' '.join('%g' % t for t in tc)))
        for s, a in zip(seeds, arr):
            for n, v in zip(names, a): f.write('%d %s %s\n' % (s, n.replace(' ', ' ').replace('-', ' '), ' '.join('%.8E' % q for q in v)))
tc = np.array(tc); ntc = len(tc); lowd = (tc < 1e-4 - 1e-12).astype(float); cnt = {p: len(X[p]) for p in parts}
M = {p: (X[p].mean(0), X[p].std(0, ddof=1)/math.sqrt(len(X[p]))) for p in parts}
def comb(spec):
    m = 0; v = 0
    for p, w in spec: m = m + M[p][0]*w; v = v + (M[p][1]*w)**2
    return m, np.sqrt(v)
SPEC = {'NLO': [('b1', 1), ('lo', 1)], 'NNLO': [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')] + [('d2', lowd)],
        'NNLO(no d2)': [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')]}
C = {k: comb(v) for k, v in SPEC.items()}
nlo = np.abs(C["NLO"][0][:, -1])
active = [i for i, n in enumerate(names) if abs(nlo[i]) > 0]
hdr = '%-18s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc)
L = ['J7 %s, window (x,Q2) = (0.01, 400), ZEUS jets; dsigma/dx dQ2 [pb/GeV2] (the unit of the RESULT line; the ZCELLS header text says pb); seeds per part: ' % S + ' '.join('%s:%d' % (p, cnt[p]) for p in parts),
     'bookkeeping (dirs, done+parsed, with non-finite): ' + '; '.join('%s %s' % (k, v) for k, v in INFO.items())]
J = {'tau_cut': list(tc), 'seeds': cnt, 'bins': {}}
for i in active:
    L += ['', '== ' + names[i], hdr]; J['bins'][names[i]] = {}
    for lab in SPEC:
        m, e = C[lab]; L += ['%-18s' % lab + ''.join('%11.5g' % v for v in m[i]), '%-18s' % '  +-' + ''.join('%11.2g' % v for v in e[i])]
        J['bins'][names[i]][lab] = [list(m[i]), list(e[i])]
open('%s/J7_%s_main.txt' % (R, S), 'w').write('\n'.join(L) + '\n')
Lp = []
for i in active:
    Lp += ['== ' + names[i], hdr]
    for p in parts: Lp += ['%-14s' % p + ''.join('%11.4g' % v for v in M[p][0][i]), '%-14s' % '  +-' + ''.join('%11.2g' % v for v in M[p][1][i])]
    Lp.append('')
open('%s/J7_%s_parts.txt' % (R, S), 'w').write('\n'.join(Lp) + '\n')
wins = [(5e-3, 5e-4), (2e-3, 2e-4), (1e-3, 1e-4), (1e-4, 1e-5)]
Ls = ['%s: NNLO change per unit ln tau_cut, (N(a)-N(b))/ln(a/b), a > b, in %% of the bin NLO at tau_cut 1e-5 (own b1+lo); +- sqrt|ea^2-eb^2|/ln; (pull)' % S,
      '%-30s ' % 'bin' + ' '.join('%-24s' % ('%.0e->%.0e' % w) for w in wins)]
sl = {}
for lab in ('NNLO', 'NNLO(no d2)'):
    m, e = C[lab]
    for i in active:
        cs = []; js = []
        for a, b in wins:
            ja = int(np.argmin(abs(tc - a))); jb = int(np.argmin(abs(tc - b))); d = m[i, ja] - m[i, jb]; ed = math.sqrt(abs(e[i, ja]**2 - e[i, jb]**2)); Lg = math.log(a/b)
            v = [100*d/Lg/nlo[i], 100*ed/Lg/abs(nlo[i]), d/ed if ed > 0 else 0]; js.append(v); cs.append('%+8.2f+-%-6.2f(%+.1f)' % tuple(v))
        Ls.append('%-30s ' % ('%s %s' % (lab, names[i])) + ' '.join('%-24s' % c for c in cs)); sl['%s %s' % (lab, names[i])] = js
J['slopes'] = sl; J['nlo_den'] = {names[i]: nlo[i] for i in active}
open('%s/J7_%s_slopes.txt' % (R, S), 'w').write('\n'.join(Ls) + '\n')
json.dump(J, open('%s/J7_%s.json' % (R, S), 'w'), indent=1)
print('\n'.join(L[:2])); print('\n'.join(Ls))
