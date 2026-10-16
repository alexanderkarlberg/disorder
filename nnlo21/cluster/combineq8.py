#!/usr/bin/env python3
"""Q8 (UPDATE 16 Oct): CH8 quark channel (mask 1, ours8j/CH8q, all parts) + the Q8 top-up of r, rcorr, d2 (ours8j/Q8, seeds 21001-).
  combineq8.py -> results/Q8_{main,slopes,parts}.txt, Q8.json, results/q8/Q8_<part>.zcell.gz (raw ZCELL rows of every new seed)
Sets: CH8q = CH8 q seeds only; Q8new = Q8 seeds for r, rcorr, d2 + CH8 q for the other parts; Q8 = equal-weight means over old + new seeds;
Q8f = Q8 with kp from bin11 (40c013a, the kp_conv CF fix; ours8j/Q8/kp11, the CH8 q kp seeds 7001-7200 rerun). Q8f is the result; the others use bin10's kp.
NNLO = b2+vi+kp+r+rcorr (+d2 for tau_cut < 1e-4), NLO = b1+lo, as combine8j.py. Slopes (N(a)-N(b))/ln(a/b) with two errors:
  e1 = sqrt|ea^2-eb^2|/ln (the combine8j definition, used in all earlier reports) and
  e2 = sqrt(sum_parts var_seed(N_p(a)-N_p(b))/N_p)/ln (per-seed difference, the definition of the Q8 sizing; d2 enters as its value below 1e-4).
In 1e-9 pb/GeV2 per unit ln and in % of the CH8 q bin |NLO(1e-5)|."""
import glob, gzip, json, math, os, re
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs/ours8j'; H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results'
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
def load(S, p):
    ds = sorted(glob.glob('%s/%s/%s/s*/' % (P, S, p)), key=lambda s: int(re.search(r's(\d+)/$', s).group(1)))
    dd, seeds = [], []
    for d in ds:
        if not os.path.exists(d + 'done'): continue
        x = rows(d + 'run.log')
        if x[2] is not None: dd.append(x); seeds.append(int(re.search(r's(\d+)/$', d).group(1)))
    a = np.array([x[2] for x in dd])
    return dd[0][0], dd[0][1], a, seeds, (len(ds), len(dd), int((~np.isfinite(a).all(axis=(1, 2))).sum()))
parts = ['b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2']; NEW = ['r', 'rcorr', 'd2']
X = {}; INFO = {}; names = tc0 = None
os.makedirs(R + '/q8', exist_ok=True)
for p in parts:
    tc, nm, a, s, info = load('CH8q', p); X['old', p] = a; INFO['CH8q ' + p] = info
    if names is None: names, tc0 = nm, tc
    assert nm == names and tc == tc0
for p in NEW:
    tc, nm, a, s, info = load('Q8', p); assert nm == names and tc == tc0; X['new', p] = a; INFO['Q8 ' + p] = info
    with gzip.open('%s/q8/Q8_%s.zcell.gz' % (R, p), 'wt') as f:
        f.write('# Q8 part %s (quark channel, PDFMASK 1): raw ZCELL rows [pb/GeV2] per seed\n# tau_cut %s\n# columns: seed observable lo hi values(tau_cut columns)\n' % (p, ' '.join('%g' % t for t in tc)))
        for sd, row in zip(s, a):
            for n, v in zip(names, row): f.write('%d %s %s\n' % (sd, n.replace('-', ' '), ' '.join('%.8E' % q for q in v)))
tc, nm, a, s, info = load('Q8', 'kp11'); assert nm == names and tc == tc0; X['fix', 'kp'] = a; INFO['Q8 kp11 (bin11)'] = info
SK = s; SKO = load('CH8q', 'kp')[3]; assert SK == SKO, 'kp11 seeds differ from CH8q kp'
tc = np.array(tc0); lowd = (tc < 1e-4 - 1e-12).astype(float)
SETS = {'CH8q': {p: X['old', p] for p in parts}, 'Q8new': {p: (X['new', p] if p in NEW else X['old', p]) for p in parts},
        'Q8': {p: (np.concatenate([X['old', p], X['new', p]]) if p in NEW else X['old', p]) for p in parts}}
SETS['Q8f'] = dict(SETS['Q8'], kp=X['fix', 'kp'])
NN = [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')] + [('d2', lowd)]
def comb(D, spec):
    m = 0; v = 0
    for p, w in spec: a = D[p]; m = m + a.mean(0)*w; v = v + (a.std(0, ddof=1)/math.sqrt(len(a))*w)**2
    return m, np.sqrt(v)
C = {s: comb(D, NN) for s, D in SETS.items()}; NLO = comb(SETS['CH8q'], [('b1', 1), ('lo', 1)])
CELLS = ['total 0-0', 'ptavg 8-15', 'q2 2000-5000', 'q2 5000-20000']; ci = {c: names.index(c) for c in CELLS}
TS = [2e-4, 1e-4, 3e-5, 1e-5]; jt = {t: int(np.argmin(abs(tc - t))) for t in TS}
J = {'tau_cut': list(tc), 'seeds': INFO, 'nnlo': {}, 'slopes': {}, 'decomp': {}}
L = ['Q8: quark channel (PDFMASK 1), J8 window; NNLO [pb/GeV2]; CH8q = old seeds only, Q8new = new r/rcorr/d2 seeds (other parts CH8q), Q8 = old + new',
     'bookkeeping (dirs, done+parsed, with non-finite): ' + '; '.join('%s %s' % (k, v) for k, v in INFO.items()), '']
hdr = '%-10s' % 'tau_cut' + ''.join('%12.0e' % t for t in tc)
for c in CELLS:
    i = ci[c]; L += ['== ' + c, hdr, '%-10s' % 'NLO' + ''.join('%12.5g' % v for v in NLO[0][i])]
    for s in SETS:
        m, e = C[s]; L += ['%-10s' % s + ''.join('%12.5g' % v for v in m[i]), '%-10s' % '  +-' + ''.join('%12.2g' % v for v in e[i])]
        J['nnlo']['%s %s' % (s, c)] = [list(m[i]), list(e[i])]
    m0, e0 = C['Q8'][0][i], C['Q8'][1][i]; dk = (X['fix', 'kp'] - X['old', 'kp'])[:, i, :]
    L += ['%-10s' % 'dkp' + ''.join('%12.4g' % v for v in dk.mean(0)), '%-10s' % '  +-' + ''.join('%12.2g' % v for v in dk.std(0, ddof=1)/math.sqrt(len(dk)))]
    pu = (C['Q8new'][0][i] - C['CH8q'][0][i])/np.sqrt(C['Q8new'][1][i]**2 + C['CH8q'][1][i]**2)
    L += ['%-10s' % 'pull n-o' + ''.join('%12.1f' % v for v in pu), '']
open(R + '/Q8_main.txt', 'w').write('\n'.join(L) + '\n')
WINS = [(5e-3, 5e-4), (1e-3, 1e-4), (2e-4, 1e-4), (1e-4, 3e-5), (1e-4, 1e-5), (3e-5, 1e-5)]
def dpart(a, p, ja, jb, i):
    """per-seed contribution of part p to N(a) - N(b) in cell i"""
    wa = lowd[ja] if p == 'd2' else 1; wb = lowd[jb] if p == 'd2' else 1
    return a[:, i, ja]*wa - a[:, i, jb]*wb
Ls = ['Q8 NNLO slopes per unit ln tau_cut, 1e-9 pb/GeV2: value +- e1 [+- e2] (value/e2) {%% of the CH8 q bin |NLO(1e-5)|}; e1 = sqrt|ea^2-eb^2|/ln, e2 = per-seed difference', '']
for c in CELLS:
    i = ci[c]; den = abs(NLO[0][i, -1]); Ls.append('== %s (|NLO(1e-5)| = %.4g e-9)' % (c, den*1e9))
    for s, D in SETS.items():
        m, e = C[s]; cs = []; js = []
        for a, b in WINS:
            ja = int(np.argmin(abs(tc - a))); jb = int(np.argmin(abs(tc - b))); Lg = math.log(a/b); d = m[i, ja] - m[i, jb]
            e1 = math.sqrt(abs(e[i, ja]**2 - e[i, jb]**2)); e2 = math.sqrt(sum(dpart(D[p], p, ja, jb, i).var(ddof=1)/len(D[p]) for p, _ in NN))
            v = [1e9*d/Lg, 1e9*e1/Lg, 1e9*e2/Lg, d/e2, 100*d/Lg/den]; js.append(v)
            cs.append('%+7.2f+-%5.2f[%5.2f](%+5.1f){%+6.1f%%}' % tuple(v))
        Ls.append('%-6s ' % s + ' '.join('%-40s' % ('%.0e->%.0e: %s' % (a, b, x)) for (a, b), x in zip(WINS, cs))); J['slopes']['%s %s' % (s, c)] = js
    Ls.append('')
open(R + '/Q8_slopes.txt', 'w').write('\n'.join(Ls) + '\n')
# per-part decomposition of N(a) - N(b) for the low windows
Lp = ['Q8 per-part decomposition of N(a) - N(b) (1e-9 pb/GeV2; not divided by ln): mean +- per-seed-difference error (N seeds); set Q8 (old + new); also CH8q and Q8new for r, rcorr, d2; kp as bin10 (Q8) and bin11 (Q8f); the sum is Q8f', '']
for c in CELLS:
    i = ci[c]; Lp.append('== ' + c)
    for a, b in [(1e-4, 1e-5), (1e-4, 3e-5), (3e-5, 1e-5), (2e-4, 1e-4)]:
        ja = jt[a]; jb = jt[b]; row = []; tot = 0; tv = 0
        for p, _ in NN:
            for s in (('CH8q', 'Q8new', 'Q8') if p in NEW else ('Q8', 'Q8f') if p == 'kp' else ('Q8',)):
                x = dpart(SETS[s][p], p, ja, jb, i); mm = x.mean()*1e9; ee = x.std(ddof=1)/math.sqrt(len(x))*1e9
                if s == 'Q8f' or (s == 'Q8' and p != 'kp'): tot += mm; tv += ee**2
                row.append('%s%s %+.2f+-%.2f (%d)' % (p, '' if s == 'Q8' else '[%s]' % s, mm, ee, len(x))); J['decomp']['%s %s %g-%g %s' % (c, p, a, b, s)] = [mm, ee, len(x)]
        x = dpart(X['fix', 'kp'] - X['old', 'kp'], 'kp', ja, jb, i)*1e9; row.append('dkp(bin11-bin10, same seeds) %+.3f+-%.3f' % (x.mean(), x.std(ddof=1)/math.sqrt(len(x))))
        Lp.append('  %.0e -> %.0e: sum (Q8f) %+.2f +- %.2f | ' % (a, b, tot, math.sqrt(tv)) + '; '.join(row))
    Lp.append('')
open(R + '/Q8_parts.txt', 'w').write('\n'.join(Lp) + '\n')
json.dump(J, open(R + '/Q8.json', 'w'), indent=1, default=float)
print('\n'.join(L[:2])); print('\n'.join(Ls)); print('\n'.join(Lp))
