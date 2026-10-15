#!/usr/bin/env python3
"""G8 (UPDATE 15 Oct): J8 window (x, Q2) = (0.1, 5000), vi/kp/r with nlo31's `flat` generator (runs ours8j/G8, seeds 6001-) against J8 (psmc, ours8j/PH9).
  combineg8.py -> results/G8_{parts,main,slopes}.txt, G8.json, results/g8/G8_<part>.zcell.gz (raw ZCELL rows of every done flat seed)
Sets: J8 = all parts psmc (PH9); G8 = vi, kp, r flat + b2, rcorr, d2 (and b1, lo for the NLO) of J8; G8h = kp, r flat + vi psmc (+ the same J8 parts).
Plain equal-weight means over done seeds, error = sd/sqrt(N), parts in quadrature; NNLO = b2+vi+kp+r+rcorr (+d2 for tau_cut < 1e-4), as combine8j.py.
Per part and cell (total, ptavg 8-15): flat and psmc means, pull (flat - psmc)/sqrt(e1^2+e2^2), per-seed sd, 1.48 MAD, median, fraction beyond 20 MAD.
Slopes as combine8j.py: (N(a)-N(b))/ln(a/b) in % of the J8 bin NLO at 1e-5, error sqrt|ea^2-eb^2|/ln, and also in 1e-9 pb/GeV2 per unit ln."""
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
    return dd[0][0], dd[0][1], np.array([x[2] for x in dd]), seeds, len(ds)
parts = ['b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2']; FL = ['vi', 'kp', 'r']
X = {}; INFO = {}
for p in parts:
    tc, names, a, s, n = load('PH9', p); X['J8', p] = a; INFO['J8 ' + p] = (n, len(a))
os.makedirs(R + '/g8', exist_ok=True)
for p in FL:
    tc2, names2, a, s, n = load('G8', p); assert names2 == names and tc2 == tc
    X['G8', p] = a; INFO['G8 ' + p] = (n, len(a))
    with gzip.open('%s/g8/G8_%s.zcell.gz' % (R, p), 'wt') as f:
        f.write('# G8 part %s (flat generator): raw ZCELL rows [pb/GeV2] per seed\n# tau_cut %s\n# columns: seed observable lo hi values(tau_cut columns)\n' % (p, ' '.join('%g' % t for t in tc)))
        for sd, row in zip(s, a):
            for nm, v in zip(names, row): f.write('%d %s %s\n' % (sd, nm.replace('-', ' '), ' '.join('%.8E' % q for q in v)))
tc = np.array(tc); lowd = (tc < 1e-4 - 1e-12).astype(float)
M = {k: (v.mean(0), v.std(0, ddof=1)/math.sqrt(len(v))) for k, v in X.items()}
SRC = {'J8': {p: 'J8' for p in parts}, 'G8': {p: ('G8' if p in FL else 'J8') for p in parts},
       'G8h': {p: ('G8' if p in ('kp', 'r') else 'J8') for p in parts}}
def comb(src, spec):
    m = 0; v = 0
    for p, w in spec: mm, ee = M[src[p], p]; m = m + mm*w; v = v + (ee*w)**2
    return m, np.sqrt(v)
NN = [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')] + [('d2', lowd)]
C = {s: comb(SRC[s], NN) for s in SRC}; NLO = comb(SRC['J8'], [('b1', 1), ('lo', 1)])
CELLS = ['total 0-0', 'ptavg 8-15']; ci = {c: names.index(c) for c in CELLS}
hdr = '%-10s' % 'tau_cut' + ''.join('%12.0e' % t for t in tc)
# per-part comparison
Lp = ['G8 per part: flat (G8) vs psmc (J8); units pb/GeV2; seeds ' + ', '.join('%s %d/%d' % (k, v[1], v[0]) for k, v in INFO.items()), '']
J = {'tau_cut': list(tc), 'seeds': {k: v[1] for k, v in INFO.items()}, 'parts': {}, 'nnlo': {}, 'slopes': {}}
for p in FL:
    for c in CELLS:
        i = ci[c]; Lp += ['== %s, %s' % (p, c), hdr]; J['parts']['%s %s' % (p, c)] = {}
        for k in ('J8', 'G8'):
            a = X[k, p][:, i, :]; med = np.median(a, 0); mad = 1.4826*np.median(abs(a - med), 0)
            fr = (abs(a - med) > 20*mad).mean(0)
            st = {'mean': a.mean(0), 'err': a.std(0, ddof=1)/math.sqrt(len(a)), 'sd': a.std(0, ddof=1), 'mad': mad, 'median': med, 'f20': fr}
            for q in ('mean', 'err', 'sd', 'mad', 'median', 'f20'):
                Lp.append('%-10s' % ('%s %s' % (k, q)) + ''.join(('%12.4g' if q != 'f20' else '%12.4f') % v for v in st[q]))
            J['parts']['%s %s' % (p, c)][k] = {q: list(v) for q, v in st.items()}
        pu = (M['G8', p][0][i] - M['J8', p][0][i])/np.sqrt(M['G8', p][1][i]**2 + M['J8', p][1][i]**2)
        Lp.append('%-10s' % 'pull G-J' + ''.join('%12.1f' % v for v in pu)); Lp.append(''); J['parts']['%s %s' % (p, c)]['pull'] = list(pu)
open(R + '/G8_parts.txt', 'w').write('\n'.join(Lp) + '\n')
Lm = ['G8 NNLO (pb/GeV2); J8 = psmc; G8 = flat vi, kp, r; G8h = flat kp, r with psmc vi; other parts from J8 in all three', '']
for c in CELLS:
    i = ci[c]; Lm += ['== ' + c, hdr, '%-10s' % 'NLO J8' + ''.join('%12.5g' % v for v in NLO[0][i])]
    for s in SRC:
        m, e = C[s]; Lm += ['%-10s' % s + ''.join('%12.5g' % v for v in m[i]), '%-10s' % '  +-' + ''.join('%12.2g' % v for v in e[i])]
        J['nnlo']['%s %s' % (s, c)] = [list(m[i]), list(e[i])]
    for s in ('G8', 'G8h'):
        pu = (C[s][0][i] - C['J8'][0][i])/np.sqrt(C[s][1][i]**2 + C['J8'][1][i]**2)
        Lm.append('%-10s' % ('pull %s-J8' % s) + ''.join('%12.1f' % v for v in pu))
    Lm.append('')
open(R + '/G8_main.txt', 'w').write('\n'.join(Lm) + '\n')
wins = [(5e-3, 5e-4), (2e-3, 2e-4), (1e-3, 1e-4), (1e-4, 1e-5)]
Ls = ['G8 NNLO slopes per unit ln tau_cut: % of the J8 bin |NLO(1e-5)| +- err (pull) | absolute in 1e-9 pb/GeV2 +- err', '%-20s ' % 'set cell' + ' '.join('%-36s' % ('%.0e->%.0e' % w) for w in wins)]
for c in CELLS:
    i = ci[c]; den = abs(NLO[0][i, -1])
    for s in SRC:
        m, e = C[s]; cs = []; js = []
        for a, b in wins:
            ja = int(np.argmin(abs(tc - a))); jb = int(np.argmin(abs(tc - b))); d = m[i, ja] - m[i, jb]; ed = math.sqrt(abs(e[i, ja]**2 - e[i, jb]**2)); Lg = math.log(a/b)
            v = [100*d/Lg/den, 100*ed/Lg/den, d/ed if ed > 0 else 0, 1e9*d/Lg, 1e9*ed/Lg]; js.append(v)
            cs.append('%+7.1f+-%-5.1f(%+5.1f) | %+6.1f+-%-5.1f' % tuple(v))
        Ls.append('%-20s ' % ('%s %s' % (s, c)) + ' '.join('%-36s' % x for x in cs)); J['slopes']['%s %s' % (s, c)] = js
open(R + '/G8_slopes.txt', 'w').write('\n'.join(Ls) + '\n')
json.dump(J, open(R + '/G8.json', 'w'), indent=1, default=float)
print('\n'.join(Ls))
