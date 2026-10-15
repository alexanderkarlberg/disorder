#!/usr/bin/env python3
"""CH8 (UPDATE 15 Oct (b)): the J8 window (x, Q2) = (0.1, 5000) split by the incoming parton, bin10 PDFMASK 1 (quarks + antiquarks,
runs ours8j/CH8q, seeds 7001-) and 2 (gluon, CH8g, seeds 8001-), every part with the same mask; J8 = ours8j/PH9 (unmasked, psmc).
  combinech8.py -> results/CH8_{main,slopes,pulls}.txt, CH8.json, results/ch8/CH8{q,g}_<part>.zcell.gz (raw ZCELL rows of every done seed)
Plain equal-weight means over done seeds, error = sd/sqrt(N), parts in quadrature; NNLO = b2+vi+kp+r+rcorr (+d2 for tau_cut < 1e-4),
NLO = b1+lo, as combine8j.py. q+g: sum of the channel means, errors in quadrature (independent seeds).
Slopes as combine8j.py: (N(a)-N(b))/ln(a/b) in % of the set's own |NLO(1e-5)| in the bin, error sqrt|ea^2-eb^2|/ln, and absolute in
1e-9 pb/GeV2 per unit ln. Pulls ((q+g) - J8)/sqrt(e1^2+e2^2) per cell and tau_cut for NLO and NNLO."""
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
parts = ['b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2']; DIR = {'q': 'CH8q', 'g': 'CH8g', 'J8': 'PH9'}
M = {}; INFO = {}; tc0 = names = None
os.makedirs(R + '/ch8', exist_ok=True)
for s, S in DIR.items():
    for p in parts:
        tc, nm, a, seeds, info = load(S, p); INFO['%s %s' % (s, p)] = info
        if names is None: tc0, names = tc, nm
        assert nm == names and tc == tc0
        M[s, p] = (a.mean(0), a.std(0, ddof=1)/math.sqrt(len(a)))
        if s == 'J8': continue
        with gzip.open('%s/ch8/%s_%s.zcell.gz' % (R, S, p), 'wt') as f:
            f.write('# CH8 %s (PDFMASK %d) part %s: raw ZCELL rows [pb/GeV2] per seed\n# tau_cut %s\n# columns: seed observable lo hi values(tau_cut columns)\n'
                    % (S, 1 if s == 'q' else 2, p, ' '.join('%g' % t for t in tc)))
            for sd, row in zip(seeds, a):
                for n, v in zip(names, row): f.write('%d %s %s\n' % (sd, n.replace('-', ' '), ' '.join('%.8E' % q for q in v)))
tc = np.array(tc0); lowd = (tc < 1e-4 - 1e-12).astype(float)
SPEC = {'NLO': [('b1', 1), ('lo', 1)], 'NNLO': [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')] + [('d2', lowd)]}
C = {}
for s in DIR:
    for lab, sp in SPEC.items():
        m = 0; v = 0
        for p, w in sp: m = m + M[s, p][0]*w; v = v + (M[s, p][1]*w)**2
        C[s, lab] = (m, np.sqrt(v))
for lab in SPEC: C['q+g', lab] = (C['q', lab][0] + C['g', lab][0], np.sqrt(C['q', lab][1]**2 + C['g', lab][1]**2))
SETS = ['q', 'g', 'q+g', 'J8']
active = [i for i in range(len(names)) if abs(C['J8', 'NLO'][0][i, -1]) > 0]
hdr = '%-14s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc)
J = {'tau_cut': list(tc), 'seeds': INFO, 'bins': {}, 'slopes': {}, 'pulls': {}, 'nlo_den': {}}
L = ['CH8, window (x,Q2) = (0.1, 5000), ZEUS jets (+ ZCUTS rows); dsigma/dx dQ2 [pb/GeV2]; q = PDFMASK 1 (quarks + antiquarks), g = PDFMASK 2, J8 = unmasked',
     'bookkeeping (dirs, done+parsed, with non-finite): ' + '; '.join('%s %s' % (k, v) for k, v in INFO.items())]
for i in active:
    L += ['', '== ' + names[i], hdr]; J['bins'][names[i]] = {}
    for lab in SPEC:
        for s in SETS:
            m, e = C[s, lab]; L += ['%-14s' % ('%s %s' % (lab, s)) + ''.join('%11.5g' % v for v in m[i]), '%-14s' % '  +-' + ''.join('%11.2g' % v for v in e[i])]
            J['bins'][names[i]]['%s %s' % (lab, s)] = [list(m[i]), list(e[i])]
open(R + '/CH8_main.txt', 'w').write('\n'.join(L) + '\n')
Lp = ['CH8 pulls ((q+g) - J8)/sqrt(e1^2+e2^2) per cell and tau_cut', '%-22s' % 'cell' + ''.join('%9.0e' % t for t in tc)]
allp = []
for lab in SPEC:
    for i in active:
        pu = (C['q+g', lab][0][i] - C['J8', lab][0][i])/np.sqrt(C['q+g', lab][1][i]**2 + C['J8', lab][1][i]**2)
        Lp.append('%-22s' % ('%s %s' % (lab, names[i])) + ''.join('%9.1f' % v for v in pu)); J['pulls']['%s %s' % (lab, names[i])] = list(pu); allp += list(pu)
allp = np.array(allp)
Lp += ['', 'all %d pulls: mean %+.2f, rms %.2f, max |pull| %.1f, fraction |pull| > 2: %.3f, > 3: %.3f' % (len(allp), allp.mean(), math.sqrt((allp**2).mean()), abs(allp).max(), (abs(allp) > 2).mean(), (abs(allp) > 3).mean())]
open(R + '/CH8_pulls.txt', 'w').write('\n'.join(Lp) + '\n')
wins = [(5e-3, 5e-4), (2e-3, 2e-4), (1e-3, 1e-4), (1e-4, 1e-5)]
Ls = ['CH8 NNLO slopes per unit ln tau_cut, (N(a)-N(b))/ln(a/b): %% of the set\'s own bin |NLO(1e-5)| +- err (pull) | absolute in 1e-9 pb/GeV2 +- err',
      '%-26s %10s ' % ('set cell', '|NLO| 1e-9') + ' '.join('%-36s' % ('%.0e->%.0e' % w) for w in wins)]
for i in active:
    for s in SETS:
        m, e = C[s, 'NNLO']; den = abs(C[s, 'NLO'][0][i, -1]); J['nlo_den']['%s %s' % (s, names[i])] = den; cs = []; js = []
        for a, b in wins:
            ja = int(np.argmin(abs(tc - a))); jb = int(np.argmin(abs(tc - b))); d = m[i, ja] - m[i, jb]; ed = math.sqrt(abs(e[i, ja]**2 - e[i, jb]**2)); Lg = math.log(a/b)
            v = [100*d/Lg/den if den > 0 else 0, 100*ed/Lg/den if den > 0 else 0, d/ed if ed > 0 else 0, 1e9*d/Lg, 1e9*ed/Lg]; js.append(v)
            cs.append('%+7.1f+-%-5.1f(%+5.1f) | %+6.2f+-%-5.2f' % tuple(v))
        Ls.append('%-26s %10.4g ' % ('%s %s' % (s, names[i]), den*1e9) + ' '.join('%-36s' % x for x in cs)); J['slopes']['%s %s' % (s, names[i])] = js
    Ls.append('')
open(R + '/CH8_slopes.txt', 'w').write('\n'.join(Ls) + '\n')
json.dump(J, open(R + '/CH8.json', 'w'), indent=1, default=float)
print('\n'.join(L[:2])); print(Lp[-1]); print('\n'.join(Ls))
