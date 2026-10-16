#!/usr/bin/env python3
"""bin6 round (13 Oct): G6W / G6P (fixed point x=0.1, Q2=5000, mode 1, tau_zQ CELL rows) and Ex6 (E P2B extras, VTARGET=0).
  combine6.py G6W|G6P|Ex6  -> results/SET_6_{main,parts,slopes}.txt, SET_6.json
Plain equal-weight means over done seeds (error = seed scatter/sqrt(N)), parts added in quadrature. Units for G6: dsigma/dx dQ2 [pb/GeV2];
DISENT reference = NLO coefficient / x (as combine_tcut). NLO = b1+lo; NNLO = b2+vi+kp+r+rcorr (+d2 for tau_cut < 1e-4, as combine5).
Ex6: bin5 plain parts (ours5/E) + Ex6 extras (lo-x, vi-x, kp-x, r-x); NNLO_P2B = NNLO + vi-x + kp-x + r-x, NLO_P2B = NLO + lo-x;
NNLO_P2B(bin5 extras) is the same with the old bin5 extras (ours5/E/*-x) for comparison. NNLOJET E: results/nnj5/E (as combine5)."""
import glob, json, math, os, re, sys
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'; H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results'
K11 = os.environ.get('KP11') == '1'  # KREDO: read kp11 / kp-x11 (bin11) instead of kp / kp-x, write to results/kp11
R0 = R  # inputs (NNLOJET references, earlier results) stay in results/
if K11: R = R + '/kp11'; os.makedirs(R, exist_ok=True)
def kd(p): return p + '11' if K11 and p in ('kp', 'kp-x') else p
S = sys.argv[1]; XF = 0.1
WREF = {(0.05, 0.1): (7.606e-4, 3.1e-6), (0.1, 0.2): (5.873e-4, 1.4e-6), (0.2, 0.3): (2.374e-4, 6.4e-7), (0.3, 0.4): (1.212e-4, 3.9e-7),
        (0.4, 0.5): (7.177e-5, 2.7e-7), (0.05, 0.5): (1.7783e-3, 2.8e-6)}   # CLUSTER-INSTRUCTIONS lines 33-36 (MPP runs/ccw1): already in combine_tcut units (= 'DISENT NLO x 1/x'), used as given (they agree with ours to <1%; dividing by x would be off by 10)
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def rows(f):
    tc, nm, r = None, [], []
    for l in open(f):
        w = l.split()
        if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
        elif l.startswith(' CELL '): nm.append(('tauzQ', float(w[1]), float(w[2]))); r.append([ff(x) for x in w[3:]])
    return tc, nm, (np.array(r) if r else None)
INFO = {}
def load(d, part):
    ds = sorted(glob.glob('%s/%s/s*/' % (d, kd(part)))); fs = [x + 'run.log' for x in ds if os.path.exists(x + 'done') and os.path.exists(x + 'run.log')]
    dd = [rows(f) for f in fs]; ok = [x for x in dd if x[2] is not None]
    arr = np.array([x[2] for x in ok]); nf = int((~np.isfinite(arr).all(axis=(1, 2))).sum())
    INFO[d.split('runs/')[1] + '/' + part] = (len(ds), len(fs), len(ok), nf)
    return ok[0][0], ok[0][1], arr
def stat(x): n = len(x); return x.mean(0), x.std(0, ddof=1)/math.sqrt(n)
if S == 'Ex6':
    src = {p: P + '/ours5/E' for p in ('b0', 'b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2')}
    src.update({p: P + '/ours6/Ex6' for p in ('lo-x', 'vi-x', 'kp-x', 'r-x')})
    parts = list(src)
    for p in ('lo-x', 'vi-x', 'kp-x', 'r-x'): src[p + '5'] = P + '/ours5/E'
    X = {}
    for p in src:
        tc, names, X[p] = load(src[p], p[:-1] + '' if p.endswith('5') and p[:-1] in ('lo-x', 'vi-x', 'kp-x', 'r-x') else p)
    parts = list(src)
else:
    parts = ['b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2']; X = {}
    for p in parts: tc, names, X[p] = load('%s/ours6/%s' % (P, S), p)
tc = np.array(tc); ntc = len(tc); nb = len(names); lowd = (tc < 1e-4 - 1e-12).astype(float)
M = {p: stat(X[p]) for p in parts}; cnt = {p: len(X[p]) for p in parts}
def comb(spec):
    m = 0; v = 0
    for p, w in spec: m = m + M[p][0]*w; v = v + (M[p][1]*w)**2
    return m, np.sqrt(v)
SPEC = {'NLO': [('b1', 1), ('lo', 1)], 'NNLO': [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')] + [('d2', lowd)],
        'NNLO(no d2)': [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')]}
if S == 'Ex6':
    SPEC = {'LO': [('b0', 1)], **SPEC}
    SPEC['NLO_P2B'] = SPEC['NLO'] + [('lo-x', 1)]
    SPEC['NNLO_P2B'] = SPEC['NNLO'] + [('vi-x', 1), ('kp-x', 1), ('r-x', 1)]
    SPEC['NNLO_P2Bextra'] = [('vi-x', 1), ('kp-x', 1), ('r-x', 1)]
    SPEC['NNLO_P2B(bin5 extras)'] = SPEC['NNLO'] + [('vi-x5', 1), ('kp-x5', 1), ('r-x5', 1)]
    SPEC['NNLO_P2Bextra(bin5 extras)'] = [('vi-x5', 1), ('kp-x5', 1), ('r-x5', 1)]
    SPEC['NLO_P2B(bin5 extras)'] = SPEC['NLO'] + [('lo-x5', 1)]
# references
ref = {}
if S == 'G6W':
    for k, (v, e) in WREF.items(): ref[('tauzQ',) + k] = (v, e)
elif S == 'G6P':
    fs = sorted(glob.glob(P + '/ours6/DISENT-ph/s*/tau2_nlo.dat')); fs = [f for f in fs if os.path.exists(os.path.dirname(f) + '/done')]
    def lt(f):
        w = open(f).read().split(); nt, nb_ = int(w[1]), int(w[2]); i = 3 + nt + 2*nb_
        return np.array([float(v.replace('D', 'E')) for v in w[i:i + nb_]]), [float(v) for v in w[3 + nt:3 + nt + nb_]], [float(v) for v in w[3 + nt + nb_:3 + nt + 2*nb_]]
    L = [lt(f) for f in fs]; a = np.array([x[0] for x in L]); lo_, hi_ = L[0][1], L[0][2]
    m, e = a.mean(0), a.std(0, ddof=1)/math.sqrt(len(a))
    INFO['DISENT-ph'] = (len(glob.glob(P + '/ours6/DISENT-ph/s*/')), len(fs), len(a), int((~np.isfinite(a).all(1)).sum()))
    for i in range(len(lo_)): ref[('tauzQ', round(lo_[i], 2), round(hi_[i], 2))] = (m[i]/XF, e[i]/XF)
elif S == 'Ex6':
    N = {}
    for o in ('LO', 'NLO', 'NNLO'):
        fn = '%s/nnj5/E/%s.cross.dat' % (R0, o); w = [float(v) for v in open(fn).read().split('\n')[1].split()]
        N[(o, 0.05, 0.5)] = (w[0]/1000, w[1]/1000)
        for l in open('%s/nnj5/E/%s.dis_tauzq.dat' % (R0, o)):
            if l.startswith('#') or not l.strip(): continue
            w = [float(v) for v in l.split()]; N[(o, w[0], w[2])] = (w[3]*(w[2] - w[0])/1000, w[4]*(w[2] - w[0])/1000)
    ref = N
lines = []; J = {'tau_cut': list(tc), 'seeds': cnt, 'bins': {}}
def pr(s=''): lines.append(s)
hdr = '%-22s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc)
pr('SET %s: seeds per part: ' % S + ' '.join('%s:%d' % (p, cnt[p]) for p in parts))
pr('seed bookkeeping (dirs, done, parsed, with non-finite cells): ' + '; '.join('%s %s' % (k, v) for k, v in INFO.items()))
def refkey(o, nm):
    if S == 'Ex6': return (o,) + nm[1:]
    return nm if o == 'NLO' else None
for i, nm in enumerate(names):
    k = '%s %g-%g' % nm; J['bins'][k] = {}; pr(); pr('== ' + k); pr(hdr)
    for lab, sp in SPEC.items():
        m, e = comb(sp); pr('%-22s' % lab + ''.join('%11.5g' % v for v in m[i])); pr('%-22s' % '  +-' + ''.join('%11.2g' % v for v in e[i]))
        J['bins'][k][lab] = [list(m[i]), list(e[i])]
        o = lab.split('_')[0].split('(')[0]
        if S == 'Ex6': r = ref.get((o,) + nm[1:]) if lab in ('LO', 'NLO', 'NNLO', 'NLO_P2B', 'NNLO_P2B') else None
        else: r = ref.get(nm) if lab == 'NLO' else None
        if r:
            pr('%-22s %.5g +- %.2g  pulls: ' % ('  ref ' + o, r[0], r[1]) + ' '.join('%+.1f' % ((m[i][c] - r[0])/math.hypot(e[i][c], r[1])) for c in range(ntc)))
            J['bins'][k]['ref_' + o] = list(r)
open('%s/%s_6_main.txt' % (R, S), 'w').write('\n'.join(lines) + '\n')
L = []
for i, nm in enumerate(names):
    L.append('== %s %g-%g' % nm); L.append('%-14s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc))
    for p in parts:
        L.append('%-14s' % p + ''.join('%11.4g' % v for v in M[p][0][i])); L.append('%-14s' % '  +-' + ''.join('%11.2g' % v for v in M[p][1][i]))
    L.append('')
open('%s/%s_6_parts.txt' % (R, S), 'w').write('\n'.join(L) + '\n')
# slopes: (N(a)-N(b))/ln(a/b), in % of the bin NLO (G6: DISENT NLO / x) or of the bin LO (Ex6, as the bin5 report)
wins = [(5e-3, 5e-4), (2e-3, 2e-4), (1e-3, 1e-4), (1e-4, 1e-5)]
unit = 'bin LO (b0 of bin5)' if S == 'Ex6' else "bin NLO (DISENT / x)"
L = ['NNLO change per unit ln tau_cut, (N(a) - N(b))/ln(a/b), a > b, in %% of the %s; +- from sqrt|ea^2-eb^2| (tau_cut points share events); (pull)' % unit,
     '%-30s ' % 'bin' + ' '.join('%-22s' % ('%.0e->%.0e' % w) for w in wins)]
sl = {}
labs = ['NNLO', 'NNLO(no d2)'] + (['NNLO_P2B', 'NNLO_P2B(bin5 extras)', 'NNLO_P2Bextra', 'NNLO_P2Bextra(bin5 extras)'] if S == 'Ex6' else [])
for lab in labs:
    m, e = comb(SPEC[lab])
    for i, nm in enumerate(names):
        den = M['b0'][0][i, 0] if S == 'Ex6' else ref[nm][0]
        cs = []; js = []
        for a, b in wins:
            ja = int(np.argmin(abs(tc - a))); jb = int(np.argmin(abs(tc - b)))
            d = m[i, ja] - m[i, jb]; ed = math.sqrt(abs(e[i, ja]**2 - e[i, jb]**2)); Lg = math.log(a/b)
            cs.append('%+7.2f+-%-6.2f(%+.1f)' % (100*d/Lg/den, 100*ed/Lg/abs(den), d/ed if ed > 0 else 0)); js.append([100*d/Lg/den, 100*ed/Lg/abs(den), d/ed if ed > 0 else 0])
        L.append('%-30s ' % ('%s %g-%g' % (lab, nm[1], nm[2])) + ' '.join('%-22s' % c for c in cs)); sl[lab + ' %g-%g' % nm[1:]] = js
J['slopes'] = sl
open('%s/%s_6_slopes.txt' % (R, S), 'w').write('\n'.join(L) + '\n')

if S == 'Ex6':   # E item 6 window (as combine5): largest tau_cut with 1e-4 <= tau_cut <= 5e-3 x lower bin edge, per bin; seed-level window sum
    edges = [n[1] for n in names[:-1]]; idx = []
    for lo in edges:
        c = [j for j in range(ntc) if tc[j] >= 1e-4 - 1e-12 and tc[j] <= 5e-3*lo + 1e-12]; idx.append(max(c, key=lambda j: tc[j]))
    W = ['== E item 6 window per bin (tau_cut per bin: %s); pull a-b vs next smaller tau_cut; last pull vs NNLOJET' % ', '.join('%g-%g: %.0e' % (n[1], n[2], tc[j]) for n, j in zip(names, idx))]
    for lab in ('NNLO', 'NNLO_P2B', 'NNLO_P2B(bin5 extras)'):
        m, e = comb(SPEC[lab]); W.append(lab); rs = []
        for b, nm in enumerate(names[:-1]):
            j = idx[b]; jb = j + 1; r = ref[('NNLO',) + nm[1:]]
            pull = (m[b, j] - m[b, jb])/math.sqrt(abs(e[b, j]**2 - e[b, jb]**2))
            W.append('  bin %.2f-%.2f tau_cut %.0e: %10.5g +- %-7.3g (smaller %.0e: %10.5g +- %-7.3g pull %+.1f)  NNLOJET %.5g +- %.2g pull %+.1f'
                     % (nm[1], nm[2], tc[j], m[b, j], e[b, j], tc[jb], m[b, jb], e[b, jb], pull, r[0], r[1], (m[b, j] - r[0])/math.hypot(e[b, j], r[1])))
        ms = 0; vs = 0
        for p, w in SPEC[lab]:
            wv = np.ones(ntc)*w if np.isscalar(w) else w
            sel = np.array([[X[p][s, b, idx[b]]*wv[idx[b]] for b in range(len(idx))] for s in range(len(X[p]))]).sum(1)
            a, b_ = stat(sel); ms = ms + a; vs = vs + b_**2
        r = ref[('NNLO', 0.05, 0.5)]
        W.append('  SUM of window values: %.5g +- %.3g   NNLOJET cross %.5g +- %.2g (pull %+.1f)' % (ms, math.sqrt(vs), r[0], r[1], (ms - r[0])/math.hypot(math.sqrt(vs), r[1])))
    open('%s/%s_6_window.txt' % (R, S), 'w').write('\n'.join(W) + '\n'); print('\n'.join(W))
json.dump(J, open('%s/%s_6.json' % (R, S), 'w'), indent=1)
print('\n'.join(lines[:2])); print('\n'.join(L))
