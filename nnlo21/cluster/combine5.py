#!/usr/bin/env python3
"""Combine one bin5 validation set (E: CELL tau_zQ rows; Z, W: ZCELL rows) and compare with NNLOJET.
  combine5.py SET   (E|Z|W)  -> results/SET5_{main,parts,slopes}.txt, SET5.json
Plain equal-weight means (result), 1%-trimmed means only as a diagnostic. NNLO = b2+vi+kp+r+rcorr (+d2 for tau_cut<1e-4);
NLO = b1+lo; E also P2B = plain + the -x (P2B - plain) terms. NNLOJET: results/nnj5/SET/ (combine_nnlojet.py)."""
import glob, json, math, os, re, sys
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'; H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results'
S = sys.argv[1]
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def rows(f):
    tc, nm, r = None, [], []
    for l in open(f):
        w = l.split()
        if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
        elif l.startswith(' ZCELL '): nm.append((w[1], float(w[2]), float(w[3]))); r.append([ff(x) for x in w[4:]])
        elif l.startswith(' CELL '): nm.append(('tauzQ', float(w[1]), float(w[2]))); r.append([ff(x) for x in w[3:]])
    return tc, nm, (np.array(r) if r else None)
def load(part):
    fs = [f for f in sorted(glob.glob('%s/ours5/%s/%s/s*/run.log' % (P, S, part))) if os.path.exists(os.path.dirname(f) + '/done')]
    d = [rows(f) for f in fs]; d = [x for x in d if x[2] is not None]
    return d[0][0], d[0][1], np.array([x[2] for x in d])
def stat(x):  # x: seeds x ...
    n = len(x); return x.mean(0), x.std(0, ddof=1)/math.sqrt(n)
parts = ['b0', 'b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2'] + (['lo-x', 'vi-x', 'kp-x', 'r-x'] if S == 'E' else [])
X = {}
for p in parts:
    tc, names, X[p] = load(p)
tc = np.array(tc); nb = len(names); ntc = len(tc); lowd = tc < 1e-4 - 1e-12
cnt = {p: len(X[p]) for p in parts}
M = {p: stat(X[p]) for p in parts}
def comb(spec):   # spec: list of (part, weight array over tau_cut columns or scalar); returns mean, err arrays bins x tc
    m = 0; v = 0
    for p, w in spec:
        m = m + M[p][0]*w; v = v + (M[p][1]*w)**2
    return m, np.sqrt(v)
dm = lowd.astype(float)
SPEC = {'LO': [('b0', 1)], 'NLO': [('b1', 1), ('lo', 1)],
        'NNLO': [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')] + [('d2', dm)],
        'NNLO(no d2)': [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')],
        'NNLO(no corr)': [(p, 1) for p in ('b2', 'vi', 'kp', 'r')]}
if S == 'E':
    SPEC['NLO_P2B'] = SPEC['NLO'] + [('lo-x', 1)]
    SPEC['NNLO_P2B'] = SPEC['NNLO'] + [('vi-x', 1), ('kp-x', 1), ('r-x', 1)]
    SPEC['NNLO_P2Bextra'] = [('vi-x', 1), ('kp-x', 1), ('r-x', 1)]
# NNLOJET
def refs():
    out = {}
    for o in ('LO', 'NLO', 'NNLO'):
        d = {}
        d['total'] = None
        fn = '%s/nnj5/%s/%s.cross.dat' % (R, S, o)
        w = [float(v) for v in open(fn).read().split('\n')[1].split()]
        d[('total', 0.0, 0.0) if S != 'E' else ('tauzQ', 0.05, 0.5)] = (w[0]/1000, w[1]/1000)
        obs = ['dis_tauzq'] if S == 'E' else ['q2', 'ptavg_12', 'm12']
        sx = 0; sv = 0
        for ob in obs:
            for l in open('%s/nnj5/%s/%s.%s.dat' % (R, S, o, ob)):
                if l.startswith('#') or not l.strip(): continue
                w = [float(v) for v in l.split()]; lo, hi = w[0], w[2]
                name = {'dis_tauzq': 'tauzQ', 'ptavg_12': 'ptavg'}.get(ob, ob)
                d[(name, lo, hi)] = (w[3]*(hi - lo)/1000, w[4]*(hi - lo)/1000)
                if S == 'E': sx += w[3]*(hi - lo)/1000; sv += (w[4]*(hi - lo)/1000)**2
        if S == 'E':   # the total of the cells: sum of the bins (errors independent per bin, uncorrelated parts) vs cross
            d['sumbins'] = (sx, math.sqrt(sv))
        out[o] = d
    return out
ref = refs()
def key(nm): return nm
lines = []; J = {'tau_cut': list(tc), 'bins': {}}
def pr(s=''): lines.append(s)
hdr = '%-14s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc)
pr('SET %s: seeds per part: ' % S + ' '.join('%s:%d' % (p, cnt[p]) for p in parts))
for i, nm in enumerate(names):
    k = '%s %g-%g' % nm; J['bins'][k] = {}
    pr(); pr('== ' + k); pr(hdr)
    for lab, sp in SPEC.items():
        m, e = comb(sp)
        pr('%-14s' % lab + ''.join('%11.5g' % v for v in m[i])); pr('%-14s' % '  +-' + ''.join('%11.2g' % v for v in e[i]))
        J['bins'][k][lab] = [list(m[i]), list(e[i])]
        o = lab.split('_')[0].split('(')[0]
        if lab in ('LO', 'NLO', 'NNLO', 'NLO_P2B', 'NNLO_P2B') and o in ref and nm in ref[o]:
            r = ref[o][nm]
            pr('%-14s %.5g +- %.2g  pulls: ' % ('NNLOJET ' + o, r[0], r[1]) + ' '.join('%+.1f' % ((m[i][c] - r[0])/math.hypot(e[i][c], r[1])) for c in range(ntc)))
            J['bins'][k]['nnlojet_' + o] = list(r)
# E: per-bin window value
win = {}
if S == 'E':
    edges = [n[1] for n in names[:-1]]
    idx = []
    for lo in edges:
        c = [j for j in range(ntc) if tc[j] >= 1e-4 - 1e-12 and tc[j] <= 5e-3*lo + 1e-12]
        idx.append(max(c, key=lambda j: tc[j]))
    for lab in ('NNLO', 'NNLO_P2B'):
        sp = SPEC[lab]
        # seed-level window sum per part for a proper total error
        ms = 0; vs = 0
        rowsw = []
        for p, w in sp:
            wv = np.ones(ntc)*w if np.isscalar(w) else w
            sel = np.array([[X[p][s, b, idx[b]]*wv[idx[b]] for b in range(len(idx))] for s in range(len(X[p]))])
            mm, ee = stat(sel); ms = ms + mm; vs = vs + ee**2
            tots = sel.sum(1); a, b_ = stat(tots)
            rowsw.append((a, b_))
        W = (ms, np.sqrt(vs), sum(a for a, b in rowsw), math.sqrt(sum(b*b for a, b in rowsw)))
        win[lab] = W
    pr(); pr('== E item 6: per-bin value at the largest tau_cut with 1e-4 <= tau_cut <= 5e-3 x (lower edge); pull vs next smaller tau_cut sqrt|sa^2-sb^2|')
    for lab in ('NNLO', 'NNLO_P2B'):
        m, e = comb(SPEC[lab]); pr(lab); tot = 0
        for b, nm in enumerate(names[:-1]):
            j = idx[b]; jb = j + 1
            pull = (m[b, j] - m[b, jb])/math.sqrt(abs(e[b, j]**2 - e[b, jb]**2)) if e[b, jb] != e[b, j] else float('nan')
            o = lab.split('_')[0]; r = ref[o][nm]
            pr('  bin %.2f-%.2f tau_cut %.0e: %10.5g +- %-7.3g (smaller %.0e: %10.5g +- %-7.3g pull %+.1f)  NNLOJET %.5g +- %.2g pull %+.1f'
               % (nm[1], nm[2], tc[j], m[b, j], e[b, j], tc[jb], m[b, jb], e[b, jb], pull, r[0], r[1], (m[b, j] - r[0])/math.hypot(e[b, j], r[1])))
        r = ref[lab.split('_')[0]]['sumbins']; v, ev = win[lab][2], win[lab][3]
        pr('  SUM of window values: %.5g +- %.3g  NNLOJET sum of bins %.5g +- %.2g (pull %+.1f); NNLOJET cross %.5g +- %.2g' % (v, ev, r[0], r[1], (v - r[0])/math.hypot(ev, r[1]), *ref[lab.split('_')[0]][('tauzQ', 0.05, 0.5)]))
        J['window_' + lab] = {'sum': [v, ev], 'bins': [[float(m[b, idx[b]]), float(e[b, idx[b]]), float(tc[idx[b]])] for b in range(len(idx))]}
open('%s/%s5_main.txt' % (R, S), 'w').write('\n'.join(lines) + '\n')
# item 8: parts, slopes
L = []
for i, nm in enumerate(names):
    L.append('== %s %g-%g' % nm); L.append('%-14s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc))
    for p in parts:
        L.append('%-14s' % p + ''.join('%11.4g' % v for v in M[p][0][i])); L.append('%-14s' % '  +-' + ''.join('%11.2g' % v for v in M[p][1][i]))
    L.append('')
open('%s/%s5_parts.txt' % (R, S), 'w').write('\n'.join(L) + '\n')
wins = [(5e-3, 5e-4), (2e-3, 2e-4), (1e-3, 1e-4)]
L = ['NNLO change per unit ln tau_cut, (N(a) - N(b))/ln(a/b) with a > b, in %% of the bin LO (b0); pull = (N(a)-N(b))/sqrt|ea^2-eb^2|; lines: NNLO (b2+vi+kp+r+rcorr) and also P2B for E',
     '%-22s %12s %12s %12s' % ('bin', '5e-3->5e-4', '2e-3->2e-4', '1e-3->1e-4')]
sl = {}
labs = ['NNLO'] + (['NNLO_P2B'] if S == 'E' else [])
for lab in labs:
    m, e = comb(SPEC[lab]); lo = M['b0'][0]
    for i, nm in enumerate(names):
        cs = []; js = []
        for a, b in wins:
            ja = int(np.argmin(abs(tc - a))); jb = int(np.argmin(abs(tc - b)))
            d = m[i, ja] - m[i, jb]; ed = math.sqrt(abs(e[i, ja]**2 - e[i, jb]**2)); L_ = math.log(a/b)
            cs.append('%+6.2f+-%-5.2f(%+.1f)' % (100*d/L_/lo[i, 0], 100*ed/L_/abs(lo[i, 0]), d/ed if ed > 0 else 0)); js.append([100*d/L_/lo[i, 0], 100*ed/L_/abs(lo[i, 0])])
        L.append('%-22s ' % ('%s %s %g-%g' % ((lab[5:] or 'NNLO',) + nm)) + ' '.join(cs)); sl[lab + ' %s %g-%g' % nm] = js
J['slopes'] = sl
open('%s/%s5_slopes.txt' % (R, S), 'w').write('\n'.join(L) + '\n')
json.dump(J, open('%s/%s5.json' % (R, S), 'w'), indent=1)
print('\n'.join(L[:20]))
