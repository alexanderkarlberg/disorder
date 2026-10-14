#!/usr/bin/env python3
"""Ex7 (14 Oct), UPDATE 13 Oct (e) item 2 / UPDATE (after REPORT-13) item 3. E set (CELL tau_zQ rows), P2B extras with ZFIX=3 P2BSLICE=1 P2BEXTRA=1.
  combine7e.py  -> results/Ex7_7_{main,parts,slopes,window,scatter,consistency}.txt, Ex7_7.json
Plain parts (b0, b1, lo, b2, vi, kp, r, rcorr, d2): bin5 runs ours5/E (as combine6). Extras (lo-x, vi-x, kp-x, r-x) from
  Ex6  = ours6/Ex6 (VTARGET=0; 200/250/200/2000 seeds), Ex7 = ours7/Ex7 (VTARGET=-1, seeds 5001-; r-x 2000) + ours7/Ex7-pilot/r-x (9101-9300) = 2200 r-x seeds,
  bin5 = ours5/E extras (lo-x 50, vi-x 62, kp-x 50, r-x 875).
Variants of the extras (NNLO_P2B = NNLO + vi-x + kp-x + r-x, NLO_P2B = NLO + lo-x; NNLO_P2Bextra = vi-x+kp-x+r-x):
  Ex7        : Ex7 alone, all four parts
  Ex6+Ex7    : per cell (bin x tau_cut) inverse-variance mean of the Ex6 and Ex7 part means (error 1/sqrt(sum 1/e^2)), per part
  Ex7+b5lvk  : r-x from Ex7 alone; lo-x, vi-x, kp-x = per-cell inverse-variance mean of Ex7 and the bin5 extras
  Ex6, bin5  : for comparison (as combine6's 'NNLO_P2B' and 'bin5 extras')
Plain equal-weight seed means per part (error = scatter/sqrt(N)), parts in quadrature. NNLOJET E: results/nnj5/E (unchanged by N7)."""
import glob, json, math, os, re, sys
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'; H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results'
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
def load(dirs, part, tag):
    fs = []
    for d in dirs: fs += [x + 'run.log' for x in sorted(glob.glob('%s/%s/s*/' % (d, part))) if os.path.exists(x + 'done') and os.path.exists(x + 'run.log')]
    dd = [rows(f) for f in fs]; ok = [x for x in dd if x[2] is not None]
    arr = np.array([x[2] for x in ok]); nf = int((~np.isfinite(arr).all(axis=(1, 2))).sum())
    INFO[tag] = (len(fs), len(ok), nf)
    return ok[0][0], ok[0][1], arr
E5 = P + '/ours5/E'; X = {}
for p in ('b0', 'b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2'): tc, names, X[p] = load([E5], p, 'bin5 ' + p)
XP = ('lo-x', 'vi-x', 'kp-x', 'r-x'); SRC = {}
for p in XP:
    _, _, SRC[('5', p)] = load([E5], p, 'bin5 ' + p)
    _, _, SRC[('6', p)] = load([P + '/ours6/Ex6'], p, 'Ex6 ' + p)
    _, _, SRC[('7', p)] = load([P + '/ours7/Ex7'] + ([P + '/ours7/Ex7-pilot'] if p == 'r-x' else []), p, 'Ex7 ' + p)
tc = np.array(tc); ntc = len(tc); nb = len(names); lowd = (tc < 1e-4 - 1e-12).astype(float)
def stat(x): n = len(x); return x.mean(0), x.std(0, ddof=1)/math.sqrt(n)
def ivar(a, b):
    (ma, ea), (mb, eb) = a, b
    wa = np.where(ea > 0, 1/np.where(ea > 0, ea, 1)**2, 0); wb = np.where(eb > 0, 1/np.where(eb > 0, eb, 1)**2, 0); s = wa + wb
    m = np.where(s > 0, (ma*wa + mb*wb)/np.where(s > 0, s, 1), (ma + mb)/2); e = np.where(s > 0, 1/np.sqrt(np.where(s > 0, s, 1)), 0)
    return m, e
M = {p: stat(X[p]) for p in X}; S = {k: stat(v) for k, v in SRC.items()}
VAR = {'Ex7': {p: S[('7', p)] for p in XP},
       'Ex6+Ex7': {p: ivar(S[('6', p)], S[('7', p)]) for p in XP},
       'Ex7+b5lvk': {**{p: ivar(S[('7', p)], S[('5', p)]) for p in ('lo-x', 'vi-x', 'kp-x')}, 'r-x': S[('7', 'r-x')]},
       'Ex6': {p: S[('6', p)] for p in XP}, 'bin5': {p: S[('5', p)] for p in XP}}
def comb(spec):
    m = 0; v = 0
    for mm, ee, w in spec: m = m + mm*w; v = v + (ee*w)**2
    return m, np.sqrt(v)
base = lambda ps: [(M[p][0], M[p][1], 1) for p in ps]
NNLO = base(('b2', 'vi', 'kp', 'r', 'rcorr')) + [(M['d2'][0], M['d2'][1], lowd)]; NLO = base(('b1', 'lo'))
def variant(v):
    ex = VAR[v]; sp = {}
    sp['NLO_P2B'] = NLO + [(*ex['lo-x'], 1)]
    sp['NNLO_P2B'] = NNLO + [(*ex[p], 1) for p in ('vi-x', 'kp-x', 'r-x')]
    sp['NNLO_P2Bextra'] = [(*ex[p], 1) for p in ('vi-x', 'kp-x', 'r-x')]
    return sp
SPEC = {'LO': base(('b0',)), 'NLO': NLO, 'NNLO': NNLO}
for v in VAR:
    for k, s in variant(v).items(): SPEC['%s [%s]' % (k, v)] = s
N = {}
for o in ('LO', 'NLO', 'NNLO'):
    fn = '%s/nnj5/E/%s.cross.dat' % (R, o); w = [float(v) for v in open(fn).read().split('\n')[1].split()]; N[(o, 0.05, 0.5)] = (w[0]/1000, w[1]/1000)
    for l in open('%s/nnj5/E/%s.dis_tauzq.dat' % (R, o)):
        if l.startswith('#') or not l.strip(): continue
        w = [float(v) for v in l.split()]; N[(o, w[0], w[2])] = (w[3]*(w[2] - w[0])/1000, w[4]*(w[2] - w[0])/1000)
lines = []; J = {'tau_cut': list(tc), 'seeds': {'%s %s' % k: len(v) for k, v in SRC.items()}, 'bins': {}}
def pr(s=''): lines.append(s)
hdr = '%-30s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc)
pr('Ex7 report: seeds per extra (set part): ' + ' '.join('%s %s:%d' % (k[0], k[1], len(v)) for k, v in SRC.items()) + '; plain parts (bin5 E): ' + ' '.join('%s:%d' % (p, len(X[p])) for p in X))
pr('seed bookkeeping (dirs, parsed, with non-finite cells): ' + '; '.join('%s %s' % (k, v) for k, v in INFO.items()))
for i, nm in enumerate(names):
    k = '%s %g-%g' % nm; J['bins'][k] = {}; pr(); pr('== ' + k); pr(hdr)
    for lab, sp in SPEC.items():
        m, e = comb(sp); pr('%-30s' % lab + ''.join('%11.5g' % v for v in m[i])); pr('%-30s' % '  +-' + ''.join('%11.2g' % v for v in e[i]))
        J['bins'][k][lab] = [list(m[i]), list(e[i])]
        o = lab.split('_')[0].split(' ')[0]; r = N.get((o,) + nm[1:]) if (lab in ('LO', 'NLO', 'NNLO') or lab.startswith(('NLO_P2B [', 'NNLO_P2B ['))) else None
        if r:
            pr('%-30s %.5g +- %.2g  pulls: ' % ('  NNLOJET ' + o, r[0], r[1]) + ' '.join('%+.1f' % ((m[i][c] - r[0])/math.hypot(e[i][c], r[1])) for c in range(ntc)))
            J['bins'][k]['ref_' + o] = list(r)
open('%s/Ex7_7_main.txt' % R, 'w').write('\n'.join(lines) + '\n')
Lp = []
for i, nm in enumerate(names):
    Lp += ['== %s %g-%g' % nm, '%-14s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc)]
    for p in X: Lp += ['%-14s' % p + ''.join('%11.4g' % v for v in M[p][0][i]), '%-14s' % '  +-' + ''.join('%11.2g' % v for v in M[p][1][i])]
    for v in VAR:
        for p in XP: Lp += ['%-14s' % ('%s[%s]' % (p, v)) + ''.join('%11.4g' % q for q in VAR[v][p][0][i]), '%-14s' % '  +-' + ''.join('%11.2g' % q for q in VAR[v][p][1][i])]
    Lp.append('')
open('%s/Ex7_7_parts.txt' % R, 'w').write('\n'.join(Lp) + '\n')
wins = [(5e-3, 5e-4), (2e-3, 2e-4), (1e-3, 1e-4), (1e-4, 1e-5)]
L = ['NNLO change per unit ln tau_cut, (N(a)-N(b))/ln(a/b), a > b, in % of the bin LO (b0 of bin5); +- from sqrt|ea^2-eb^2| (tau_cut points share events); (pull)',
     '%-42s ' % 'bin' + ' '.join('%-24s' % ('%.0e->%.0e' % w) for w in wins)]
sl = {}
labs = ['NNLO'] + ['NNLO_P2B [%s]' % v for v in VAR] + ['NNLO_P2Bextra [%s]' % v for v in VAR]
for lab in labs:
    m, e = comb(SPEC[lab])
    for i, nm in enumerate(names):
        den = M['b0'][0][i, 0]; cs = []; js = []
        for a, b in wins:
            ja = int(np.argmin(abs(tc - a))); jb = int(np.argmin(abs(tc - b))); d = m[i, ja] - m[i, jb]; ed = math.sqrt(abs(e[i, ja]**2 - e[i, jb]**2)); Lg = math.log(a/b)
            cs.append('%+8.2f+-%-6.2f(%+.1f)' % (100*d/Lg/den, 100*ed/Lg/abs(den), d/ed if ed > 0 else 0)); js.append([100*d/Lg/den, 100*ed/Lg/abs(den), d/ed if ed > 0 else 0])
        L.append('%-42s ' % ('%s %g-%g' % (lab, nm[1], nm[2])) + ' '.join('%-24s' % c for c in cs)); sl[lab + ' %g-%g' % nm[1:]] = js
J['slopes'] = sl
open('%s/Ex7_7_slopes.txt' % R, 'w').write('\n'.join(L) + '\n')
# E item 6 window
edges = [n[1] for n in names[:-1]]; idx = []
for lo in edges:
    c = [j for j in range(ntc) if tc[j] >= 1e-4 - 1e-12 and tc[j] <= 5e-3*lo + 1e-12]; idx.append(max(c, key=lambda j: tc[j]))
W = ['== E item 6 window per bin (tau_cut per bin: %s); pull a-b vs next smaller tau_cut; last pull vs NNLOJET' % ', '.join('%g-%g: %.0e' % (n[1], n[2], tc[j]) for n, j in zip(names, idx)),
     'SUM of window values: seed-level sum over the five bins per part (per set), the two sets of a combined part joined by inverse variance of that sum']
def wsum(arr, w): return stat(np.array([[arr[s, b, idx[b]]*w[idx[b]] for b in range(len(idx))] for s in range(len(arr))]).sum(1))
def wsum_part(p, w=1.0): return wsum(X[p], np.ones(ntc)*w)
SETS = {'Ex7': {p: [('7',)] for p in XP}, 'Ex6+Ex7': {p: [('6',), ('7',)] for p in XP}, 'Ex7+b5lvk': {**{p: [('7',), ('5',)] for p in ('lo-x', 'vi-x', 'kp-x')}, 'r-x': [('7',)]},
        'Ex6': {p: [('6',)] for p in XP}, 'bin5': {p: [('5',)] for p in XP}}
def xsum(v, p):
    parts_ = [wsum(SRC[(s[0], p)], np.ones(ntc)) for s in SETS[v][p]]
    if len(parts_) == 1: return parts_[0]
    a, b = parts_; wa, wb = 1/a[1]**2, 1/b[1]**2; return (a[0]*wa + b[0]*wb)/(wa + wb), 1/math.sqrt(wa + wb)
for lab in ['NNLO'] + ['NNLO_P2B [%s]' % v for v in VAR]:
    m, e = comb(SPEC[lab]); W.append(lab); o = 'NNLO'
    for b, nm in enumerate(names[:-1]):
        j = idx[b]; jb = j + 1; r = N[(o,) + nm[1:]]; pull = (m[b, j] - m[b, jb])/math.sqrt(abs(e[b, j]**2 - e[b, jb]**2))
        W.append('  bin %.2f-%.2f tau_cut %.0e: %10.5g +- %-7.3g (smaller %.0e: %10.5g +- %-7.3g pull %+.1f)  NNLOJET %.5g +- %.2g pull %+.1f'
                 % (nm[1], nm[2], tc[j], m[b, j], e[b, j], tc[jb], m[b, jb], e[b, jb], pull, r[0], r[1], (m[b, j] - r[0])/math.hypot(e[b, j], r[1])))
    ms = 0; vs = 0
    for p, w in [(q, 1.0) for q in ('b2', 'vi', 'kp', 'r', 'rcorr')] + [('d2', lowd)]:
        a, b_ = wsum(X[p], np.ones(ntc)*w); ms += a; vs += b_**2
    if '[' in lab:
        v = lab.split('[')[1][:-1]
        for p in ('vi-x', 'kp-x', 'r-x'): a, b_ = xsum(v, p); ms += a; vs += b_**2
    r = N[('NNLO', 0.05, 0.5)]
    W.append('  SUM of window values: %.5g +- %.3g   NNLOJET cross %.5g +- %.2g (pull %+.1f)' % (ms, math.sqrt(vs), r[0], r[1], (ms - r[0])/math.hypot(math.sqrt(vs), r[1])))
    J.setdefault('window', {})[lab] = [ms, math.sqrt(vs)]
open('%s/Ex7_7_window.txt' % R, 'w').write('\n'.join(W) + '\n')
# per-seed scatter: sd and 1.48 MAD per part and tau_cut (total row 0.05-0.5), bin5, Ex6, Ex7 side by side; outlier fraction on final per-seed values
Sc = ['Per-seed scatter of the total row (bin 0.05-0.5), pb: sd, 1.48 MAD, median, fraction of seeds with |x - median| > 20 x 1.48 MAD (final per-seed values, not iteration lines)',
      '%-8s %-5s %6s' % ('part', 'set', 'N') + ''.join('| tau_cut %-5.0e sd / 1.48MAD / median / frac' % t if False else '' for t in []) ]
tot = nb - 1
for p in XP + ('r', 'rcorr', 'vi', 'kp', 'lo'):
    for tag, arr in ((('bin5', SRC[('5', p)]), ('Ex6', SRC[('6', p)]), ('Ex7', SRC[('7', p)])) if p in XP else (('bin5', X[p]),)):
        for j, t in enumerate(tc):
            if t not in (1e-3, 1e-4, 1e-5): continue
            x = arr[:, tot, j]; md = np.median(x); mad = 1.4826*np.median(np.abs(x - md))
            Sc.append('%-6s %-5s N=%5d tau_cut %.0e  sd %10.4g  1.48MAD %10.4g  median %10.4g  frac>20MAD %.4f  max|x-med|/1.48MAD %.3g' % (p, tag, len(x), t, x.std(ddof=1), mad, md, np.mean(np.abs(x - md) > 20*mad) if mad > 0 else float('nan'), np.abs(x - md).max()/mad if mad > 0 else float('nan')))
open('%s/Ex7_7_scatter.txt' % R, 'w').write('\n'.join(Sc) + '\n')
# consistency of the three sources of each extra part (same integrand): total row, difference of means and pull sqrt(e1^2+e2^2)
Cn = ['Consistency of the extra parts, total row (bin 0.05-0.5), pb: Ex7 - Ex6 and Ex7 - bin5 (diff +- sqrt(e1^2+e2^2), pull); errors are seed scatter/sqrt(N) (Ex6 and bin5 r-x/vi-x errors are tail-driven)',
      '%-6s %-9s ' % ('part', 'compare') + ''.join('%-24s' % ('%.0e' % t) for t in tc)]
for p in XP:
    for o, key in (('Ex6', '6'), ('bin5', '5')):
        a, b = S[('7', p)], S[(key, p)]; cs = []
        for j in range(ntc):
            d = a[0][tot, j] - b[0][tot, j]; ee = math.hypot(a[1][tot, j], b[1][tot, j]); cs.append('%+9.4g+-%-7.3g(%+.1f)' % (d, ee, d/ee if ee > 0 else 0))
        Cn.append('%-6s Ex7-%-5s ' % (p, o) + ''.join('%-24s' % c for c in cs))
open('%s/Ex7_7_consistency.txt' % R, 'w').write('\n'.join(Cn) + '\n')
json.dump(J, open('%s/Ex7_7.json' % R, 'w'), indent=1)
print('\n'.join(W)); print('\n'.join(Sc[:3]))
