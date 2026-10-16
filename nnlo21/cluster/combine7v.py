#!/usr/bin/env python3
"""J7v (15 Oct): J7 photon, bin9 + ZCUTS=1, seeds 3001- ($P/runs/ours7j/PHv9: vi, kp, r with VTARGET=-1 VEGAS_EQUAL=2; b2, rcorr, d2, b1, lo reruns) next to J7 PH (bin6, seeds 1001-).
  combine7v.py -> results/J7v_{main,main_zeus,parts,slopes,consistency,scatter}.txt, J7v.json, results/j7/PHv9_<part>.zcell.gz
Sets: PH = old J7 PH; V9 = all parts from PHv9 (ZEUS rows and the six ZCUTS rows et4/noeta/et4ne x total/ptavg 8-15); Z = J7v of UPDATE (c): vi, kp, r from PHv9, b2/rcorr/d2/b1/lo of J7 PH (ZEUS rows only).
Plain equal-weight means over done seeds, error = seed scatter/sqrt(N), parts in quadrature; slopes and pulls as combine7j.py.
Consistency: per part (vi, kp, r), bin, tau_cut: J7v - J7 PH +- sqrt(e1^2+e2^2) (pull). Scatter: per-seed sd, 1.48 MAD, median, fraction beyond 20 x 1.48 MAD."""
import glob, gzip, json, math, os, re
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs/ours7j'; H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results'
K11 = os.environ.get('KP11') == '1'  # KREDO: read kp11 / kp-x11 (bin11) instead of kp / kp-x, write to results/kp11
R0 = R  # inputs (NNLOJET references, earlier results) stay in results/
if K11: R = R + '/kp11'; os.makedirs(R, exist_ok=True)
def kd(p): return p + '11' if K11 and p in ('kp', 'kp-x') else p
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
INFO = {}
def load(S, p, tag):
    ds = sorted(glob.glob('%s/%s/%s/s*/' % (P, S, kd(p))), key=lambda s: int(re.search(r's(\d+)/$', s).group(1)))
    dd, seeds = [], []
    for d in ds:
        if not os.path.exists(d + 'done'): continue
        x = rows(d + 'run.log')
        if x[2] is not None: dd.append(x); seeds.append(int(re.search(r's(\d+)/$', d).group(1)))
    arr = np.array([x[2] for x in dd]); INFO[tag] = (len(ds), len(dd), int((~np.isfinite(arr).all(axis=(1, 2))).sum()))
    return dd[0][0], dd[0][1], arr, seeds
PARTS = ('b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2')
A = {}; V = {}; seedsV = {}
for p in PARTS:
    tc, namesA, A[p], _ = load('PH', p, 'PH ' + p)
    tcv, names, V[p], seedsV[p] = load('PHv9', p, 'PHv9 ' + p); assert list(tcv) == list(tc) and names[:len(namesA)] == namesA
nz = len(namesA)
os.makedirs(R + '/j7', exist_ok=True)
for p in PARTS:
    with gzip.open('%s/j7/PHv9_%s.zcell.gz' % (R, p), 'wt') as f:
        f.write('# J7v part %s: raw ZCELL rows [pb] per seed (bin9, ZCUTS=1; vi/kp/r with VTARGET=-1 VEGAS_EQUAL=2; seeds 3001-)\n# tau_cut %s\n# columns: seed observable lo hi values(tau_cut columns)\n' % (p, ' '.join('%g' % t for t in tc)))
        for s, a in zip(seedsV[p], V[p]):
            for n, v in zip(names, a): f.write('%d %s %s\n' % (s, n.replace('-', ' '), ' '.join('%.8E' % q for q in v)))
tc = np.array(tc); lowd = (tc < 1e-4 - 1e-12).astype(float)
st = lambda x: (x.mean(0), x.std(0, ddof=1) / math.sqrt(len(x)))
MA = {p: st(A[p]) for p in A}; MV = {p: st(V[p]) for p in V}
# Z: vi, kp, r of the new runs (ZEUS rows) + J7 PH's old other parts
MZ = dict(MA); MZ.update({p: (MV[p][0][:nz], MV[p][1][:nz]) for p in ('vi', 'kp', 'r')})
def comb(M, spec):
    m = 0; v = 0
    for p, w in spec: m = m + M[p][0] * w; v = v + (M[p][1] * w) ** 2
    return m, np.sqrt(v)
SPEC = {'NLO': [('b1', 1), ('lo', 1)], 'NNLO': [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')] + [('d2', lowd)],
        'NNLO(no d2)': [(p, 1) for p in ('b2', 'vi', 'kp', 'r', 'rcorr')]}
C = {'PH': {k: comb(MA, v) for k, v in SPEC.items()}, 'Z': {k: comb(MZ, v) for k, v in SPEC.items()}, 'V9': {k: comb(MV, v) for k, v in SPEC.items()}}
NR = {'PH': namesA, 'Z': namesA, 'V9': names}
nlo = {k: np.abs(C[k]['NLO'][0][:, -1]) for k in C}
hdr = '%-24s' % 'tau_cut' + ''.join('%11.0e' % t for t in tc)
cnt = {k: {p: len(M[p][0]) and (len(A[p]) if k == 'PH' else len(V[p])) for p in PARTS} for k, M in (('PH', MA), ('V9', MV))}
L = ['J7v, window (x,Q2) = (0.01, 400), photon; dsigma/dx dQ2 [pb/GeV2]. Sets: PH = J7 PH (bin6); Z = J7v of UPDATE (c): vi, kp, r from the new bin9 runs, other parts from J7 PH; V9 = all parts from the new bin9 runs',
     'seeds PH: ' + ' '.join('%s:%d' % kv for kv in cnt['PH'].items()) + '; seeds V9: ' + ' '.join('%s:%d' % kv for kv in cnt['V9'].items()),
     'bookkeeping (dirs, done+parsed, non-finite): ' + '; '.join('%s %s' % kv for kv in INFO.items())]
J = {'tau_cut': list(tc), 'seeds': cnt, 'bins': {}}
def block(tag, i, labs):
    out = []
    for lab in labs:
        m, e = C[tag][lab]; out += ['%-24s' % ('%s %s' % (lab, tag)) + ''.join('%11.5g' % v for v in m[i]), '%-24s' % '  +-' + ''.join('%11.2g' % v for v in e[i])]
        J['bins'].setdefault(NR[tag][i], {})['%s %s' % (lab, tag)] = [list(m[i]), list(e[i])]
    return out
labs = ('NLO', 'NNLO', 'NNLO(no d2)')
Lz = list(L)
for i, n in enumerate(names):
    if nlo['V9'][i] > 0:
        L += ['', '== ' + n, hdr] + block('V9', i, labs)
for i, n in enumerate(namesA):
    if nlo['PH'][i] > 0:
        Lz += ['', '== ' + n + ' (ZEUS rows)', hdr] + block('Z', i, labs) + block('PH', i, labs)
open(R + '/J7v_main.txt', 'w').write('\n'.join(L) + '\n'); open(R + '/J7v_main_zeus.txt', 'w').write('\n'.join(Lz) + '\n')
Lp = []
for i, n in enumerate(names):
    if nlo['V9'][i] <= 0: continue
    Lp += ['== ' + n + ' (parts, V9 = bin9 reruns; PH = J7 PH; ZEUS rows only for PH)', hdr]
    for p in PARTS:
        Lp += ['%-24s' % (p + ' V9') + ''.join('%11.4g' % v for v in MV[p][0][i]), '%-24s' % '  +-' + ''.join('%11.2g' % v for v in MV[p][1][i])]
        if i < nz: Lp += ['%-24s' % (p + ' PH') + ''.join('%11.4g' % v for v in MA[p][0][i]), '%-24s' % '  +-' + ''.join('%11.2g' % v for v in MA[p][1][i])]
    Lp.append('')
open(R + '/J7v_parts.txt', 'w').write('\n'.join(Lp) + '\n')
wins = [(5e-3, 5e-4), (2e-3, 2e-4), (1e-3, 1e-4), (1e-4, 1e-5)]
Ls = ['NNLO change per unit ln tau_cut, (N(a)-N(b))/ln(a/b), a > b, in % of |NLO| of the same row at tau_cut 1e-5 (NLO of the same set); +- sqrt|ea^2-eb^2|/ln; (pull). Sets PH, Z, V9 as in J7v_main.txt',
      '%-40s ' % 'row' + ' '.join('%-24s' % ('%.0e->%.0e' % w) for w in wins)]
sl = {}
for lab in ('NNLO', 'NNLO(no d2)'):
    for tag in ('PH', 'Z', 'V9'):
        m, e = C[tag][lab]
        for i, n in enumerate(NR[tag]):
            if nlo[tag][i] <= 0 or (tag == 'V9' and False): continue
            cs = []; js = []
            for a, b in wins:
                ja = int(np.argmin(abs(tc - a))); jb = int(np.argmin(abs(tc - b))); d = m[i, ja] - m[i, jb]; ed = math.sqrt(abs(e[i, ja] ** 2 - e[i, jb] ** 2)); Lg = math.log(a / b)
                v = [100 * d / Lg / nlo[tag][i], 100 * ed / Lg / nlo[tag][i], d / ed if ed > 0 else 0]; js.append(v); cs.append('%+8.2f+-%-6.2f(%+.1f)' % tuple(v))
            Ls.append('%-40s ' % ('%s %s %s' % (lab, n, tag)) + ' '.join('%-24s' % c for c in cs)); sl['%s %s %s' % (lab, n, tag)] = js
J['slopes'] = sl
open(R + '/J7v_slopes.txt', 'w').write('\n'.join(Ls) + '\n')
# consistency: new bin9 part means minus J7 PH, ZEUS rows
Lc = ['New (bin9) minus J7 PH per part, ZEUS rows; part values x 1e6 [pb/GeV2]; diff +- sqrt(e1^2+e2^2) (pull). vi, kp, r: J7v (VTARGET=-1) - J7 PH; b2, rcorr, d2, b1, lo: reruns - J7 PH (same integrand, ZCUTS only adds rows); NNLO: Z - PH and V9 - PH',
      '%-30s' % 'row part' + ''.join('%-22s' % ('%.0e' % t) for t in tc)]
cons = {}
def pc(i, a, b, lab):
    (m1, e1), (m2, e2) = a, b
    d = (m1[i] - m2[i]) * 1e6; ed = np.hypot(e1[i], e2[i]) * 1e6; pu = np.where(ed > 0, d / np.where(ed > 0, ed, 1), 0)
    Lc.append('%-30s' % lab + ''.join('%-22s' % ('%+.3g+-%.2g(%+.1f)' % (x, y, z)) for x, y, z in zip(d, ed, pu))); cons[lab] = [list(d), list(ed), list(pu)]
for i, n in enumerate(namesA):
    if nlo['PH'][i] <= 0: continue
    for p in ('vi', 'kp', 'r', 'b2', 'rcorr', 'd2', 'b1', 'lo'): pc(i, MV[p][:2] if False else (MV[p][0][:nz], MV[p][1][:nz]), MA[p], '%s %s' % (n, p))
    pc(i, C['Z']['NNLO'], C['PH']['NNLO'], '%s NNLO(Z-PH)' % n); pc(i, (C['V9']['NNLO'][0][:nz], C['V9']['NNLO'][1][:nz]), C['PH']['NNLO'], '%s NNLO(V9-PH)' % n)
J['consistency'] = cons
open(R + '/J7v_consistency.txt', 'w').write('\n'.join(Lc) + '\n')
Lsc = ['Per-seed scatter, units 1e-6 pb/GeV2 (value per seed of the part): sd, 1.48 MAD, median, fraction with |x - median| > 20 x 1.48 MAD; N; mean +- sd/sqrtN',
       '%-6s %-9s %-10s %-7s %5s %12s %10s %10s %10s %10s %9s' % ('part', 'set', 'cell', 'tau_cut', 'N', 'mean', 'sd/sqrtN', 'sd', '1.48MAD', 'median', 'frac>20M')]
sc = {}
for p in ('vi', 'kp', 'r'):
    for tag, D in (('J7v', V[p]), ('J7 PH', A[p])):
        for cell in ('total 0-0', 'ptavg 8-15'):
            ic = names.index(cell)
            for t in (1e-3, 1e-4):
                j = int(np.argmin(abs(tc - t))); x = D[:, ic, j] * 1e6; n = len(x); md = np.median(x); sd = x.std(ddof=1); mad = 1.4826 * np.median(np.abs(x - md)); fr = np.mean(np.abs(x - md) > 20 * mad)
                Lsc.append('%-6s %-9s %-10s %-7g %5d %12.5g %10.4g %10.4g %10.4g %10.4g %9.4f' % (p, tag, cell, t, n, x.mean(), sd / math.sqrt(n), sd, mad, md, fr))
                sc['%s %s %s %g' % (p, tag, cell, t)] = dict(n=n, mean=x.mean(), err=sd / math.sqrt(n), sd=sd, mad=mad, med=md, frac=fr)
J['scatter'] = sc
open(R + '/J7v_scatter.txt', 'w').write('\n'.join(Lsc) + '\n')
json.dump(J, open(R + '/J7v.json', 'w'), indent=1)
print('\n'.join(L[:3])); print('\n'.join(Lsc))
