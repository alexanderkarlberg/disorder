#!/usr/bin/env python3
"""UPDATE 16 Oct (c)/(d): per set, with the bin11 kp / kp-x (KREDO), from the per-seed outputs:
 1. NNLO slopes (N(a)-N(b))/ln(a/b) for the Q8 windows with e1 = sqrt|ea^2-eb^2|/ln and e2 = per-seed difference
    (sqrt(sum_parts var_seed(N_p(a)-N_p(b))/N_p)/ln, d2 weighted by tau_cut < 1e-4 as in the NNLO);
 2. the plateau P = mean of the NNLO (and of the NLO = b1+lo) at tau_cut 1e-4, 3e-5, 1e-5: per part the per-seed average over the three
    tau_cut (d2 with its tau_cut < 1e-4 weight, i.e. 0 at 1e-4), mean over seeds, error sd/sqrt(N), summed over parts in quadrature;
    pulls of P and of N(tau_cut) against NNLOJET where the set has a reference (Z, W, E, B: results/kp11/<set>.json; N7 W, Z: results/nnj7).
Self-check: the per-tau_cut NNLO and NLO means per bin are compared with the committed results/kp11 JSON of the set (printed, max |rel diff|).
  plateau.py [SET ...]  -> results/kp11/PLATEAU_<SET>.txt, PLATEAU.json (merged), PLATEAU_tables.md (summary tables, all sets in PLATEAU.json)
B: the archived run logs are read from $BPROD/prod/B (see combine_b4.sh), its kp11 from $P/runs/prod/B/kp11."""
import glob, json, math, os, re, sys
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'; H = os.path.dirname(os.path.abspath(__file__)); R0 = H + '/results'; R = R0 + '/kp11'
BP = os.environ.get('BPROD', P)
TP = [1e-4, 3e-5, 1e-5]; WINS = [(5e-3, 5e-4), (1e-3, 1e-4), (2e-4, 1e-4), (1e-4, 3e-5), (1e-4, 1e-5), (3e-5, 1e-5)]
def kd(p): return p + '11' if p in ('kp', 'kp-x') else p
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def rows(f):
    tc, nm, r = None, [], []
    for l in open(f):
        w = l.split()
        if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
        elif l.startswith(' ZCELL '): nm.append('%s %g-%g' % (w[1], float(w[2]), float(w[3]))); r.append([ff(x) for x in w[4:]])
        elif l.startswith(' CELL ') and tc and len(w) == len(tc) + 3: nm.append('tauzQ %g-%g' % (float(w[1]), float(w[2]))); r.append([ff(x) for x in w[3:]])
    return tc, nm, (np.array(r) if r else None)
def load(dirs):
    """per-seed array (seeds x rows x tau_cut) over every done seed of the part directories dirs; rows truncated to the common first n"""
    out, tc, nm = [], None, None
    for d in dirs:
        for s in sorted(glob.glob('%s/s*/' % d)):
            if not (os.path.exists(s + 'done') and os.path.exists(s + 'run.log')): continue
            t, n, a = rows(s + 'run.log')
            if a is None: continue
            if nm is None: tc, nm = t, n
            k = min(len(nm), len(n)); assert t == tc and n[:k] == nm[:k], s; nm = nm[:k]; out.append(a)
    k = len(nm); return tc, nm, np.array([a[:k] for a in out])
J8, J7 = P + '/ours8j', P + '/ours7j'
STD = ['b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2']
def std(d): return {p: [d + '/' + kd(p)] for p in STD}
def pd(ds, p): return [d + '/' + kd(p) for d in ds]
SETS = {'J8': (std(J8 + '/PH9'), 'J8_PH9.json', lambda b: b['NNLO'], lambda b: b['NLO']),
        'CH8q': (std(J8 + '/CH8q'), 'CH8.json', lambda b: b['NNLO q'], lambda b: b['NLO q']),
        'CH8g': (std(J8 + '/CH8g'), 'CH8.json', lambda b: b['NNLO g'], lambda b: b['NLO g']),
        'Q8f': ({**std(J8 + '/CH8q'), **{p: pd([J8 + '/CH8q', J8 + '/Q8'], p) for p in ('r', 'rcorr', 'd2')}}, None, None, None),
        'J7PH': (std(J7 + '/PH'), 'J7_PH.json', lambda b: b['NNLO'], lambda b: b['NLO']),
        'J7W': (std(J7 + '/W'), 'J7_W.json', lambda b: b['NNLO'], lambda b: b['NLO']),
        'J7vV9': (std(J7 + '/PHv9'), None, None, None),
        'J7vZ': ({**std(J7 + '/PH'), **{p: pd([J7 + '/PHv9'], p) for p in ('vi', 'kp', 'r')}}, None, None, None),
        'Z': (std(P + '/ours5/Z'), 'Z5.json', lambda b: b['NNLO'], lambda b: b['NLO']),
        'W': (std(P + '/ours5/W'), 'W5.json', lambda b: b['NNLO'], lambda b: b['NLO']),
        'E': (std(P + '/ours5/E'), 'E5.json', lambda b: b['NNLO'], lambda b: b['NLO']),
        'E_P2B_bin5': ({**std(P + '/ours5/E'), **{p: pd([P + '/ours5/E'], p) for p in ('lo-x', 'vi-x', 'kp-x', 'r-x')}}, 'E5.json', lambda b: b['NNLO_P2B'], lambda b: b['NLO_P2B']),
        'E_P2B_Ex7': ({**std(P + '/ours5/E'), **{p: pd([P + '/ours7/Ex7'], p) for p in ('lo-x', 'vi-x', 'kp-x')}, 'r-x': pd([P + '/ours7/Ex7', P + '/ours7/Ex7-pilot'], 'r-x')},
                      'Ex7_7.json', lambda b: b['NNLO_P2B [Ex7]'], lambda b: b['NLO_P2B [Ex7]']),
        'B': ({'b1': [BP + '/prod/B/b1'], 'lo': [BP + '/prod/B/lo'], 'b2': [BP + '/prod/B/b2'], 'vi': [BP + '/prod/B/vi'], 'kp': [P + '/prod/B/kp11'],
               'r': [BP + '/prod/B/r-4'], 'rcorr': [BP + '/prod/B/rcorr-4']}, 'B_4_plain.json', lambda b: b['NNLO'], lambda b: b['NLO']),
        'G6P': (std(P + '/ours6/G6P'), 'G6P_6.json', lambda b: b['NNLO'], lambda b: b['NLO']),
        'G6W': (std(P + '/ours6/G6W'), 'G6W_6.json', lambda b: b['NNLO'], lambda b: b['NLO'])}
def nnj7(S, o):
    d = '%s/nnj7/%s' % (R0, S); out = {}
    w = [float(v) for v in open('%s/%s.cross.dat' % (d, o)).read().split('\n')[1].split()]; out['total 0-0'] = (w[0]/1000, w[1]/1000)
    for ob in ('q2', 'ptavg_12', 'm12'):
        for l in open('%s/%s.%s.dat' % (d, o, ob)):
            if l.startswith('#') or not l.strip(): continue
            w = [float(v) for v in l.split()]; lo, hi = w[0], w[2]
            out['%s %g-%g' % ({'ptavg_12': 'ptavg'}.get(ob, ob), lo, hi)] = (w[3]*(hi - lo)/1000, w[4]*(hi - lo)/1000)
    return out
JF = R + '/PLATEAU.json'; JA = json.load(open(JF)) if os.path.exists(JF) else {}
for S in (sys.argv[1:] or list(SETS)):
    src, jf, gnn, gnl = SETS[S]; X = {}; tc = names = None; nseed = {}
    for p, dirs in src.items():
        t, n, a = load(dirs); X[p] = a; nseed[p] = len(a)
        if names is None: tc, names = t, n
        k = min(len(names), len(n)); assert t == tc and n[:k] == names[:k], (S, p); names = names[:k]
    nb = len(names); X = {p: a[:, :nb] for p, a in X.items()}; tca = np.array(tc); lowd = (tca < 1e-4 - 1e-12).astype(float)
    xp = [p for p in ('lo-x',) if p in X]; xnn = [p for p in ('vi-x', 'kp-x', 'r-x') if p in X]
    SPEC = {'NLO': [(p, 1.0) for p in ['b1', 'lo'] + xp], 'NNLO': [(p, 1.0) for p in ['b2', 'vi', 'kp', 'r', 'rcorr'] + xnn] + ([('d2', lowd)] if 'd2' in X else [])}
    jt = [int(np.argmin(abs(tca - t))) for t in TP]; assert all(abs(tca[j] - t) < 1e-12 for j, t in zip(jt, TP))
    O = {'tau_cut': tc, 'seeds': nseed, 'bins': {}}; L = ['%s: seeds %s' % (S, ' '.join('%s:%d' % kv for kv in nseed.items()))]
    ref = json.load(open(R + '/' + jf))['bins'] if jf else {}
    chk = []
    for i, n in enumerate(names):
        O['bins'][n] = {}
        for lab, sp in SPEC.items():
            m = sum(X[p][:, i, :].mean(0)*w for p, w in sp); e = np.sqrt(sum((X[p][:, i, :].std(0, ddof=1)*w)**2/len(X[p]) for p, w in sp))
            if n in ref:
                c = (gnn if lab == 'NNLO' else gnl)(ref[n])[0]; chk.append(max(abs(a - b)/max(abs(b), 1e-300) for a, b in zip(m, c)) if abs(np.array(c)).max() > 0 else 0)
            def wt(w, j): return w[j] if isinstance(w, np.ndarray) else w
            ys = {p: sum(X[p][:, i, j]*wt(w, j) for j in jt)/3 for p, w in sp}
            Pm = sum(y.mean() for y in ys.values()); Pe = math.sqrt(sum(y.var(ddof=1)/len(y) for y in ys.values()))
            sl = []
            for a, b in WINS:
                ja = int(np.argmin(abs(tca - a))); jb = int(np.argmin(abs(tca - b))); Lg = math.log(a/b); d = m[ja] - m[jb]
                e1 = math.sqrt(abs(e[ja]**2 - e[jb]**2))/Lg
                e2 = math.sqrt(sum((X[p][:, i, ja]*wt(w, ja) - X[p][:, i, jb]*wt(w, jb)).var(ddof=1)/len(X[p]) for p, w in sp))/Lg
                sl.append(['%.0e->%.0e' % (a, b), d/Lg, e1, e2])
            O['bins'][n][lab] = {'mean': list(m), 'err': list(e), 'P': [Pm, Pe], 'slopes': sl}
    # references
    R7 = {}
    if S in ('W', 'Z'): R7 = {o: nnj7(S, o) for o in ('NLO', 'NNLO')}
    for n, b in O['bins'].items():
        for lab in ('NLO', 'NNLO'):
            refs = {}
            if n in ref:
                for k in ('nnlojet_' + lab, 'ref_' + lab):
                    if k in ref[n]: refs['NNLOJET' if S not in ('G6P', 'G6W') else 'DISENT/x'] = ref[n][k]
            if R7 and n in R7[lab]: refs['N7'] = list(R7[lab][n])
            b[lab]['refs'] = {}
            for rn, (rv, re_) in refs.items():
                m, e = np.array(b[lab]['mean']), np.array(b[lab]['err']); Pm, Pe = b[lab]['P']
                b[lab]['refs'][rn] = {'ref': [rv, re_], 'pull_P': (Pm - rv)/math.hypot(Pe, re_), 'pulls': list((m - rv)/np.sqrt(e**2 + re_**2))}
    O['check_max_rel'] = max(chk) if chk else None
    JA[S] = O
    # text
    L.append('self-check vs results/kp11/%s: max |rel diff| of NLO/NNLO per tau_cut and bin = %s' % (jf, O['check_max_rel']))
    L.append('P = mean over tau_cut 1e-4, 3e-5, 1e-5 (per-seed average per part); slopes per unit ln: value +- e1 [e2] (value/e2); pulls vs refs at every tau_cut and of P')
    for n, b in O['bins'].items():
        L.append('== ' + n)
        for lab in ('NLO', 'NNLO'):
            x = b[lab]; L.append('  %-5s P = %.6g +- %.3g   N(1e-4, 3e-5, 1e-5) = %s' % (lab, x['P'][0], x['P'][1], ', '.join('%.5g +- %.2g' % (x['mean'][j], x['err'][j]) for j in jt)))
            for rn, r in x['refs'].items():
                L.append('        vs %s %.6g +- %.3g: pull(P) %+.2f; pulls %s' % (rn, r['ref'][0], r['ref'][1], r['pull_P'], ' '.join('%+.1f' % v for v in r['pulls'])))
        L.append('  NNLO slopes: ' + '  '.join('%s %+.4g+-%.2g[%.2g](%+.1f)' % (w, v, e1, e2, v/e2 if e2 > 0 else 0) for w, v, e1, e2 in O['bins'][n]['NNLO']['slopes']))
    open('%s/PLATEAU_%s.txt' % (R, S), 'w').write('\n'.join(L) + '\n'); print(L[0][:200]); print(L[1])
json.dump(JA, open(JF, 'w'), indent=1, default=float)
