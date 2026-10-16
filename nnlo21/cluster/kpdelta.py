#!/usr/bin/env python3
"""KREDO (UPDATE 16 Oct (b)) item 3: dkp = bin11 - bin10 kp (or kp-x) per cell, from the same seeds.
  kpdelta.py LABEL OLDDIR [LABEL OLDDIR ...]   (NEWDIR = OLDDIR + '11')  -> results/kp11/KPDELTA_<LABEL>.txt, KPDELTA.json (merged)
Rows: every ZCELL row (jet sets) or, without ZCELL rows, the RESULT-like lines 'tau_cut' table rows the run prints (see rows()).
Per cell and tau_cut: mean of the per-seed difference, error sd/sqrt(N); and its slope (d(a)-d(b))/ln(a/b) per unit ln for the
windows of the reports, error from the per-seed difference of the slope."""
import glob, json, math, os, re, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results/kp11'; os.makedirs(R, exist_ok=True)
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def rows(f):
    """ZCELL rows if present, else the lines that start with a name and carry one number per tau_cut (the tau_cut table)"""
    tc, nm, r, alt = None, [], [], []
    for l in open(f):
        w = l.split()
        if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
        elif l.startswith(' ZCELL '): nm.append('%s %g-%g' % (w[1], float(w[2]), float(w[3]))); r.append([ff(x) for x in w[4:]])
        elif tc and len(w) == len(tc) + 1 and not l.startswith(' ZCELL'):
            try: alt.append((w[0], [ff(x) for x in w[1:]]))
            except ValueError: pass
    if not r and alt:
        seen = {}
        for n, v in alt:
            k = n if n not in seen else '%s#%d' % (n, seen[n]); seen[n] = seen.get(n, 0) + 1; nm.append(k); r.append(v)
    return tc, nm, (np.array(r) if r else None)
def seeds(d):
    out = {}
    for s in glob.glob(d + '/s*/'):
        if os.path.exists(s + 'done') and os.path.exists(s + 'run.log'): out[int(re.search(r's(\d+)/$', s).group(1))] = s + 'run.log'
    return out
WINS = [(5e-3, 5e-4), (2e-3, 2e-4), (1e-3, 1e-4), (1e-4, 1e-5)]
JF = R + '/KPDELTA.json'; J = json.load(open(JF)) if os.path.exists(JF) else {}
a = sys.argv[1:]
for lab, old in zip(a[::2], a[1::2]):
    old = old.rstrip('/'); new = old + '11'; so, sn = seeds(old), seeds(new); com = sorted(set(so) & set(sn))
    D = []; tc = names = None
    for s in com:
        t1, n1, r1 = rows(so[s]); t2, n2, r2 = rows(sn[s])
        if names is None: tc, names = t1, n1
        assert t1 == t2 == tc and n1 == n2 == names, (lab, s); D.append(r2 - r1)
    D = np.array(D); m = D.mean(0); e = D.std(0, ddof=1)/math.sqrt(len(D)); tca = np.array(tc)
    L = ['dkp = bin11 - bin10 (same seeds), %s: old %s, new %s; %d common seeds (old done %d, new done %d); units of the run output' % (lab, old, new, len(com), len(so), len(sn)),
         '%-22s' % 'tau_cut' + ''.join('%12.0e' % t for t in tc)]
    J[lab] = {'tau_cut': tc, 'n': len(com), 'cells': {}}
    for i, n in enumerate(names):
        L += ['%-22s' % n + ''.join('%12.4g' % v for v in m[i]), '%-22s' % '  +-' + ''.join('%12.2g' % v for v in e[i])]
        sl = []
        for x, y in WINS:
            ja = int(np.argmin(abs(tca - x))); jb = int(np.argmin(abs(tca - y)))
            if abs(tca[ja] - x) > 1e-12 or abs(tca[jb] - y) > 1e-12: continue
            q = (D[:, i, ja] - D[:, i, jb])/math.log(x/y); sl.append(['%.0e->%.0e' % (x, y), q.mean(), q.std(ddof=1)/math.sqrt(len(q))])
        L.append('%-22s' % '  slope/ln' + '  '.join('%s %+.4g+-%.2g' % tuple(s) for s in sl))
        J[lab]['cells'][n] = {'mean': list(m[i]), 'err': list(e[i]), 'slopes': sl}
    open('%s/KPDELTA_%s.txt' % (R, lab), 'w').write('\n'.join(L) + '\n'); print('\n'.join(L[:2] + L[2:5]))
json.dump(J, open(JF, 'w'), indent=1, default=float)
