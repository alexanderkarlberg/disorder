#!/usr/bin/env python3
"""X8 (UPDATE 14 Oct item 2): total row (CELL 0.05 0.50) of the extra parts lo-x, vi-x, kp-x, r-x at tau_cut 2e-2 ... 2e-4:
X8 (ours7/X8, VTARGET=53, seeds 8001-), Ex7 (ours7/Ex7 5001- + ours7/Ex7-pilot r-x 9101-9300, VTARGET=-1), Ex6 (ours6/Ex6, VTARGET=0),
bin5 (ours5/E extras). Plain equal-weight means, error = sd/sqrt(N); pulls = difference / sqrt(e1^2 + e2^2).
  x8.py  -> results/X8_total.txt (numbers), results/X8_total.json, REPORT table on stdout (markdown)"""
import glob, json, math, os, re
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'; H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results'
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def total(f):
    tc, v = None, None
    for l in open(f):
        w = l.split()
        if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
        elif l.startswith(' CELL ') and w[1] == '0.05' and w[2] == '0.50': v = [ff(x) for x in w[3:]]
    return tc, v
def load(dirs, part):
    a, n_dir, tc = [], 0, None
    for d in dirs:
        for s in sorted(glob.glob('%s/%s/s*/' % (d, part))):
            if not os.path.exists(s + 'done'): continue
            n_dir += 1; t, v = total(s + 'run.log')
            if v is None: continue
            tc = t; a.append(v)
    return tc, (np.array(a) if a else np.zeros((0, 1))), n_dir
SETS = {'X8': [P + '/ours7/X8'], 'Ex7': [P + '/ours7/Ex7', P + '/ours7/Ex7-pilot'], 'Ex6': [P + '/ours6/Ex6'], 'bin5': [P + '/ours5/E']}
XP = ('lo-x', 'vi-x', 'kp-x', 'r-x'); D = {}; TC = None; info = []
for s, dirs in SETS.items():
    for p in XP:
        tc, a, nd = load(dirs if not (s == 'Ex7' and p != 'r-x') else dirs[:1], p)
        if tc is None: tc = TC
        TC = TC or tc; assert tc == TC
        D[(s, p)] = a; info.append('%s %s: %d done dirs, %d parsed, %d with non-finite cell' % (s, p, nd, len(a), int((~np.isfinite(a).all(1)).sum())))
cols = [j for j, t in enumerate(TC) if 2e-4 - 1e-12 <= t <= 2e-2 + 1e-12]
def st(x):
    n = len(x); md = np.median(x); sd = x.std(ddof=1); return dict(n=n, mean=x.mean(), err=sd/math.sqrt(n), sd=sd, mad=1.4826*np.median(np.abs(x - md)), med=md)
J = {'tau_cut': [TC[j] for j in cols], 'stats': {}, 'pulls': {}}; L = []; M = []
for p in XP:
    J['stats'][p] = {}; J['pulls'][p] = {}
    M += ['', '#### %s (total row 0.05-0.50, pb)' % p, '', '| tau_cut | set | N | mean +- sd/sqrtN | per-seed sd | 1.48 MAD | median |', '|---|---|---|---|---|---|---|']
    for j in cols:
        t = TC[j]; S = {s: st(D[(s, p)][:, j]) for s in SETS}; J['stats'][p]['%g' % t] = S
        for s, r in S.items():
            M.append('| %.0e | %s | %d | %.5g +- %.3g | %.4g | %.4g | %.4g |' % (t, s, r['n'], r['mean'], r['err'], r['sd'], r['mad'], r['med']))
            L.append('%-5s tau_cut %.0e %-4s N=%4d mean %12.6g err %10.4g sd %11.5g 1.48MAD %10.4g median %11.5g' % (p, t, s, r['n'], r['mean'], r['err'], r['sd'], r['mad'], r['med']))
        J['pulls'][p]['%g' % t] = {o: (S['X8']['mean'] - S[o]['mean'])/math.hypot(S['X8']['err'], S[o]['err']) for o in ('Ex7', 'Ex6', 'bin5')}
    M += ['', '| pull (X8 - other)/sqrt(e1^2+e2^2) | ' + ' | '.join('%.0e' % TC[j] for j in cols) + ' |', '|---|' + '---|'*len(cols)]
    for o in ('Ex7', 'Ex6', 'bin5'): M.append('| X8 - %s | ' % o + ' | '.join('%+.1f' % J['pulls'][p]['%g' % TC[j]][o] for j in cols) + ' |')
# sum vi-x + kp-x + r-x (NNLO extras), means and errors in quadrature
M += ['', '#### vi-x + kp-x + r-x (sum of the part means, errors in quadrature; total row, pb)', '', '| tau_cut | ' + ' | '.join('X8 | Ex7 | Ex6 | bin5'.split(' | ')) + ' | pull X8-Ex7 | X8-Ex6 | X8-bin5 |', '|---|---|---|---|---|---|---|---|']
J['sum'] = {}
for j in cols:
    v = {s: (sum(J['stats'][p]['%g' % TC[j]][s]['mean'] for p in XP[1:]), math.sqrt(sum(J['stats'][p]['%g' % TC[j]][s]['err']**2 for p in XP[1:]))) for s in SETS}
    pl = [(v['X8'][0] - v[o][0])/math.hypot(v['X8'][1], v[o][1]) for o in ('Ex7', 'Ex6', 'bin5')]; J['sum']['%g' % TC[j]] = {'sets': v, 'pulls': pl}
    M.append('| %.0e | ' % TC[j] + ' | '.join('%.4g +- %.3g' % v[s] for s in SETS) + ' | ' + ' | '.join('%+.1f' % x for x in pl) + ' |')
hdr = ['X8 report numbers: total row (CELL 0.05 0.50) of the extra parts, tau_cut 2e-2..2e-4; plain equal-weight means, err = sd/sqrt(N)'] + info + ['']
open(R + '/X8_total.txt', 'w').write('\n'.join(hdr + L) + '\n'); json.dump(J, open(R + '/X8_total.json', 'w'), indent=1)
open(R + '/X8_tables.md', 'w').write('\n'.join(M) + '\n')
print('\n'.join(info)); print('\n'.join(M))
