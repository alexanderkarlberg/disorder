#!/usr/bin/env python3
"""UPDATE 17 Oct: error budget per part for the ZEUS-wide sets (analysis of existing outputs, kp from bin11).
  budget.py  -> results/kp11/BUDGET.txt, BUDGET.json, BUDGET_tables.md
Per set and part: N seeds, mean CPU per job (user+sys from time.log, core-h, as cpu_summary.py), and for the total row (E, G6: tau_zQ 0.05-0.5)
the per-seed sd and 1.48 MAD of Q = N(1e-4), N(1e-5), N(1e-4) - N(1e-5), each part with its NNLO weight (d2 enters only below 1e-4, so its
N(1e-4) is 0 and its difference is -d2(1e-5)); c = sd*sqrt(CPU) and its share of the sum over parts; K = (sum c)^2 [pb^2 core-h], the cost
constant of the optimal allocation (N_p proportional to sd_p/sqrt(CPU_p)): error eps costs K/eps^2 core-h; the core-h for an error of 10% of
NNLOJET's value at 1e-5 (W 0.148, Z 12.9, B 11.1, E 59.5 pb). Per part also the fraction of seeds beyond 20 x 1.48 MAD from the median and
the share of sum (x - mean)^2 they carry. NLO parts (b1, lo; E P2B also lo-x) as a separate group.
B: archived run.log / time.log from $BPROD/prod/B (see combine_b4.sh), kp11 on /ptmp."""
import glob, json, math, os, re, sys
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'; H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results/kp11'
BP = os.environ.get('BPROD', P); TQ = [1e-4, 1e-5]
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
def total(f, cell):
    tc = None
    for l in open(f):
        w = l.split()
        if l.startswith(' tau_cut'): tc = [float(x) for x in w[1:]]
        elif cell == 'total' and l.startswith(' ZCELL total '): return tc, np.array([ff(x) for x in w[4:]])
        elif cell == 'tauzq' and l.startswith(' CELL ') and abs(float(w[1]) - 0.05) < 1e-9 and abs(float(w[2]) - 0.5) < 1e-9: return tc, np.array([ff(x) for x in w[3:]])
    return tc, None
def cpu(s):
    t = open(s + 'time.log').read()
    return (float(re.search(r'User time \(seconds\): ([\d.]+)', t).group(1)) + float(re.search(r'System time \(seconds\): ([\d.]+)', t).group(1)))/3600
def load(dirs, cell):
    v, c, tc = [], [], None
    for d in dirs:
        for s in sorted(glob.glob(d + '/s*/')):
            if not (os.path.exists(s + 'done') and os.path.exists(s + 'run.log')): continue
            t, a = total(s + 'run.log', cell)
            if a is None: continue
            tc = tc or t; assert t == tc, s; v.append(a); c.append(cpu(s))
    return tc, np.array(v), np.array(c)
o5 = P + '/ours5'
def std(d, kp='kp11'): return {'b2': [d + '/b2'], 'vi': [d + '/vi'], 'kp': [d + '/' + kp], 'r': [d + '/r'], 'rcorr': [d + '/rcorr'], 'd2': [d + '/d2'], 'b1': [d + '/b1'], 'lo': [d + '/lo']}
EX7 = {'vi-x': [P + '/ours7/Ex7/vi-x'], 'kp-x': [P + '/ours7/Ex7/kp-x11'], 'r-x': [P + '/ours7/Ex7/r-x', P + '/ours7/Ex7-pilot/r-x'], 'lo-x': [P + '/ours7/Ex7/lo-x']}
SETS = {'W5': (std(o5 + '/W'), 'total', 0.148), 'Z5': (std(o5 + '/Z'), 'total', 12.9),
        'B': ({'b2': [BP + '/prod/B/b2'], 'vi': [BP + '/prod/B/vi'], 'kp': [P + '/prod/B/kp11'], 'r': [BP + '/prod/B/r-4'], 'rcorr': [BP + '/prod/B/rcorr-4'],
               'b1': [BP + '/prod/B/b1'], 'lo': [BP + '/prod/B/lo']}, 'total', 11.1),
        'E5': (std(o5 + '/E'), 'tauzq', 59.5), 'E5+Ex7 (P2B)': ({**std(o5 + '/E'), **EX7}, 'tauzq', 59.5), 'G6P': (std(P + '/ours6/G6P'), 'tauzq', None)}
NNP = ['b2', 'vi', 'kp', 'r', 'rcorr', 'd2', 'vi-x', 'kp-x', 'r-x']; NLP = ['b1', 'lo', 'lo-x']
JA = {}; L = []; md = []
for S, (src, cell, ref) in SETS.items():
    D = {}; tc = None
    for p, dirs in src.items():
        t, v, c = load(dirs, cell); D[p] = (v, c); tc = tc or t
    tca = np.array(tc); j4, j5 = [int(np.argmin(abs(tca - t))) for t in TQ]
    J = {'parts': {}, 'groups': {}}
    for grp, plist in (('NNLO', NNP), ('NLO', NLP)):
        ps = [p for p in plist if p in D]; rowsq = {}
        for p in ps:
            v, c = D[p]; w4 = 0.0 if p == 'd2' else 1.0
            Q = {'N(1e-4)': v[:, j4]*w4, 'N(1e-5)': v[:, j5], 'N(1e-4)-N(1e-5)': v[:, j4]*w4 - v[:, j5]}
            J['parts'][p] = {'N': len(v), 'cpu': float(c.mean()), 'q': {}}
            for q, x in Q.items():
                sd = float(x.std(ddof=1)) if len(x) > 1 else 0.0; med = np.median(x); mad = float(1.4826*np.median(abs(x - med)))
                out = abs(x - med) > 20*mad if mad > 0 else np.zeros(len(x), bool); dev = (x - x.mean())**2
                J['parts'][p]['q'][q] = {'sd': sd, 'mad': mad, 'c': sd*math.sqrt(c.mean()), 'f20': float(out.mean()), 'var20': float(dev[out].sum()/dev.sum()) if dev.sum() > 0 else 0.0}
        for q in ('N(1e-4)', 'N(1e-5)', 'N(1e-4)-N(1e-5)'):
            cs = sum(J['parts'][p]['q'][q]['c'] for p in ps); K = cs**2
            for p in ps: J['parts'][p]['q'][q]['share'] = J['parts'][p]['q'][q]['c']/cs if cs > 0 else 0
            J['groups']['%s %s' % (grp, q)] = {'K': K, 'coreh_10pc': (K/(0.1*ref)**2 if ref else None)}
    JA[S] = J
    # markdown table
    for grp, plist in (('NNLO', NNP), ('NLO', NLP)):
        ps = [p for p in plist if p in D]
        for q in ('N(1e-4)', 'N(1e-5)', 'N(1e-4)-N(1e-5)'):
            g = J['groups']['%s %s' % (grp, q)]; x = J['parts']
            g['coreh_used'] = sum(x[p]['N']*x[p]['cpu'] for p in ps); g['err_now'] = math.sqrt(sum(x[p]['q'][q]['sd']**2/x[p]['N'] for p in ps))
            g['coreh_opt_same_err'] = g['K']/g['err_now']**2 if g['err_now'] > 0 else None
    md += ['### %s (%s row; %s, CPU in core-h per job)' % (S, 'total' if cell == 'total' else 'tau_zQ 0.05-0.5', 'pb/GeV^2' if S.startswith('G6') else 'pb'), '',
           '| part | N | CPU | Q | sd | 1.48 MAD | sd/MAD | sd*sqrt(CPU) | share | f(>20 MAD) | var share of those |', '|---|---|---|---|---|---|---|---|---|---|---|']
    for grp, plist in (('NNLO', NNP), ('NLO', NLP)):
        for p in [p for p in plist if p in D]:
            x = J['parts'][p]
            for q, y in x['q'].items():
                md.append('| %s | %d | %.3f | %s | %.3g | %.3g | %s | %.3g | %.1f%% | %.4f | %.1f%% |' % (p, x['N'], x['cpu'], q, y['sd'], y['mad'], ('%.1f' % (y['sd']/y['mad'])) if y['mad'] > 0 else '-',
                          y['c'], 100*y['share'], y['f20'], 100*y['var20']))
        for q in ('N(1e-4)', 'N(1e-5)', 'N(1e-4)-N(1e-5)'):
            g = J['groups']['%s %s' % (grp, q)]
            md.append('| **%s K** | | | %s | | | | K = %.3g (unit^2 core-h) | | %s | now: %.3g core-h for error %.3g; optimal allocation for the same error %.3g core-h |' % (grp, q, g['K'], ('core-h for 10%% of %.3g pb: %.3g' % (ref, g['coreh_10pc'])) if g['coreh_10pc'] is not None else '(no NNLOJET reference)', g['coreh_used'], g['err_now'], g['coreh_opt_same_err']))
    md.append('')
    print(S, 'done', {p: len(D[p][0]) for p in D})
json.dump(JA, open(R + '/BUDGET.json', 'w'), indent=1)
open(R + '/BUDGET_tables.md', 'w').write('\n'.join(md) + '\n')
