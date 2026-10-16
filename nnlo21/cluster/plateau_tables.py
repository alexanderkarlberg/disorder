#!/usr/bin/env python3
"""markdown tables for the UPDATE 16 Oct (c) addendum from results/kp11/PLATEAU.json -> results/kp11/PLATEAU_tables.md"""
import json, os
H = os.path.dirname(os.path.abspath(__file__)); A = json.load(open(H + '/results/kp11/PLATEAU.json'))
ORDER = ['J8', 'CH8q', 'CH8g', 'Q8f', 'J7PH', 'J7W', 'J7vV9', 'J7vZ', 'Z', 'W', 'B', 'E', 'E_P2B_bin5', 'E_P2B_Ex7', 'G6P', 'G6W']
def tot(S): b = A[S]['bins']; return 'total 0-0' if 'total 0-0' in b else 'tauzQ 0.05-0.5'
def sc(S): return 1e9 if S.startswith(('J8', 'CH8', 'Q8', 'J7')) else (1e6 if S.startswith('G6') else 1)
UN = {1e9: '1e-9 pb/GeV^2', 1e6: '1e-6 pb/GeV^2', 1: 'pb'}
o = ['### 1. NNLO slopes per unit ln tau_cut, total row (J8/CH8/Q8f/J7 in 1e-9 pb/GeV^2, G6 in 1e-6 pb/GeV^2, Z/W/B/E in pb): value +- e1 [e2] (value/e2)', '']
W = [w for w, *_ in A['J8']['bins']['total 0-0']['NNLO']['slopes']]
o += ['| set | ' + ' | '.join(W) + ' |', '|---|' + '---|'*len(W)]
for S in ORDER:
    f = sc(S); o.append('| %s | ' % S + ' | '.join('%+.3g +- %.2g [%.2g] (%+.1f)' % (v*f, e1*f, e2*f, v/e2 if e2 else 0) for w, v, e1, e2 in A[S]['bins'][tot(S)]['NNLO']['slopes']) + ' |')
o += ['', '### 2. Plateau P (mean over tau_cut 1e-4, 3e-5, 1e-5), total row; pull of P and pulls at 1e-4 / 3e-5 / 1e-5 against the reference', '',
      '| set | order | P | N(1e-4) | N(3e-5) | N(1e-5) | reference | pull P | pulls 1e-4 / 3e-5 / 1e-5 |', '|---|---|---|---|---|---|---|---|---|']
tc = A['J8']['tau_cut']; jt = [tc.index(t) for t in (1e-4, 3e-5, 1e-5)]
for S in ORDER:
    f = sc(S); b = A[S]['bins'][tot(S)]
    for lab in ('NLO', 'NNLO'):
        x = b[lab]; refs = x['refs'] or {'': None}
        for rn, r in refs.items():
            o.append('| %s | %s | %.4g +- %.2g | %s | %s |' % (S, lab, x['P'][0]*f, x['P'][1]*f, ' | '.join('%.4g +- %.2g' % (x['mean'][j]*f, x['err'][j]*f) for j in jt),
                     '%s %.4g +- %.2g | %+.2f | %s' % (rn, r['ref'][0]*f, r['ref'][1]*f, r['pull_P'], ' / '.join('%+.1f' % r['pulls'][j] for j in jt)) if r else '- | - | -'))
o += ['', '### 3. Per bin: NNLO plateau P with its pull, and the pulls at 1e-4 / 3e-5 / 1e-5 (Z, W against bin5 NNLOJET and N7; B, E against bin5 NNLOJET)', '']
for S in ('Z', 'W', 'B', 'E', 'E_P2B_bin5', 'E_P2B_Ex7'):
    o += ['#### %s (pb)' % S, '', '| bin | P (NNLO) | ref | pull P | pulls 1e-4 / 3e-5 / 1e-5 | P (NLO) | NLO ref | pull P (NLO) |', '|---|---|---|---|---|---|---|---|']
    for n, b in A[S]['bins'].items():
        for rn in (b['NNLO']['refs'] or {'': None}):
            r = b['NNLO']['refs'].get(rn); rl = b['NLO']['refs'].get(rn)
            o.append('| %s | %.4g +- %.2g | %s | %s | %s | %.4g +- %.2g | %s | %s |' % (n, b['NNLO']['P'][0], b['NNLO']['P'][1], ('%s %.4g +- %.2g' % (rn, *r['ref'])) if r else '-',
                     ('%+.2f' % r['pull_P']) if r else '-', ' / '.join('%+.1f' % r['pulls'][j] for j in jt) if r else '-', b['NLO']['P'][0], b['NLO']['P'][1],
                     ('%.4g +- %.2g' % tuple(rl['ref'])) if rl else '-', ('%+.2f' % rl['pull_P']) if rl else '-'))
    o.append('')
open(H + '/results/kp11/PLATEAU_tables.md', 'w').write('\n'.join(o) + '\n'); print('\n'.join(o[:40]))
