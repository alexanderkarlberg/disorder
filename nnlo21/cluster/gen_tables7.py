#!/usr/bin/env python3
"""markdown tables for REPORT-2026-10-14 from results/{N7_W,N7_Z,J7_PH,J7_W,Ex7_7}.json.  gen_tables7.py N7 W|Z | J7 | Ex7"""
import json, math, os, sys
R = os.path.dirname(os.path.abspath(__file__)) + '/results'; K = sys.argv[1]
def hd(tc): return '| tau_cut | ' + ' | '.join('%.0e' % t for t in tc) + ' |\n|---|' + '---|'*len(tc)
if K == 'N7':
    S = sys.argv[2]; J = json.load(open('%s/N7_%s.json' % (R, S))); tc = J['tau_cut']
    for k in ('total 0-0', 'ptavg 30-60') + tuple(x for x in J['bins'] if x not in ('total 0-0', 'ptavg 30-60')):
        b = J['bins'][k]; m, e = b['ours']; r, er = b['nnlojet']; ro, ero = b['nnlojet_old']
        print('\n#### %s %s: NNLOJET NNLO %.5g +- %.2g pb (before top-up %.5g +- %.2g)\n' % (S, k, r, er, ro, ero)); print(hd(tc))
        print('| ours NNLO | ' + ' | '.join('%.4g +- %.2g' % (a, c) for a, c in zip(m, e)) + ' |')
        print('| pull | ' + ' | '.join('%+.1f' % p for p in b['pull']) + ' |')
        print('\ncrossing tau_cut: %s' % (', '.join('%.2g' % c for c in b['crossing']) or 'none'))
if K == 'J7':
    P = json.load(open(R + '/J7_PH.json')); W = json.load(open(R + '/J7_W.json')); tc = P['tau_cut']
    for k in P['bins']:
        print('\n#### bin %s (dsigma/dx dQ2, 1e-6 pb/GeV2)\n' % k); print(hd(tc))
        for nm, J in (('PH', P), ('W', W)):
            for lab in ('NLO', 'NNLO'):
                m, e = J['bins'][k][lab]; print('| %s %s | ' % (nm, lab) + ' | '.join('%.4g +- %.2g' % (a*1e6, c*1e6) for a, c in zip(m, e)) + ' |')
    print('\n### slopes, % of |NLO(1e-5)| per unit ln tau_cut, windows 5e-3->5e-4 | 2e-3->2e-4 | 1e-3->1e-4 | 1e-4->1e-5; (pull); last column ratio W/PH\n')
    for lab in ('NNLO', 'NNLO(no d2)'):
        print('\n**%s**\n' % lab); print('| bin | boson | 5e-3->5e-4 | 2e-3->2e-4 | 1e-3->1e-4 | 1e-4->1e-5 |\n|---|---|---|---|---|---|')
        for k in P['bins']:
            sp = P['slopes']['%s %s' % (lab, k)]; sw = W['slopes']['%s %s' % (lab, k)]
            print('| %s | PH | ' % k + ' | '.join('%+.1f +- %.1f (%+.1f)' % tuple(x) for x in sp) + ' |')
            print('| | W | ' + ' | '.join('%+.1f +- %.1f (%+.1f)' % tuple(x) for x in sw) + ' |')
            rr = []
            for a, b in zip(sw, sp):
                if abs(b[0]) > 0 and abs(a[0]) > 0:
                    r = a[0]/b[0]; er = abs(r)*math.hypot(a[1]/a[0], b[1]/b[0]); rr.append('%.2f +- %.2f' % (r, er))
                else: rr.append('-')
            print('| | W/PH | ' + ' | '.join(rr) + ' |')
if K == 'Ex7':
    J = json.load(open(R + '/Ex7_7.json')); tc = J['tau_cut']; B = J['bins']
    VS = ('Ex7', 'Ex6+Ex7', 'Ex7+b5lvk', 'Ex6', 'bin5')
    def row(nm, v): return '| %s | ' % nm + ' | '.join('%.4g +- %.2g' % (a, c) for a, c in zip(*v)) + ' |'
    def pull(v, r): return '| pull vs NNLOJET %.5g +- %.2g | ' % tuple(r) + ' | '.join('%+.1f' % ((a - r[0])/math.hypot(c, r[1])) for a, c in zip(*v)) + ' |'
    for k in B:
        if k.endswith('0.05-0.5') or True:
            t = B[k]; print('\n#### bin %s\n' % k.replace('tauzQ ', '')); print(hd(tc))
            print(row('NLO (plain)', t['NLO']))
            for v in VS[:3]: print(row('NLO_P2B [%s]' % v, t['NLO_P2B [%s]' % v]))
            print(pull(t['NLO_P2B [Ex7]'], t['ref_NLO']))
            print(row('NNLO (plain)', t['NNLO']))
            for v in VS: print(row('NNLO_P2B [%s]' % v, t['NNLO_P2B [%s]' % v]))
            for v in VS[:3]: print(pull(t['NNLO_P2B [%s]' % v], t['ref_NNLO']) .replace('pull vs', 'pull [%s] vs' % v))
