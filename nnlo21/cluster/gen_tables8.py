#!/usr/bin/env python3
"""markdown tables for REPORT-2026-10-15-J8 from results/{J8_PH9,J7_PH}.json.  gen_tables8.py
Per bin: J8 NLO and NNLO (1e-9 pb/GeV2); slopes J8 next to J7 PH (% of the bin's own |NLO(1e-5)|); for the total and the ZCUTS rows also the
absolute slope (1e-9 pb/GeV2 per unit ln tau_cut) and the slope in % of the ZEUS total NLO(1e-5) (one common denominator for all variants)."""
import json, os
R = os.path.dirname(os.path.abspath(__file__)) + '/results'
J = json.load(open(R + '/J8_PH9.json')); P = json.load(open(R + '/J7_PH.json')); tc = J['tau_cut']
W = ['5e-3->5e-4', '2e-3->2e-4', '1e-3->1e-4', '1e-4->1e-5']
def hd(tc): return '| tau_cut | ' + ' | '.join('%.0e' % t for t in tc) + ' |\n|---|' + '---|'*len(tc)
def fs(x): return '%+.1f +- %.1f (%+.1f)' % tuple(x)
print('### J8 per bin (dsigma/dx dQ2, 1e-9 pb/GeV2)')
for k in J['bins']:
    print('\n#### %s\n' % k); print(hd(tc))
    for lab in ('NLO', 'NNLO', 'NNLO(no d2)'):
        m, e = J['bins'][k][lab]; print('| %s | ' % lab + ' | '.join('%.4g +- %.2g' % (a*1e9, c*1e9) for a, c in zip(m, e)) + ' |')
for lab in ('NNLO', 'NNLO(no d2)'):
    print('\n### slopes, %s, %% of the bin |NLO(1e-5)| per unit ln tau_cut; (pull); J7 PH where the bin exists there\n' % lab)
    print('| bin | set | ' + ' | '.join(W) + ' |\n|---|---|' + '---|'*len(W))
    for k in J['bins']:
        print('| %s | J8 | ' % k + ' | '.join(fs(x) for x in J['slopes']['%s %s' % (lab, k)]) + ' |')
        kk = '%s %s' % (lab, k)
        if kk in P['slopes']: print('| | J7 PH | ' + ' | '.join(fs(x) for x in P['slopes'][kk]) + ' |')
    if lab == 'NNLO':
        kk = '%s q2 250-500' % lab; print('| q2 250-500 | J7 PH | ' + ' | '.join(fs(x) for x in P['slopes'][kk]) + ' |')
den0 = J['nlo_den']['total 0-0']
for lab in ('NNLO', 'NNLO(no d2)'):
    print('\n### ZEUS rows vs ZCUTS rows, %s: slope in 1e-9 pb/GeV2 per unit ln tau_cut [and in %% of the ZEUS total |NLO(1e-5)| = %.4g e-9]; (pull)\n' % (lab, den0*1e9))
    print('| bin | NLO(1e-5) [1e-9] | ' + ' | '.join(W) + ' |\n|---|---|' + '---|'*len(W))
    for k in ('total 0-0', 'et4 0-0', 'noeta 0-0', 'et4ne 0-0', 'ptavg 8-15', 'et4 8-15', 'noeta 8-15', 'et4ne 8-15'):
        d = J['nlo_den'][k]; c = []
        for s, e, p in J['slopes']['%s %s' % (lab, k)]:
            a = s*d/100; ea = e*d/100; c.append('%+.2f +- %.2f [%+.1f%%] (%+.1f)' % (a*1e9, ea*1e9, 100*a/den0, p))
        print('| %s | %.4g | ' % (k, d*1e9) + ' | '.join(c) + ' |')
