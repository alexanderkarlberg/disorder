#!/usr/bin/env python3
"""markdown tables for REPORT-2026-10-13 from results/{G6W,G6P,Ex6}_6.json and *_6_parts.txt.  gen_tables6.py SET"""
import json, math, re, sys, os
R = os.path.dirname(os.path.abspath(__file__)) + '/results'; S = sys.argv[1]
J = json.load(open('%s/%s_6.json' % (R, S))); tc = J['tau_cut']
hdr = '| tau_cut | ' + ' | '.join('%.0e' % t for t in tc) + ' |\n|---|' + '---|'*len(tc)
def ve(v, e, sc): return '%.4g +- %.2g' % (v/sc, e/sc)
def parse_parts():
    d = {}; cur = None
    for l in open('%s/%s_6_parts.txt' % (R, S)):
        if l.startswith('=='): cur = l[3:].strip().replace('tauzQ ', ''); d[cur] = {}; last = None
        elif cur and l.strip() and not l.startswith('tau_cut'):
            w = l.split()
            if w[0] == '+-': d[cur][last] = (d[cur][last], [float(x) for x in w[1:]])
            else: last = w[0]; d[cur][last] = [float(x) for x in w[1:]]
    return d
if S != 'Ex6':
    sc = 1e-6; print('Units: 1e-6 pb/GeV^2 (dsigma/dx dQ^2), value +- error of the mean.\n')
    P = parse_parts()
    print('### total bin 0.05-0.5, every part\n'); print(hdr)
    for p in ('b1', 'lo', 'b2', 'vi', 'kp', 'r', 'rcorr', 'd2'):
        v, e = P['0.05-0.5'][p]; print('| %s | ' % p + ' | '.join(ve(a, b, sc) for a, b in zip(v, e)) + ' |')
    b = J['bins']['tauzQ 0.05-0.5']
    for lab in ('NLO', 'NNLO', 'NNLO(no d2)'):
        print('| **%s** | ' % lab + ' | '.join(ve(a, c, sc) for a, c in zip(*b[lab])) + ' |')
    r = b['ref_NLO']; pulls = [(a - r[0])/math.hypot(c, r[1]) for a, c in zip(*b['NLO'])]
    print('| DISENT NLO %.5g +- %.2g: pull | ' % (r[0]/sc, r[1]/sc) + ' | '.join('%+.1f' % x for x in pulls) + ' |')
    print('\n### per bin: NLO (pull vs DISENT) and NNLO\n'); print(hdr)
    for k, v in J['bins'].items():
        nm = k.replace('tauzQ ', '')
        print('| NLO %s | ' % nm + ' | '.join(ve(a, c, sc) for a, c in zip(*v['NLO'])) + ' |')
        r = v['ref_NLO']; print('| pull vs DISENT %.5g +- %.2g | ' % (r[0]/sc, r[1]/sc) + ' | '.join('%+.1f' % ((a - r[0])/math.hypot(c, r[1])) for a, c in zip(*v['NLO'])) + ' |')
        print('| NNLO %s | ' % nm + ' | '.join(ve(a, c, sc) for a, c in zip(*v['NNLO'])) + ' |')
else:
    P = parse_parts(); b = J['bins']; ref = lambda v, k: v['ref_' + k]
    def row(name, v): return '| %s | ' % name + ' | '.join('%.4g +- %.2g' % (a, c) for a, c in zip(*v)) + ' |'
    def pull(v, r): return '| pull vs NNLOJET %.5g +- %.2g | ' % tuple(r) + ' | '.join('%+.1f' % ((a - r[0])/math.hypot(c, r[1])) for a, c in zip(*v)) + ' |'
    t = b['tauzQ 0.05-0.5']
    print('### Ex6 extras per part, total bin 0.05-0.5 (pb); bin5 extras (875/62/50/50 seeds) for comparison\n'); print(hdr)
    for p in ('lo-x', 'vi-x', 'kp-x', 'r-x', 'lo-x5', 'vi-x5', 'kp-x5', 'r-x5'):
        v, e = P['0.05-0.5'][p]; print(row(p.replace('5', ' (bin5)') if p.endswith('5') else p + ' (Ex6)', (v, e)))
    print('\n### total bin, pb\n'); print(hdr)
    for lab in ('NLO', 'NLO_P2B', 'NLO_P2B(bin5 extras)'): print(row(lab, t[lab]))
    print(pull(t['NLO_P2B'], ref(t, 'NLO')))
    for lab in ('NNLO', 'NNLO_P2Bextra', 'NNLO_P2B', 'NNLO_P2B(bin5 extras)'):
        print(row(lab, t[lab]))
        if lab in ('NNLO', 'NNLO_P2B'): print(pull(t[lab], ref(t, 'NNLO')))
    for k, v in b.items():
        if k.endswith('0.05-0.5'): continue
        print('\n#### bin %s\n' % k.replace('tauzQ ', '')); print(hdr)
        for lab in ('NNLO', 'NNLO_P2Bextra', 'NNLO_P2B', 'NNLO_P2B(bin5 extras)'): print(row(lab, v[lab]))
        print(pull(v['NNLO_P2B'], v['ref_NNLO']))
