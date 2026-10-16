#!/usr/bin/env python3
"""N7 (14 Oct): NNLOJET W and Z after the top-up (results/nnj7/{W,Z}, combine_nnlojet.py over $P/runs/nnj-{W,Z}/prod, all done seeds, equal weights)
against ours bin5 (results/{W,Z}5.json, NNLO incl. d2 below 1e-4, as combine5).
  combine7n.py W|Z -> results/N7_SET_main.txt (per bin: ours, NNLOJET, pulls at every tau_cut, crossing, compatible range), N7_SET.json
Crossing tau_cut: ours - NNLOJET changes sign between two neighbouring tau_cut points (linear interpolation in ln tau_cut); 'none' if no sign change.
Compatible: |pull| < 2 (and < 1) listed as tau_cut ranges."""
import json, math, os, sys
H = os.path.dirname(os.path.abspath(__file__)); R = H + '/results'; S = sys.argv[1]
R0 = R  # KREDO: KP11=1 reads ours from results/kp11/{W,Z}5.json (bin11 kp) and writes to results/kp11; NNLOJET stays in results/
if os.environ.get('KP11') == '1': R = R + '/kp11'
J = json.load(open('%s/%s5.json' % (R, S))); tc = J['tau_cut']
def ref(d, o):
    out = {}
    w = [float(v) for v in open('%s/%s.cross.dat' % (d, o)).read().split('\n')[1].split()]; out['total 0-0'] = (w[0]/1000, w[1]/1000)
    for ob in ('q2', 'ptavg_12', 'm12'):
        for l in open('%s/%s.%s.dat' % (d, o, ob)):
            if l.startswith('#') or not l.strip(): continue
            w = [float(v) for v in l.split()]; lo, hi = w[0], w[2]
            out['%s %g-%g' % ({'ptavg_12': 'ptavg'}.get(ob, ob), lo, hi)] = (w[3]*(hi - lo)/1000, w[4]*(hi - lo)/1000)
    return out
new = {o: ref('%s/nnj7/%s' % (R0, S), o) for o in ('LO', 'NLO', 'NNLO')}; old = ref('%s/nnj5/%s' % (R0, S), 'NNLO')
def rng(flags):
    out, cur = [], None
    for t, f in zip(tc, flags):
        if f and cur is None: cur = [t, t]
        elif f: cur[1] = t
        elif cur: out.append(cur); cur = None
    if cur: out.append(cur)
    return ', '.join('%.0e..%.0e' % (a, b) if a != b else '%.0e' % a for a, b in out) or 'none'
L = ['NNLOJET %s (top-up, results/nnj7/%s) vs ours bin5 (NNLO incl. d2 below 1e-4); pb; pull = (ours - NNLOJET)/sqrt(e_ours^2 + e_NNLOJET^2)' % (S, S),
     'tau_cut            ' + ''.join('%11.0e' % t for t in tc)]
out = {'tau_cut': tc, 'bins': {}}
for k, v in J['bins'].items():
    m, e = v['NNLO']; r, er = new['NNLO'][k]; ro, ero = old[k]
    p = [(a - r)/math.hypot(b, er) for a, b in zip(m, e)]
    cr = []
    for j in range(len(tc) - 1):
        d1, d2 = m[j] - r, m[j + 1] - r
        if d1*d2 < 0:
            f = d1/(d1 - d2); cr.append(math.exp(math.log(tc[j]) + f*(math.log(tc[j + 1]) - math.log(tc[j]))))
    L += ['', '== %s   NNLOJET NNLO %.5g +- %.3g (before top-up %.5g +- %.3g); NLO %.5g +- %.2g; LO %.6g' % (k, r, er, ro, ero, new['NLO'][k][0], new['NLO'][k][1], new['LO'][k][0]),
          'ours NNLO          ' + ''.join('%11.5g' % x for x in m), '  +-               ' + ''.join('%11.2g' % x for x in e),
          'pull               ' + ''.join('%+11.1f' % x for x in p),
          'crossing tau_cut: ' + (', '.join('%.2g' % c for c in cr) or 'none') + ';  |pull|<2: ' + rng([abs(x) < 2 for x in p]) + ';  |pull|<1: ' + rng([abs(x) < 1 for x in p])]
    out['bins'][k] = dict(ours=[m, e], nnlojet=[r, er], nnlojet_old=[ro, ero], pull=p, crossing=cr)
open('%s/N7_%s_main.txt' % (R, S), 'w').write('\n'.join(L) + '\n'); json.dump(out, open('%s/N7_%s.json' % (R, S), 'w'), indent=1)
print('\n'.join(L[:12]))
# summary over bins: worst-bin pull, mean and rms of the pulls per tau_cut; per bin the new and old NNLOJET value and the |pull|<2 tau_cut set
B = out['bins']; Ls = ['%s: per tau_cut, over the %d bins: worst-bin pull, mean pull, rms pull (bins are correlated through the totals)' % (S, len(B))]
for j, t in enumerate(tc):
    w = max(B, key=lambda k: abs(B[k]['pull'][j])); ps = [B[k]['pull'][j] for k in B]
    Ls.append('  %.0e worst %-14s %+.1f  mean %+.2f rms %.2f' % (t, w, B[w]['pull'][j], sum(ps)/len(ps), math.sqrt(sum(p*p for p in ps)/len(ps))))
Ls.append('per bin: NNLOJET NNLO new (old) pb; crossing tau_cut; tau_cut with |pull|<2')
for k, b in B.items():
    Ls.append('  %-14s %9.4g +- %-7.2g (%9.4g +- %-7.2g) cross %-26s |pull|<2: %s' % (k, *b['nnlojet'], *b['nnlojet_old'], ','.join('%.2g' % c for c in b['crossing']) or 'none', ' '.join('%.0e' % t for t, p in zip(tc, b['pull']) if abs(p) < 2)))
open('%s/N7_%s_summary.txt' % (R, S), 'w').write('\n'.join(Ls) + '\n')
