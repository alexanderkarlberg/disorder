#!/usr/bin/env python3
"""Combine the ZEUS-like dijet runs of tau_2-sliced NNLO DIS 2+1 (mode 2) and
compare with NNLOJET (nnlo21/validation/nnlojet_epLJJ_zeus2j.run).

Inputs: outputs of `sliced21 b0|b1|b2 ... zeus <tables>` (below the cut) and
`nlo31 lo|vi|kp|r ... 2 0 0 psmc` (above the cut), with ZCELL lines
(observable bin, then one value per tau_cut), sigma per bin in pb. Seeds of a
part are combined with equal weights (error = seed scatter/sqrt(N)), as in
combine_tcut.py.

Prints per observable bin and tau_cut:
  LO   = b0                         (tau_cut independent)
  NLO  = b1 + lo                    (O(alpha_s^2) coefficient)
  NNLO = b2 + vi + kp + r           (O(alpha_s^3) coefficient)
and, with --nnlojet DIR, NNLOJET's numbers per order from the combined
NNLOJET result files in DIR named <order>.<observable>.dat (order = LO, NLO,
NNLO for the coefficients; NNLOJET's histogram files hold dsigma/dX in fb per
unit, the cross file sigma in fb; both converted to pb per bin).
--json FILE writes everything for plots.
Usage: combine_zeus.py [--nnlojet DIR] [--json f] files..."""
import sys, re, json, math, os
import numpy as np


def read(fn):
    part, tc, cells = None, None, {}
    for l in open(fn):
        m = re.match(r'\s*(sliced21|nlo31) part (\S+)', l)
        if m:
            part = m.group(2)
        if l.startswith(' tau_cut'):
            tc = [float(v) for v in l.split()[1:]]
        if l.startswith(' ZCELL '):
            w = l.split()
            cells[(w[1], float(w[2]), float(w[3]))] = np.array([float(v) for v in w[4:]])
    if part is None or not cells:
        return None
    return part, tc, cells


def nnlojet(d):
    """{(obs, lo, hi): (value, error)} in pb per bin, per order"""
    out = {}
    if not d:
        return out
    for order in ('LO', 'NLO', 'NNLO'):
        ref = {}
        fn = os.path.join(d, '%s.cross.dat' % order)
        if os.path.exists(fn):
            for l in open(fn):
                if l.startswith('#') or not l.strip():
                    continue
                w = [float(v) for v in l.split()]
                ref[('total', 0.0, 0.0)] = (w[0]/1000, w[1]/1000)
        for obs, name in (('q2', 'q2'), ('ptavg', 'ptavg_12'), ('m12', 'm12')):
            fn = os.path.join(d, '%s.%s.dat' % (order, name))
            if not os.path.exists(fn):
                continue
            for l in open(fn):
                if l.startswith('#') or not l.strip():
                    continue
                w = [float(v) for v in l.split()]
                lo, hi = w[0], w[2]
                ref[(obs, lo, hi)] = (w[3]*(hi - lo)/1000, w[4]*(hi - lo)/1000)
        if ref:
            out[order] = ref
    return out


def main():
    args = sys.argv[1:]
    nd, js = None, None
    while args and args[0].startswith('--'):
        o = args.pop(0)
        if o == '--nnlojet': nd = args.pop(0)
        elif o == '--json': js = args.pop(0)
    runs, tc = {}, None
    for fn in args:
        r = read(fn)
        if r is None:
            continue
        part, tc, cells = r
        runs.setdefault(part, []).append(cells)
    comb = {}
    for part, lst in runs.items():
        n = len(lst)
        comb[part] = {}
        for b in lst[0]:
            a = np.array([c[b] for c in lst])
            m = a.mean(axis=0)
            e = a.std(axis=0, ddof=1)/math.sqrt(n) if n > 1 else np.zeros_like(m)
            comb[part][b] = (m, e)
        print('%-3s %3d seeds' % (part, n))
    ref = nnlojet(nd)
    out = {'tau_cut': tc, 'bins': {}}
    for b in next(iter(comb.values())):
        def tot(parts):
            if not all(p in comb for p in parts):
                return None
            return (sum(comb[p][b][0] for p in parts), np.sqrt(sum(comb[p][b][1]**2 for p in parts)))
        lo, nlo, nnlo = tot(['b0']), tot(['b1', 'lo']), tot(['b2', 'vi', 'kp', 'r'])
        key = '%s %g-%g' % b
        print('\n%s' % key)
        for order, v in (('LO', lo), ('NLO', nlo), ('NNLO', nnlo)):
            if v is None:
                continue
            r = ref.get(order, {}).get(b)
            s = '  %-4s' % order + ' '.join('%10.5g' % x for x in v[0])
            print(s)
            print('   +- ' + ' '.join('%10.2g' % x for x in v[1]))
            if r:
                print('   NNLOJET %10.5g +- %.2g   (ours at smallest tau_cut %.5g +- %.2g: %+.1f sigma)'
                      % (r[0], r[1], v[0][-1], v[1][-1], (v[0][-1] - r[0])/math.hypot(v[1][-1], r[1])))
        out['bins'][key] = {o: ([list(v[0]), list(v[1])] if v else None) for o, v in
                            (('LO', lo), ('NLO', nlo), ('NNLO', nnlo))}
        out['bins'][key]['nnlojet'] = {o: ref[o].get(b) for o in ref}
        out['bins'][key]['parts'] = {p: [list(comb[p][b][0]), list(comb[p][b][1])] for p in comb}
    if js:
        json.dump(out, open(js, 'w'), indent=1)


if __name__ == '__main__':
    main()
