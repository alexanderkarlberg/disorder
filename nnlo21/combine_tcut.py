#!/usr/bin/env python3
"""Combine the tau_cut test of tau_2-sliced NNLO DIS 2+1 at fixed (x, Q^2).

Inputs: output files of `sliced21 b1|b2` (below the cut) and of `nlo31 lo|vi|kp|r`
in mode 1 (above the cut), each with CELL lines (tau_zQ bin, then one value per
tau_cut) in dsigma/dx dQ^2 [pb/GeV^2]. Per part, the seeds are combined with
equal weights, error = seed scatter/sqrt(N) (the r part has heavy tails;
inverse-variance weights of VEGAS errors are biased, see docs/notebook.md 3 Oct).

Prints per tau_zQ bin and tau_cut:
  NLO:  b1 + lo                (O(alpha_s^2) coefficient of 2+1; reference
                                DISENT's NLO coefficient x 1/x, option --disent)
  NNLO: b2 + vi + kp + r       (O(alpha_s^3) coefficient), which must not
                                depend on tau_cut for small tau_cut.
--json FILE writes everything for plots.
Usage: combine_tcut.py [--disent tau2_combine_output] [--x 0.01] [--json f] files..."""
import sys, re, json, math
import numpy as np


def read(fn):
    part, tc, cells = None, None, {}
    for l in open(fn):
        m = re.match(r'\s*(sliced21|nlo31) part (\S+)', l)
        if m:
            part = m.group(2)
        if l.startswith(' tau_cut'):
            tc = [float(v) for v in l.split()[1:]]
        if l.startswith(' CELL '):
            w = l.split()
            cells[(float(w[1]), float(w[2]))] = np.array([float(v) for v in w[3:]])
    if part is None or not cells:
        return None
    return part, tc, cells


def main():
    args = sys.argv[1:]
    dis, x, js = None, 0.01, None
    while args and args[0].startswith('--'):
        o = args.pop(0)
        if o == '--disent': dis = args.pop(0)
        elif o == '--x': x = float(args.pop(0))
        elif o == '--json': js = args.pop(0)
    runs = {}
    tc = None
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
    ref = {}
    if dis:
        txt = open(dis).read()
        for blk in re.split(r'\ntau_zQ in ', txt)[1:]:
            lo, hi = map(float, re.match(r'\[([\d.]+),([\d.]+)\)', blk).groups())
            m = re.search(r'NLO coefficient \(DISENT\) ([-+.\deE]+) \+- ([-+.\deE]+)', blk)
            ref[(lo, hi)] = (float(m.group(1))/x, float(m.group(2))/x)
    out = {'tau_cut': tc, 'bins': {}}
    for b in sorted(next(iter(comb.values()))):
        def tot(parts):
            if not all(p in comb for p in parts):
                return None
            m = sum(comb[p][b][0] for p in parts)
            e = np.sqrt(sum(comb[p][b][1]**2 for p in parts))
            return m, e
        nlo = tot(['b1', 'lo']); nnlo = tot(['b2', 'vi', 'kp', 'r'])
        print('\ntau_zQ in [%.2f, %.2f)' % b + ('   DISENT NLO x 1/x: %.5g +- %.2g' % ref[b] if b in ref else ''))
        print('   tau_cut          NLO (b1+lo)            NNLO (b2+vi+kp+r)       below b2       above vi+kp+r')
        for i, t in enumerate(tc):
            s = '  %8.1e' % t
            s += ('  %12.5g +- %9.2g' % (nlo[0][i], nlo[1][i])) if nlo else ' ' * 28
            s += ('  %12.5g +- %9.2g' % (nnlo[0][i], nnlo[1][i])) if nnlo else ' ' * 28
            if nnlo:
                ab = sum(comb[p][b][0][i] for p in ('vi', 'kp', 'r'))
                s += '  %12.5g  %12.5g' % (comb['b2'][b][0][i], ab)
            print(s)
        out['bins']['%.2f-%.2f' % b] = {
            'nlo': [list(nlo[0]), list(nlo[1])] if nlo else None,
            'nnlo': [list(nnlo[0]), list(nnlo[1])] if nnlo else None,
            'disent_nlo': ref.get(b),
            'parts': {p: [list(comb[p][b][0]), list(comb[p][b][1])] for p in comb}}
    if js:
        json.dump(out, open(js, 'w'), indent=1)


if __name__ == '__main__':
    main()
