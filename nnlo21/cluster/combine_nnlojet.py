#!/usr/bin/env python3
"""Combine NNLOJET production jobs of the ZEUS dijet runs.

  combine_nnlojet.py PRODDIR OUTDIR

PRODDIR/<part>/s<seed>/ hold one job each (part = LO, R, V, RRa_1, ..., VV_7).
Per part, the seeds are averaged with EQUAL weights (value = mean, error =
seed scatter/sqrt(N)); NNLOJET's own per-job errors are only used for the
diagnostic ratio scatter/quoted. Parts are summed per coefficient:
LO = LO, NLO = R + V, NNLO = RR* + RV* + VV*. OUTDIR gets <order>.cross.dat
(sigma in fb) and <order>.<obs>.dat (dsigma/dX in fb per unit, columns lo
center hi value error), the format of combine_zeus.py --nnlojet, plus
parts.txt (per part: seeds, total, errors, CPU)."""
import glob, math, os, re, sys
import numpy as np

ORDERS = {'LO': ('LO',), 'NLO': ('R', 'V'), 'NNLO': ('RR', 'RV', 'VV')}
OBS = ('cross', 'q2', 'ptavg_12', 'm12')


def read(fn):
    rows = []
    for l in open(fn):
        if l.startswith('#') or not l.strip():
            continue
        rows.append([float(v) for v in l.split()])
    return np.array(rows)


def cpu(d):
    try:
        t = open(os.path.join(d, 'time.log')).read()
        u = float(re.search(r'User time \(seconds\): ([\d.]+)', t).group(1))
        s = float(re.search(r'System time \(seconds\): ([\d.]+)', t).group(1))
        return u + s
    except (OSError, AttributeError):
        return float('nan')


def main():
    prod, out = sys.argv[1:3]
    os.makedirs(out, exist_ok=True)
    comb = {}
    rep = []
    for pd in sorted(glob.glob(os.path.join(prod, '*'))):
        part = os.path.basename(pd)
        if not os.path.isdir(pd):
            continue
        jobs = [d for d in sorted(glob.glob(os.path.join(pd, 's*'))) if os.path.exists(os.path.join(d, 'done'))]
        if not jobs:
            continue
        data = {o: [] for o in OBS}
        cpus = []
        for d in jobs:
            fs = {o: glob.glob(os.path.join(d, 'DIS.*.%s.s*.dat' % o)) for o in OBS}
            if not all(len(fs[o]) == 1 for o in OBS):
                print('skipping', d, {o: len(fs[o]) for o in OBS}); continue
            for o in OBS:
                data[o].append(read(fs[o][0]))
            cpus.append(cpu(d))
        n = len(data['cross'])
        comb[part] = {}
        for o in OBS:
            a = np.array(data[o])            # seeds x rows x cols
            col = 0 if o == 'cross' else 3
            v = a[:, :, col]
            m = v.mean(axis=0)
            e = v.std(axis=0, ddof=1)/math.sqrt(n) if n > 1 else a[0, :, col + 1]
            q = np.sqrt((a[:, :, col + 1]**2).sum(axis=0))/n
            comb[part][o] = (a[0], m, e, q)
        m, e, q = comb[part]['cross'][1][0], comb[part]['cross'][2][0], comb[part]['cross'][3][0]
        c = np.nanmean(cpus)
        rep.append('%-7s %5d seeds  sigma %14.6g +- %10.3g fb (scatter)  quoted %10.3g  ratio %5.2f  '
                   'CPU/job %8.0f s  CPU*err^2/job-equiv %10.4g' % (part, n, m, e, q, e/q if q else 0, c,
                                                                      c*n*e**2))
    for order, kinds in ORDERS.items():
        parts = [p for p in comb if re.match(r'^(%s)[ab]?(_\d+)?$' % '|'.join(kinds), p)]
        if not parts:
            continue
        for o in OBS:
            base = comb[parts[0]][o][0]
            m = sum(comb[p][o][1] for p in parts)
            e = np.sqrt(sum(comb[p][o][2]**2 for p in parts))
            with open(os.path.join(out, '%s.%s.dat' % (order, o)), 'w') as f:
                f.write('# %s %s: equal-weight seed average per part, sum of parts %s\n' % (order, o, ' '.join(parts)))
                for i in range(len(m)):
                    if o == 'cross':
                        f.write('%.10e %.10e\n' % (m[i], e[i]))
                    else:
                        f.write('%.10e %.10e %.10e %.10e %.10e\n' % (base[i, 0], base[i, 1], base[i, 2], m[i], e[i]))
        tot = sum(comb[p]['cross'][1][0] for p in parts)
        err = math.sqrt(sum(comb[p]['cross'][2][0]**2 for p in parts))
        rep.append('%-4s coefficient: %.6g +- %.3g pb (%s)' % (order, tot/1000, err/1000, ' '.join(parts)))
    open(os.path.join(out, 'parts.txt'), 'w').write('\n'.join(rep) + '\n')
    print('\n'.join(rep))


main()
