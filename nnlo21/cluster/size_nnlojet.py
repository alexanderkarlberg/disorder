#!/usr/bin/env python3
"""Size the NNLOJET production per part from the production pilot.
For each order, the jobs per part follow the optimal split for the total
(N_i ∝ sigma_i/sqrt(t_i)); then raised until every histogram bin meets its
target (relative to the order's value in that bin, with an absolute floor).
  size_nnlojet.py PRODDIR"""
import glob, os, re, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

prod = sys.argv[1]
ORD = {'LO': r'^LO$', 'NLO': r'^(R|V)$', 'NNLO': r'^(RR[ab]|RV|VV)_\d+$'}
# targets: (total relative, per-bin relative, absolute floor in fb per unit*width)
TGT = {'LO': (1e-4, 1e-4), 'NLO': (3e-3, 1e-2), 'NNLO': (2e-2, 5e-2)}
REF_NNLO_PB = float(os.environ.get('NNLO_REF_PB', 40.0))   # our estimate of the NNLO coefficient (MPP, P2B)


def rd(fn):
    return np.array([[float(v) for v in l.split()] for l in open(fn) if l.strip() and not l.startswith('#')])


parts = {}
for pd in sorted(glob.glob(prod + '/*')):
    p = os.path.basename(pd)
    js = [d for d in glob.glob(pd + '/s*') if os.path.exists(d + '/done')]
    if len(js) < 5:
        continue
    vals, cpu = [], []
    for d in js:
        row = [rd(glob.glob(d + '/DIS.*.cross.s*.dat')[0])[0, 0]]
        for o in ('q2', 'ptavg_12', 'm12'):
            h = rd(glob.glob(d + '/DIS.*.%s.s*.dat' % o)[0])
            row += list(h[:, 3]*(h[:, 2] - h[:, 0]))
        vals.append(row)
        t = open(d + '/time.log').read()
        cpu.append(float(re.search(r'User time \(seconds\): ([\d.]+)', t).group(1)))
    v = np.array(vals)
    parts[p] = (v.mean(0), v.std(0, ddof=1), np.mean(cpu), len(js))
tot = 0
for o, rx in ORD.items():
    ps = [p for p in parts if re.match(rx, p)]
    if not ps:
        continue
    val = sum(parts[p][0] for p in ps)
    if o == 'NNLO':
        val = val.copy(); val[0] = REF_NNLO_PB*1000      # total from our side; bins: |value| with floor
    tr, br = TGT[o]
    w = sum(parts[p][1][0]*math.sqrt(parts[p][2]) for p in ps)
    T = tr*abs(val[0])
    N = {p: parts[p][1][0]*w/math.sqrt(parts[p][2])/T**2 for p in ps}
    # per-bin: scale all parts of the order by a common factor until each bin meets its target
    fl = 0.02*abs(val[0])/6 if o == 'NNLO' else 0       # NNLO bins: floor = 2% of total/6
    f = 1.0
    for b in range(1, len(val)):
        e2 = sum(parts[p][1][b]**2/N[p] for p in ps)
        tb = max(br*abs(val[b]), fl)
        f = max(f, e2/tb**2)
    cost = 0
    print('%s: total %.5g fb, target %.2g; bin factor %.2f' % (o, val[0], T, f))
    for p in ps:
        n = max(int(math.ceil(N[p]*f)), parts[p][3])
        cost += n*parts[p][2]/3600
        print('   %-6s sigma/job %10.4g fb  CPU/job %5.0f s  -> %6d jobs' % (p, parts[p][1][0], parts[p][2], n))
    print('   cost %.0f core-h' % cost)
    tot += cost
print('total %.0f core-h (parts without pilot not included)' % tot)
