import math, numpy as np
from g8_pilot import *
L10 = math.log(10); T = 18e-9
tc, nm, _, _ = load('PH9', 'b2'); tc = np.array(tc); i = nm.index('total 0-0'); ja = list(tc).index(1e-3) if 1e-3 in tc else int(np.argmin(abs(tc-1e-3))); jb = int(np.argmin(abs(tc-1e-4)))
def var(A):   # per-seed variances of slope: signed (combine) and direct
    a, b = A[:, i, ja], A[:, i, jb]
    return (a.var(ddof=1) - b.var(ddof=1)) / L10**2, (a - b).var(ddof=1) / L10**2, a.std(ddof=1), b.std(ddof=1)
data = {}
for setn, part in (('PH9', 'b2'), ('PH9', 'rcorr'), ('PH9', 'vi'), ('PH9', 'kp'), ('PH9', 'r'), ('G8', 'vi'), ('G8', 'kp'), ('G8', 'r')):
    A = load(setn, part)[2]; v = var(A); data[(setn, part)] = v
    print('%-4s %-6s N=%4d  sd_a %.3g sd_b %.3g  slope var/seed signed %.3g  direct %.3g  -> slope err at N=1: %.3g (direct %.3g)' % (setn, part, len(A), v[2], v[3], v[0], v[1], math.sqrt(abs(v[0])), math.sqrt(v[1])))
cost = {'vi': cpu('G8', 'vi').mean(), 'kp': cpu('G8', 'kp').mean(), 'r': cpu('G8', 'r').mean()}
print('cost', cost)
fixed = sum(data[('PH9', p)][1] / n for p, n in (('b2', 200), ('rcorr', 200)))
print('J8 b2+rcorr slope err (direct, N=200):', math.sqrt(fixed), ' with J8 psmc vi,kp,r (240,200,1600):', math.sqrt(fixed + sum(data[('PH9', p)][1] / n for p, n in (('vi', 240), ('kp', 200), ('r', 1600)))))
for k in (1, 0):
    f = sum(abs(data[('PH9', p)][k]) / n for p, n in (('b2', 200), ('rcorr', 200)))
    rem = T**2 - f
    if rem <= 0: print('fixed part exceeds target'); continue
    ps = ('vi', 'kp', 'r'); v = {p: abs(data[('G8', p)][k]) for p in ps}
    s = sum(math.sqrt(v[p] * cost[p]) for p in ps)
    Nn = {p: math.sqrt(v[p] / cost[p]) * s / rem for p in ps}
    print('mode', 'direct' if k else 'signed', 'fixed err %.3g' % math.sqrt(f), 'N:', {p: round(Nn[p]) for p in ps}, 'core-h:', {p: round(Nn[p] * cost[p]) for p in ps}, 'total %.0f' % sum(Nn[p] * cost[p] for p in ps))
