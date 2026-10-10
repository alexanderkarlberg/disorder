#!/usr/bin/env python3
"""Our side of the 10 Oct validation sets (Z, W, E), bin5, B commands.
  mk_ours_set.py SET pilot|prod [SCALE]   -> job lists + groups lines (stdout)
Lists in $P/runs/ours5/<SET>[-pilot]/<name>/jobs.list. SCALE multiplies seed counts (prod)."""
import subprocess, sys, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
BIN = P + '/bin5'
setn, mode = sys.argv[1], sys.argv[2]
scale = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
SETS = {  # env common, tables, P2B?, time factor
  'Z': (['EW31="1 0"'], 'zeus5', True, 10),
  'W': (['EW31="2 0"'], 'zeus5w', True, 10),
  'E': (['ZFIX=3'], 'zeus5', False, 6),
}
env0, tab, p2b, tl = SETS[setn]
soft = P + '/tables/soft.tab'
ZP = P + '/tables/' + tab
MODE2 = ['1d-9', '2', '0', '0', 'psmc', '0', '0', '1d-12']
# name: (binary, part, ncall, itmx, seed0, nseeds, extra env, time-limit hours)
parts = {
 'b0': ('sliced21', 'b0', 7500000, 6, 2001, 16, [], 1),
 'b1': ('sliced21', 'b1', 3000000, 6, 2001, 200, [], 2 if tl == 6 else 4),
 'b2': ('sliced21', 'b2', 5000000, 6, 2001, 100, [], 3 if tl == 6 else 5),
 'lo': ('nlo31', 'lo', 3000000, 6, 2001, 200, ['P2B'], 1 if tl == 6 else 2),
 'vi': ('nlo31', 'vi', 1000000, 6, 2001, 250, ['P2B'], tl),
 'kp': ('nlo31', 'kp', 2000000, 6, 2001, 200, ['P2B'], 1 if tl == 6 else 2),
 'r': ('nlo31', 'r', 1000000, 6, 1001, 3500, ['P2B', 'EQ'], tl),
 'rcorr': ('nlo31', 'r', 1000000, 6, 2001, 400, ['P2B', 'EQ', 'TD'], tl),
}
runs = []   # (name, partkey, env, nseeds)
for k in parts:
    runs.append((k, k, [], 1.0))
if setn == 'E':
    for k, f in (('lo', .25), ('vi', .25), ('kp', .25), ('r', .25)):
        runs.append((k + '-x', k, ['P2BEXTRA'], f))
base = P + '/runs/ours5/' + setn + ('-pilot' if mode == 'pilot' else '')
out = []
for name, k, extra, f in runs:
    b, part, nc, it, s0, ns, flags, th = parts[k]
    env = list(env0)
    if 'P2BEXTRA' in extra: env += ['P2BSLICE=1', 'P2BEXTRA=1']
    elif 'P2B' in flags and p2b: env += ['P2BSLICE=1']
    if 'EQ' in flags: env += ['VEGAS_EQUAL=2']
    if 'TD' in flags: env += ['TECHDIFF=1', 'TECHMIN="1d-10 1d-8"']
    if mode == 'pilot':
        s0, n = 9001, 4
    else:
        n = max(1, round(ns * f * scale))
        if name == 'b0': n = ns
    if b == 'sliced21':
        cmd = [BIN + '/sliced21', part, str(nc), str(it), '{seed}', 'zeus', ZP, soft]
    else:
        cmd = [BIN + '/nlo31', part, str(nc), str(it), '{seed}'] + MODE2
    d = '%s/%s' % (base, name)
    os.makedirs(d, exist_ok=True)
    a = [S + '/mkjobs.py', d + '/jobs.list', d, '%d-%d' % (s0, s0 + n - 1), '--pattern', 'RESULT']
    for e in env: a += ['--env', e]
    subprocess.run(a + ['--'] + cmd, check=True, stdout=subprocess.DEVNULL)
    out.append('dis-%s%s-%s %s/jobs.list %d:00:00 1000MB 1 et[07-40]' % (setn, 'p' if mode == 'pilot' else '', name, d, th))
open(base + '/groups.txt', 'w').write('\n'.join(out) + '\n')
print(base + '/groups.txt', len(out), 'groups')
