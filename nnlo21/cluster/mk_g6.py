#!/usr/bin/env python3
"""12 Oct round, bin6: G6 (fixed point x=0.1 Q2=5000, mode 1, plain slicing; photon G6P, W- G6W) and Ex6 (E P2B extras, VTARGET=0).
  mk_g6.py G6W|G6P|Ex6 pilot|prod   -> job lists $P/runs/ours6/<SET>[-pilot]/<part>/jobs.list + groups.txt"""
import subprocess, sys, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
BIN = P + '/bin6'; setn, mode = sys.argv[1], sys.argv[2]
soft = P + '/tables/soft.tab'
FP = ['0.1', '5000']
M1 = ['1d-9', '1'] + FP + ['psmc', '0', '0', '1d-12']
M2 = ['1d-9', '2', '0', '0', 'psmc', '0', '0', '1d-12']
# name: (binary, part, ncall, seed0, nseeds, flags, hours)
G6 = {
 'b1': ('sliced21', 'b1', 3000000, 2001, 100, [], 4), 'b2': ('sliced21', 'b2', 5000000, 2001, 100, [], 5),
 'lo': ('nlo31', 'lo', 3000000, 2001, 100, [], 2), 'vi': ('nlo31', 'vi', 1000000, 2001, 250, [], 10),
 'kp': ('nlo31', 'kp', 2000000, 2001, 200, [], 2), 'r': ('nlo31', 'r', 1000000, 1001, 2000, ['EQ'], 10),
 'rcorr': ('nlo31', 'r', 1000000, 2001, 300, ['EQ', 'TD'], 10), 'd2': ('nlo31', 'r', 1000000, 2001, 150, ['EQ', 'TD2'], 10)}
EX = {'lo-x': ('nlo31', 'lo', 3000000, 2001, 200, [], 1), 'vi-x': ('nlo31', 'vi', 1000000, 2001, 250, [], 6),
      'kp-x': ('nlo31', 'kp', 2000000, 2001, 200, [], 1), 'r-x': ('nlo31', 'r', 1000000, 1001, 2000, ['EQ'], 6)}
if setn == 'Ex6':
    parts, env0 = EX, ['ZFIX=3', 'P2BSLICE=1', 'P2BEXTRA=1', 'VTARGET=0']; ZP = P + '/tables/zeus5'
else:
    parts, env0 = G6, ['EW31="%s"' % ('2 0' if setn == 'G6W' else '0 0')]
base = P + '/runs/ours6/' + setn + ('-pilot' if mode == 'pilot' else '')
out = []
for name, (b, part, nc, s0, ns, flags, th) in parts.items():
    env = list(env0)
    if 'EQ' in flags: env += ['VEGAS_EQUAL=2']
    if 'TD' in flags: env += ['TECHDIFF=1', 'TECHMIN="1d-10 1d-8"']
    if 'TD2' in flags: env += ['TECHDIFF=1', 'TECHMIN="1d-11 1d-9"', 'TECHREF="1d-10 1d-8"']
    if mode == 'pilot':
        if setn == 'Ex6' and name != 'r-x': continue
        s0, n = 9001, (20 if setn == 'Ex6' else 4)
    else: n = ns
    if setn == 'Ex6':
        cmd = [BIN + '/nlo31', part, str(nc), '6', '{seed}'] + M2
    elif b == 'sliced21':
        cmd = [BIN + '/sliced21', part, str(nc), '6', '{seed}'] + FP + [soft]
    else:
        cmd = [BIN + '/nlo31', part, str(nc), '6', '{seed}'] + M1
    d = '%s/%s' % (base, name); os.makedirs(d, exist_ok=True)
    a = [S + '/mkjobs.py', d + '/jobs.list', d, '%d-%d' % (s0, s0 + n - 1), '--pattern', 'RESULT']
    for e in env: a += ['--env', e]
    subprocess.run(a + ['--'] + cmd, check=True, stdout=subprocess.DEVNULL)
    tag = {'G6W': 'gw', 'G6P': 'gp', 'Ex6': 'x6'}[setn] + ('p' if mode == 'pilot' else '')
    out.append('dis-%s-%s %s/jobs.list %d:00:00 1000MB 1 et[07-40]' % (tag, name, d, 6 if mode == 'pilot' else th))
open(base + '/groups.txt', 'w').write('\n'.join(out) + '\n'); print(base + '/groups.txt', len(out), 'groups')
