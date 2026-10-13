#!/usr/bin/env python3
"""13 Oct round, bin7: Ex7 = Ex6 (E P2B extras) with VTARGET=-1, new seeds (5001-).
  mk_ex7.py pilot|prod   -> $P/runs/ours7/Ex7[-pilot]/<part>/jobs.list + groups.txt
Seeds are disjoint from Ex6 (lo/vi/kp 2001-, r 1001-3000), the bin5 E extras (1001-1875, 2001-2062) and the
Ex6 pilot (9001-9020); the pilot here uses 9101-9120."""
import subprocess, sys, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
BIN = P + '/bin7'; mode = sys.argv[1]
M2 = ['1d-9', '2', '0', '0', 'psmc', '0', '0', '1d-12']
# name: (part, ncall, seed0, nseeds, flags, hours)
EX = {'lo-x': ('lo', 3000000, 5001, 200, [], 1), 'vi-x': ('vi', 1000000, 5001, 250, [], 6),
      'kp-x': ('kp', 2000000, 5001, 200, [], 1), 'r-x': ('r', 1000000, 5001, 2000, ['EQ'], 6)}
env0 = ['ZFIX=3', 'P2BSLICE=1', 'P2BEXTRA=1', 'VTARGET=-1']
base = P + '/runs/ours7/Ex7' + ('-pilot' if mode == 'pilot' else '')
out = []
for name, (part, nc, s0, n, flags, th) in EX.items():
    if mode == 'pilot':
        if name != 'r-x': continue
        s0, n = 9101, 20
    env = list(env0) + (['VEGAS_EQUAL=2'] if 'EQ' in flags else [])
    cmd = [BIN + '/nlo31', part, str(nc), '6', '{seed}'] + M2
    d = '%s/%s' % (base, name); os.makedirs(d, exist_ok=True)
    a = [S + '/mkjobs.py', d + '/jobs.list', d, '%d-%d' % (s0, s0 + n - 1), '--pattern', 'RESULT']
    for e in env: a += ['--env', e]
    subprocess.run(a + ['--'] + cmd, check=True, stdout=subprocess.DEVNULL)
    out.append('dis-x7%s-%s %s/jobs.list %d:00:00 1000MB 1 et[07-40]' % ('p' if mode == 'pilot' else '', name, d, 6 if mode == 'pilot' else th))
open(base + '/groups.txt', 'w').write('\n'.join(out) + '\n'); print(base + '/groups.txt', len(out), 'groups')
