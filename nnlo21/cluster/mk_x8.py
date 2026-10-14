#!/usr/bin/env python3
"""14 Oct, X8 (UPDATE 14 Oct item 2): Ex7 commands and bin7 with VTARGET=53 (VEGAS aimed at the total row,
tau_cut 5e-3: cell = column + 10 x (row - 1)), seeds 8001-: lo-x 200, vi-x 250, kp-x 200, r-x 400.
  mk_x8.py  -> $P/runs/ours7/X8/<part>/jobs.list + $P/runs/ours7/X8/groups.txt"""
import subprocess, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
BIN = P + '/bin7'; M2 = ['1d-9', '2', '0', '0', 'psmc', '0', '0', '1d-12']
EX = {'lo-x': ('lo', 3000000, 200, 1, []), 'vi-x': ('vi', 1000000, 250, 6, []),
      'kp-x': ('kp', 2000000, 200, 1, []), 'r-x': ('r', 1000000, 400, 6, ['EQ'])}
base = P + '/runs/ours7/X8'; out = []
for name, (part, nc, n, th, flags) in EX.items():
    env = ['ZFIX=3', 'P2BSLICE=1', 'P2BEXTRA=1', 'VTARGET=53'] + (['VEGAS_EQUAL=2'] if 'EQ' in flags else [])
    d = '%s/%s' % (base, name); os.makedirs(d, exist_ok=True)
    a = [S + '/mkjobs.py', d + '/jobs.list', d, '8001-%d' % (8000 + n), '--pattern', 'RESULT']
    for e in env: a += ['--env', e]
    subprocess.run(a + ['--', BIN + '/nlo31', part, str(nc), '6', '{seed}'] + M2, check=True, stdout=subprocess.DEVNULL)
    out.append('dis-x8-%s %s/jobs.list %d:00:00 1000MB 1 et[07-40]' % (name, d, th))
open(base + '/groups.txt', 'w').write('\n'.join(out) + '\n'); print(base + '/groups.txt', len(out), 'groups')
