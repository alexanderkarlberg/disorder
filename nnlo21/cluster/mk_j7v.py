#!/usr/bin/env python3
"""14 Oct, J7v (UPDATE 14 Oct (c) + (d)): J7 PH commands with bin9, ZCUTS=1, seeds 3001-:
 vi 240, kp 200, r 1600 with P2BSLICE=1 VTARGET=-1 VEGAS_EQUAL=2; reruns b2 200, rcorr 200, d2 100, b1 100, lo 100 (J7 PH settings, VTARGET unset).
  mk_j7v.py  -> $P/runs/ours7j/PHv9/<part>/{jobs.list,pilot.list} + groups.txt, groups-pilot.txt
(the bin7 pilots of J7v vi/kp/r, seeds 3001-3005, sit in $P/runs/ours7j/PHv and are only a cost estimate)"""
import subprocess, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
BIN = P + '/bin9'; W = ['1d-9', '2', '0.01', '400', 'psmc', '0', '0', '1d-12']
soft = P + '/tables/soft.tab'; tab = P + '/tables/zeus5'
# name: (binary, part, ncall, n, hours, extra env)
V = ['VTARGET=-1', 'VEGAS_EQUAL=2']
PARTS = {'vi': ('nlo31', 'vi', 1000000, 240, 10, ['P2BSLICE=1'] + V), 'kp': ('nlo31', 'kp', 2000000, 200, 2, ['P2BSLICE=1'] + V),
         'r': ('nlo31', 'r', 1000000, 1600, 10, ['P2BSLICE=1'] + V),
         'b2': ('sliced21', 'b2', 5000000, 200, 5, []),
         'rcorr': ('nlo31', 'r', 1000000, 200, 10, ['P2BSLICE=1', 'VEGAS_EQUAL=2', 'TECHDIFF=1', 'TECHMIN="1d-10 1d-8"']),
         'd2': ('nlo31', 'r', 1000000, 100, 10, ['P2BSLICE=1', 'VEGAS_EQUAL=2', 'TECHDIFF=1', 'TECHMIN="1d-11 1d-9"', 'TECHREF="1d-10 1d-8"']),
         'b1': ('sliced21', 'b1', 3000000, 100, 3, []), 'lo': ('nlo31', 'lo', 3000000, 100, 2, ['P2BSLICE=1'])}
base = P + '/runs/ours7j/PHv9'; full, pil = [], []
for name, (b, part, nc, n, th, ex) in PARTS.items():
    env = ['ZFIX=2', 'EW31="0 0"', 'ZCUTS=1'] + ex
    cmd = [BIN + '/sliced21', part, str(nc), '6', '{seed}', 'zeus', tab, soft] if b == 'sliced21' else [BIN + '/nlo31', part, str(nc), '6', '{seed}'] + W
    d = '%s/%s' % (base, name); os.makedirs(d, exist_ok=True)
    for lst, k in ((d + '/jobs.list', n), (d + '/pilot.list', 5)):
        a = [S + '/mkjobs.py', lst, d, '3001-%d' % (3000 + k), '--pattern', 'RESULT']
        for e in env: a += ['--env', e]
        subprocess.run(a + ['--'] + cmd, check=True, stdout=subprocess.DEVNULL)
    full.append('dis-jv9-%s %s/jobs.list %d:00:00 1000MB 1 et[07-40]' % (name, d, th))
    pil.append('dis-jv9p-%s %s/pilot.list 6:00:00 1000MB 1 et[07-40]' % (name, d))
open(base + '/groups.txt', 'w').write('\n'.join(full) + '\n'); open(base + '/groups-pilot.txt', 'w').write('\n'.join(pil) + '\n')
print(base, len(full), 'groups')
