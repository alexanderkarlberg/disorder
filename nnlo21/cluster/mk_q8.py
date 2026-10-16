#!/usr/bin/env python3
"""16 Oct round Q8 (UPDATE 16 Oct): top-up of the CH8 quark channel (mask 1, zeus5q, J8 window and commands, bin10) in the parts that
carry the error below tau_cut 1e-4: r, rcorr, d2. Seeds 21001- for every part (disjoint from all CH8 seeds).
Counts from the CH8 q per-seed sd of N(1e-4) - N(1e-5) (total row; sizing in NOTES.md, 16 Oct) for a slope error <= 6e-9 per unit ln, +~15%.
  mk_q8.py   -> $P/runs/ours8j/Q8/<part>/jobs.list, $P/runs/ours8j/groups-q8.txt (groups dis-q8-<part>)"""
import subprocess, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
BIN = P + '/bin10'; EXCL = 'et[07-40],gt25,gt27,gt28,gt29,gt39,gt40,kt01,kt03,kt10'
CNT = dict(r=2400, rcorr=2200, d2=1600); S0 = 21001
# name: (part, ncall, extra env, hours), as mk_ch8.py
PARTS = {'r': ('r', 1000000, [], 10), 'rcorr': ('r', 1000000, ['TECHDIFF=1', 'TECHMIN="1d-10 1d-8"'], 10),
         'd2': ('r', 1000000, ['TECHDIFF=1', 'TECHMIN="1d-11 1d-9"', 'TECHREF="1d-10 1d-8"'], 10)}
base = P + '/runs/ours8j/Q8'; full = []
for name, (part, nc, extra, th) in PARTS.items():
    env = ['ZFIX=2', 'ZCUTS=1', 'EW31="0 0"', 'P2BSLICE=1', 'VEGAS_EQUAL=2'] + extra
    cmd = [BIN + '/nlo31', part, str(nc), '6', '{seed}', '1d-9', '2', '0.1', '5000', 'psmc', '1', '0', '1d-12']
    d = '%s/%s' % (base, name); os.makedirs(d, exist_ok=True)
    a = [S + '/mkjobs.py', d + '/jobs.list', d, '%d-%d' % (S0, S0 + CNT[name] - 1), '--pattern', 'RESULT']
    for e in env: a += ['--env', e]
    subprocess.run(a + ['--'] + cmd, check=True, stdout=subprocess.DEVNULL)
    full.append('dis-q8-%s %s/jobs.list %d:00:00 1000MB 1 %s' % (name, d, th, EXCL))
open(P + '/runs/ours8j/groups-q8.txt', 'w').write('\n'.join(full) + '\n')
print(base, CNT, 'core-h ~%.0f' % (CNT['r']*2.03 + CNT['rcorr']*0.61 + CNT['d2']*0.60))
