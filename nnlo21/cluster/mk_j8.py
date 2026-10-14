#!/usr/bin/env python3
"""14 Oct round J8 (UPDATE 14 Oct (b)): ZEUS jets in the window around (x, Q2) = (0.1, 5000), ZFIX=2;
UPDATE 14 Oct (d): bin9 with ZCUTS=1 (cut-design rows), run dirs $P/runs/ours8j/<SET>9 (the bin8 set <SET> was pilots only); sliced21 jobs ZWIN="0.1 5000", nlo31 window args 2 0.1 5000.
  mk_j8.py PH        -> $P/runs/ours8j/<SET>/<part>/{jobs.list,pilot.list} + <SET>/groups.txt (full) and groups-pilot.txt
PH: EW31="0 0", tables zeus5, seeds 1001-; W: EW31="2 0", tables zeus5w, seeds 5001-.
pilot.list = the first 5 seeds of each part (they are production seeds, they count)."""
import subprocess, sys, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
BIN = P + '/bin9'; setn = sys.argv[1]
soft = P + '/tables/soft.tab'
W = ['1d-9', '2', '0.1', '5000', 'psmc', '0', '0', '1d-12']
CNT = {'PH': dict(b2=200, vi=240, kp=200, r=1600, rcorr=200, d2=100, s0=1001, tab='zeus5', ew='0 0', tag='jp'),
       'W': dict(b2=150, vi=180, kp=150, r=1200, rcorr=200, d2=100, s0=5001, tab='zeus5w', ew='2 0', tag='jw')}[setn]
# name: (binary, part, ncall, flags, hours)
PARTS = {'b2': ('sliced21', 'b2', 5000000, [], 5), 'vi': ('nlo31', 'vi', 1000000, [], 10), 'kp': ('nlo31', 'kp', 2000000, [], 2),
         'r': ('nlo31', 'r', 1000000, ['EQ'], 10), 'rcorr': ('nlo31', 'r', 1000000, ['EQ', 'TD'], 10),
         'd2': ('nlo31', 'r', 1000000, ['EQ', 'TD2'], 10)}
base = '%s/runs/ours8j/%s9' % (P, setn); full, pil = [], []
for name, (b, part, nc, flags, th) in PARTS.items():
    env = ['ZFIX=2', 'ZCUTS=1', 'EW31="%s"' % CNT['ew']]
    if b == 'sliced21': env += ['ZWIN="0.1 5000"']
    if 'EQ' in flags: env += ['P2BSLICE=1', 'VEGAS_EQUAL=2'] if b == 'nlo31' else []
    if b == 'nlo31' and 'EQ' not in flags: env += ['P2BSLICE=1']
    if 'TD' in flags: env += ['TECHDIFF=1', 'TECHMIN="1d-10 1d-8"']
    if 'TD2' in flags: env += ['TECHDIFF=1', 'TECHMIN="1d-11 1d-9"', 'TECHREF="1d-10 1d-8"']
    if b == 'sliced21': cmd = [BIN + '/sliced21', part, str(nc), '6', '{seed}', 'zeus', P + '/tables/' + CNT['tab'], soft]
    else: cmd = [BIN + '/nlo31', part, str(nc), '6', '{seed}'] + W
    d = '%s/%s' % (base, name); os.makedirs(d, exist_ok=True)
    s0 = CNT['s0']
    for lst, n in ((d + '/jobs.list', CNT[name]), (d + '/pilot.list', 5)):
        a = [S + '/mkjobs.py', lst, d, '%d-%d' % (s0, s0 + n - 1), '--pattern', 'RESULT']
        for e in env: a += ['--env', e]
        subprocess.run(a + ['--'] + cmd, check=True, stdout=subprocess.DEVNULL)
    full.append('dis-j9%s-%s %s/jobs.list %d:00:00 1000MB 1 et[07-40]' % (CNT['tag'], name, d, th))
    pil.append('dis-j9%sp-%s %s/pilot.list 6:00:00 1000MB 1 et[07-40]' % (CNT['tag'], name, d))
open(base + '/groups.txt', 'w').write('\n'.join(full) + '\n'); open(base + '/groups-pilot.txt', 'w').write('\n'.join(pil) + '\n')
print(base, len(full), 'groups')
