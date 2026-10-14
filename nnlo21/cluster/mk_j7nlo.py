#!/usr/bin/env python3
"""14 Oct: NLO pieces (b1, lo) in the J7 window, which UPDATE 13 Oct (c) item 4 needs as 'the NLO' but item 2 did not list.
  mk_j7nlo.py PH|W -> $P/runs/ours7j/<SET>/{b1,lo}/jobs.list (100 seeds each, seeds as the main parts: PH 1001-, W 5001-),
  <SET>/groups-nlo.txt. Same settings as the bin5 b1/lo of the validation sets (b1: sliced21 b1 3M, lo: nlo31 lo 3M, P2BSLICE=1) with ZFIX=2, bin6, window 2 0.01 400."""
import subprocess, sys, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'; BIN = P + '/bin6'; setn = sys.argv[1]
C = {'PH': dict(s0=1001, tab='zeus5', ew='0 0', tag='jp'), 'W': dict(s0=5001, tab='zeus5w', ew='2 0', tag='jw')}[setn]
base = '%s/runs/ours7j/%s' % (P, setn); full = []
for name, th in (('b1', 3), ('lo', 2)):
    d = '%s/%s' % (base, name); os.makedirs(d, exist_ok=True)
    if name == 'b1': env = ['ZFIX=2', 'EW31="%s"' % C['ew']]; cmd = [BIN + '/sliced21', 'b1', '3000000', '6', '{seed}', 'zeus', P + '/tables/' + C['tab'], P + '/tables/soft.tab']
    else: env = ['ZFIX=2', 'EW31="%s"' % C['ew'], 'P2BSLICE=1']; cmd = [BIN + '/nlo31', 'lo', '3000000', '6', '{seed}', '1d-9', '2', '0.01', '400', 'psmc', '0', '0', '1d-12']
    a = [S + '/mkjobs.py', d + '/jobs.list', d, '%d-%d' % (C['s0'], C['s0'] + 99), '--pattern', 'RESULT']
    for e in env: a += ['--env', e]
    subprocess.run(a + ['--'] + cmd, check=True, stdout=subprocess.DEVNULL)
    full.append('dis-%s-%s %s/jobs.list %d:00:00 1000MB 1 et[07-40]' % (C['tag'], name, d, th))
open(base + '/groups-nlo.txt', 'w').write('\n'.join(full) + '\n')
