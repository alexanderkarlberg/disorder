#!/usr/bin/env python3
"""15 Oct round G8 (UPDATE 15 Oct): the J8 window (ZEUS jets, (x,Q2)=(0.1,5000), ZCUTS=1, bin9) with nlo31's argument 9 `flat` instead of `psmc`, parts vi, kp, r.
  mk_g8.py NSEEDS [NSEEDS_X4]  -> $P/runs/ours8j/G8/<part> (ncall as J8: vi 1e6, kp 2e6, r 1e6) and $P/runs/ours8j/G8x4/<part> (4x ncall), seeds 6001-
  Writes jobs.list in each part dir, groups.txt (and groups-pilot.txt) in G8/ and G8x4/. Pilot = the first 5 seeds (they count for production).
  Seeds are appended to jobs.list idempotently (mkjobs.py), so the script can be rerun with larger numbers.
  NSEEDS is a comma list vi,kp,r (e.g. 240,200,1600); use 5 for the pilot."""
import subprocess, sys, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
BIN = P + '/bin9'
W = ['1d-9', '2', '0.1', '5000', 'flat', '0', '0', '1d-12']
def cnt(a): v = [int(x) for x in a.split(',')]; return dict(zip(('vi', 'kp', 'r'), v * 3 if len(v) == 1 else v))
N = cnt(sys.argv[1]); N4 = cnt(sys.argv[2]) if len(sys.argv) > 2 else N
NC = {'vi': 1000000, 'kp': 2000000, 'r': 1000000}
TH = {'G8': dict(vi=14, kp=3, r=12), 'G8x4': dict(vi=24, kp=24, r=24)}
for setn, mult, nn in (('G8', 1, N), ('G8x4', 4, N4)):
    base = '%s/runs/ours8j/%s' % (P, setn); full, pil = [], []
    for part in ('vi', 'kp', 'r'):
        env = ['ZFIX=2', 'ZCUTS=1', 'EW31="0 0"', 'P2BSLICE=1'] + (['VEGAS_EQUAL=2'] if part == 'r' else [])
        cmd = [BIN + '/nlo31', part, str(NC[part] * mult), '6', '{seed}'] + W
        d = '%s/%s' % (base, part); os.makedirs(d, exist_ok=True)
        for lst, n in ((d + '/jobs.list', nn[part]), (d + '/pilot.list', 5)):
            a = [S + '/mkjobs.py', lst, d, '6001-%d' % (6000 + n), '--pattern', 'RESULT']
            for e in env: a += ['--env', e]
            subprocess.run(a + ['--'] + cmd, check=True, stdout=subprocess.DEVNULL)
        full.append('dis-g%s-%s %s/jobs.list %d:00:00 1000MB 1 et[07-40],gt25,gt27,gt28,gt29,gt39,gt40,kt01,kt03,kt10' % (setn[1:].lower(), part, d, TH[setn][part]))
        pil.append('dis-g%sp-%s %s/pilot.list %d:00:00 1000MB 1 et[07-40],gt25,gt27,gt28,gt29,gt39,gt40,kt01,kt03,kt10' % (setn[1:].lower(), part, d, TH[setn][part]))
    open(base + '/groups.txt', 'w').write('\n'.join(full) + '\n'); open(base + '/groups-pilot.txt', 'w').write('\n'.join(pil) + '\n')
    print(base, nn)
