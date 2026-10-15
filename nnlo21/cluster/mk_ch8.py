#!/usr/bin/env python3
"""15 Oct round CH8 (UPDATE 15 Oct (b)): the J8 window (ZEUS jets at (x, Q2) = (0.1, 5000), photon, ZFIX=2 ZCUTS=1) split by the
incoming parton: mask 1 = quarks + antiquarks (CH8q, tables zeus5q, seeds 7001-), mask 2 = gluon (CH8g, zeus5g, seeds 8001-).
bin10. Commands as mk_j8.py / mk_j8nlo.py (J8, bin9) plus the channel mask: sliced21 args 8, 9 = `0.02 <mask>`, nlo31 argument 10 = mask.
  mk_ch8.py [q] [g]   -> $P/runs/ours8j/CH8<ch>/<part>/{jobs.list,pilot.list} + CH8<ch>/groups.txt, groups-pilot.txt
pilot.list = the first 5 seeds of each part (they are production seeds, they count).
Groups: dis-ch8<ch>-<part> (production), dis-ch8<ch>p-<part> (pilot)."""
import subprocess, sys, os
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
BIN = P + '/bin10'; soft = P + '/tables/soft.tab'
EXCL = 'et[07-40],gt25,gt27,gt28,gt29,gt39,gt40,kt01,kt03,kt10'
CH = {'q': dict(mask=1, tab='zeus5q', s0=7001), 'g': dict(mask=2, tab='zeus5g', s0=8001)}
CNT = dict(b2=200, vi=240, kp=200, r=1600, rcorr=200, d2=100, b1=100, lo=100)
# name: (binary, part, ncall, flags, hours)
PARTS = {'b2': ('sliced21', 'b2', 5000000, [], 5), 'vi': ('nlo31', 'vi', 1000000, [], 10), 'kp': ('nlo31', 'kp', 2000000, [], 2),
         'r': ('nlo31', 'r', 1000000, ['EQ'], 10), 'rcorr': ('nlo31', 'r', 1000000, ['EQ', 'TD'], 10),
         'd2': ('nlo31', 'r', 1000000, ['EQ', 'TD2'], 10),
         'b1': ('sliced21', 'b1', 3000000, [], 3), 'lo': ('nlo31', 'lo', 3000000, [], 2)}
for ch in (sys.argv[1:] or ['q', 'g']):
    C = CH[ch]; m = str(C['mask']); base = '%s/runs/ours8j/CH8%s' % (P, ch); full, pil = [], []
    for name, (b, part, nc, flags, th) in PARTS.items():
        env = ['ZFIX=2', 'ZCUTS=1', 'EW31="0 0"']
        if b == 'sliced21': env += ['ZWIN="0.1 5000"']
        if 'EQ' in flags: env += ['P2BSLICE=1', 'VEGAS_EQUAL=2']
        if b == 'nlo31' and 'EQ' not in flags: env += ['P2BSLICE=1']
        if 'TD' in flags: env += ['TECHDIFF=1', 'TECHMIN="1d-10 1d-8"']
        if 'TD2' in flags: env += ['TECHDIFF=1', 'TECHMIN="1d-11 1d-9"', 'TECHREF="1d-10 1d-8"']
        if b == 'sliced21': cmd = [BIN + '/sliced21', part, str(nc), '6', '{seed}', 'zeus', P + '/tables/' + C['tab'], soft, '0.02', m]
        else: cmd = [BIN + '/nlo31', part, str(nc), '6', '{seed}', '1d-9', '2', '0.1', '5000', 'psmc', m, '0', '1d-12']
        d = '%s/%s' % (base, name); os.makedirs(d, exist_ok=True); s0 = C['s0']
        for lst, n in ((d + '/jobs.list', CNT[name]), (d + '/pilot.list', 5)):
            a = [S + '/mkjobs.py', lst, d, '%d-%d' % (s0, s0 + n - 1), '--pattern', 'RESULT']
            for e in env: a += ['--env', e]
            subprocess.run(a + ['--'] + cmd, check=True, stdout=subprocess.DEVNULL)
        full.append('dis-ch8%s-%s %s/jobs.list %d:00:00 1000MB 1 %s' % (ch, name, d, th, EXCL))
        pil.append('dis-ch8%sp-%s %s/pilot.list 6:00:00 1000MB 1 %s' % (ch, name, d, EXCL))
    open(base + '/groups.txt', 'w').write('\n'.join(full) + '\n'); open(base + '/groups-pilot.txt', 'w').write('\n'.join(pil) + '\n')
    print(base, len(full), 'groups')
if len(sys.argv) == 1:
    b = P + '/runs/ours8j'
    open(b + '/groups-ch8.txt', 'w').write(open(b + '/CH8q/groups.txt').read() + open(b + '/CH8g/groups.txt').read())
    open(b + '/groups-ch8-pilot.txt', 'w').write(open(b + '/CH8q/groups-pilot.txt').read() + open(b + '/CH8g/groups-pilot.txt').read())
