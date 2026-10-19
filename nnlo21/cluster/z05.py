#!/usr/bin/env python3
"""Z05 (UPDATE 18 Oct): N3LO 1+1 P2B, ZEUS set-up of 1803.09973, LABSET=zeus05 SCALES=3, bin12.
  z05.py init                 write $Z/<part>/cmd (command template, SEED replaced by run_z05.sh)
  z05.py submit [--dry]       submit every seed up to the part's cap ($Z/caps: 'part maxseed') that has no .out
                              and is neither pending nor running (array task = seed - 60000, job name dis-z05-<part>)
  z05.py status               per part: seeds in range, cap, done, active, failed attempts
One file per finished job: $Z/<part>/<part>_<seed>.out (see run_z05.sh). Caps (as submit.py): user total < USERCAP,
own dis-* < DISCAP. A seed with >= MAXATT attempts is left alone."""
import os, re, subprocess, sys
from collections import Counter
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; Z = P + '/Z05'; B = P + '/bin12'; T = P + '/tables'
HERE = os.path.dirname(os.path.abspath(__file__))
USERCAP = int(os.environ.get('USERCAP', 19000)); DISCAP = int(os.environ.get('DISCAP', 11500)); MAXATT = int(os.environ.get('MAXATT', 4))
EXCL = 'et[07-40],ct30,gt34,gt37,kt02,kt06,gt25,gt27,gt28,gt29,gt39,gt40,kt01,kt03,kt10'
ENV = 'export LABSET=zeus05 SCALES=3; '
RV = 'export VEGAS_EQUAL=2 VTARGET=0 VEGAS_MIX=0.2; '
M3 = '1d-9 3 0 0 psmc 0 0 1d-12'
SOFT = T + '/softG.tab'
# part: (command, first seed, last seed, time limit)
PARTS = {
    'r':   (RV + B + '/nlo31 r 2000000 6 SEED ' + M3, 60001, 78700, '8:00:00'),
    'vi':  (B + '/nlo31 vi 1000000 6 SEED ' + M3, 80001, 82729, '8:00:00'),
    'kp':  (B + '/nlo31 kp 2000000 6 SEED ' + M3, 83001, 85479, '2:00:00'),
    'b2':  (B + '/sliced21 b2 5000000 6 SEED p2b %s/z05 %s' % (T, SOFT), 86001, 86729, '10:00:00'),
    'rc1': (RV + "export TECHDIFF=1 TECHMIN='1d-10 1d-8'; " + B + '/nlo31 r 1000000 6 SEED ' + M3, 87001, 88000, '4:00:00'),
    'rc2': (RV + "export TECHDIFF=1 TECHMIN='1d-11 1d-9' TECHREF='1d-10 1d-8'; " + B + '/nlo31 r 1000000 6 SEED 1d-9 3 0 0 psmc 0 0 1d-14', 88001, 89000, '4:00:00'),
    'b0':  (B + '/sliced21 b0 10000000 6 SEED p2b none', 89001, 89200, '2:00:00'),
    'b1':  (B + '/sliced21 b1 3000000 6 SEED p2b %s/z05 %s' % (T, SOFT), 89201, 89600, '6:00:00'),
    'lo':  (B + '/nlo31 lo 20000000 6 SEED ' + M3, 89601, 90000, '4:00:00'),
}
def sh(c): return subprocess.run(c, check=True, text=True, capture_output=True).stdout
def ranges(v):
    out, i = [], 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[j + 1] == v[j] + 1: j += 1
        out.append(str(v[i]) if i == j else '%d-%d' % (v[i], v[j])); i = j + 1
    return out
def caps():
    c = {}
    if os.path.exists(Z + '/caps'):
        for l in open(Z + '/caps'):
            w = l.split()
            if len(w) == 2 and not w[0].startswith('#'): c[w[0]] = int(w[1])
    return c
def state(p):
    d = Z + '/' + p
    done = {int(m.group(1)) for f in os.listdir(d) for m in [re.match(r'%s_(\d+)\.out$' % p, f)] if m}
    att = Counter()
    if os.path.exists(d + '/attempts'):
        for l in open(d + '/attempts'):
            w = l.split()
            if w: att[int(w[0])] += 1
    return done, att
def queue():
    q = sh(['squeue', '-u', os.environ['USER'], '-h', '-r', '-o', '%j %K %t']).splitlines()
    act = {}
    for l in q:
        f = l.split()
        if f[0].startswith('dis-z05-') and f[2] in ('PD', 'R', 'CF') and f[1].isdigit(): act.setdefault(f[0][8:], set()).add(int(f[1]) + 60000)
    return q, act
cmd = sys.argv[1] if len(sys.argv) > 1 else 'status'
if cmd == 'init':
    for p, (c, a, b, t) in PARTS.items():
        os.makedirs(Z + '/' + p, exist_ok=True)
        open(Z + '/' + p + '/cmd', 'w').write('set -e; ' + ENV + c + '\n')
    print('wrote', ', '.join('%s/cmd' % p for p in PARTS))
elif cmd == 'status':
    C = caps(); q, act = queue()
    for p, (c, a, b, t) in PARTS.items():
        done, att = state(p); n = b - a + 1
        A = act.get(p, set()); fail = sum(max(0, att[s] - (1 if s in done or s in A else 0)) for s in att)
        print('%-4s seeds %5d-%5d (%5d) cap %5d  done %5d  active %5d  failed attempts %d' % (p, a, b, n, C.get(p, a - 1), len(done), len(act.get(p, ())), fail))
elif cmd == 'submit':
    dry = '--dry' in sys.argv; C = caps(); q, act = queue()
    ndis = sum(1 for l in q if l.split()[0].startswith('dis-'))
    room = min(USERCAP - len(q), DISCAP - ndis)
    # r last: the small parts first, so every part is complete early (r is 90 % of the jobs)
    for p, (c, a, b, t) in sorted(PARTS.items(), key=lambda x: x[0] == 'r'):
        cap = min(C.get(p, a - 1), b)
        done, att = state(p); A = act.get(p, set())
        todo = [s for s in range(a, cap + 1) if s not in done and s not in A and att[s] < MAXATT]
        stuck = sum(1 for s in range(a, cap + 1) if s not in done and att[s] >= MAXATT)
        msg = '%s: cap %d, done %d, active %d, stuck %d, to submit %d' % (p, cap, len(done), len(A), stuck, len(todo))
        todo = todo[:max(0, room)]
        if not todo: print(msg + ('; queue at cap' if room <= 0 else '')); continue
        rr = ranges([s - 60000 for s in todo])
        # sbatch array specs: submit in chunks of at most 200 ranges
        for k in range(0, len(rr), 200):
            spec = ','.join(rr[k:k + 200])
            sc = ['sbatch', '--parsable', '--nice=100', '--array=' + spec, '--time=' + t, '--mem=1000MB', '--cpus-per-task=1',
                  '-J', 'dis-z05-' + p, '--export=ALL,PART=' + p, '-o', '/dev/null', '--exclude=' + EXCL, HERE + '/run_z05.sh']
            if dry: print(' '.join(sc)); continue
            jid = sh(sc).strip(); msg += '; array %s' % jid
        room -= len(todo); print(msg)
        if not dry:
            with open(Z + '/submitted.log', 'a') as f: f.write(sh(['date', '+%F %T']).strip() + ' ' + msg + '\n')
