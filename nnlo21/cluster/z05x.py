#!/usr/bin/env python3
"""Z05x (UPDATE 20 Oct): cost study of the P2B extras, sharp and smoothed (ZSMEAR), bin13, LABSET=zeus05.
  z05x.py init                write $Z/<part>/cmd (command template, SEED replaced) and $Z/jobs.list
  z05x.py submit [--dry]      submit the jobs up to the caps in $Z/caps ('part maxseed') that have no .out and are neither
                              pending nor running; array task = line number of $Z/jobs.list (see run_z05x.sh)
  z05x.py status              per part: seeds, cap, done, active, failed attempts
jobs.list holds 'part seed' lines in the order rx s, rxs s, rxt s, ... seed by seed, so that any partial set gives pairs
(a sharp part and its smoothed ones with the same seed). One file per job: $Z/<part>/<part>_<seed>.out. Caps as z05.py."""
import os, re, subprocess, sys
from collections import Counter
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; Z = P + '/Z05x'; B = P + '/bin13'
HERE = os.path.dirname(os.path.abspath(__file__))
USERCAP = int(os.environ.get('USERCAP', 19000)); DISCAP = int(os.environ.get('DISCAP', 11500)); MAXATT = int(os.environ.get('MAXATT', 4))
EXCL = 'et[07-40],ct30,gt09,gt34,gt37,kt02,kt06,gt25,gt27,gt28,gt29,gt39,gt40,kt01,kt03,kt10'
ENV = 'export LABSET=zeus05; '
X = 'export P2BSLICE=1 P2BEXTRA=1 VTARGET=-1; '
RV = 'export VEGAS_EQUAL=2 VEGAS_MIX=0.2; '
M3 = '1d-9 3 0 0 psmc 0 0 1d-12'
def zs(e): return 'export ZSMEAR=%s; ' % e if e else ''
# part: (command, first seed, last seed)
PARTS = {
    'rx':   (X + RV + B + '/nlo31 r 2000000 6 SEED ' + M3, 6001, 10000),
    'rxs':  (X + RV + zs('0.1') + B + '/nlo31 r 2000000 6 SEED ' + M3, 6001, 10000),
    'rxt':  (X + RV + zs('0.3') + B + '/nlo31 r 2000000 6 SEED ' + M3, 6001, 7000),
    'vix':  (X + B + '/nlo31 vi 1000000 6 SEED ' + M3, 14001, 16000),
    'vixs': (X + zs('0.1') + B + '/nlo31 vi 1000000 6 SEED ' + M3, 14001, 16000),
    'kpx':  (X + B + '/nlo31 kp 2000000 6 SEED ' + M3, 16001, 17000),
    'kpxs': (X + zs('0.1') + B + '/nlo31 kp 2000000 6 SEED ' + M3, 16001, 17000),
}
def sh(c): return subprocess.run(c, check=True, text=True, capture_output=True).stdout
def jobs():
    L = []
    for s in sorted({s for p, (c, a, b) in PARTS.items() for s in range(a, b + 1)}):
        for p, (c, a, b) in PARTS.items():       # PARTS order: rx, rxs, rxt, ...
            if a <= s <= b: L.append((p, s))
    return L
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
    act = set()
    for l in q:
        f = l.split()
        if f[0].startswith('dis-z05x') and f[2] in ('PD', 'R', 'CF') and f[1].isdigit(): act.add(int(f[1]))
    return q, act
cmd = sys.argv[1] if len(sys.argv) > 1 else 'status'
if cmd == 'init':
    for p, (c, a, b) in PARTS.items():
        os.makedirs(Z + '/' + p, exist_ok=True)
        open(Z + '/' + p + '/cmd', 'w').write('set -e; ' + ENV + c + '\n')
    open(Z + '/jobs.list', 'w').write(''.join('%s %d\n' % j for j in jobs()))
    print('wrote', ', '.join('%s/cmd' % p for p in PARTS), 'and jobs.list (%d lines)' % len(jobs()))
elif cmd == 'status':
    C = caps(); q, act = queue(); L = jobs(); line = {j: i + 1 for i, j in enumerate(L)}
    for p, (c, a, b) in PARTS.items():
        done, att = state(p); n = b - a + 1
        A = {s for s in range(a, b + 1) if line[(p, s)] in act}
        fail = sum(max(0, att[s] - (1 if s in done or s in A else 0)) for s in att)
        print('%-5s seeds %5d-%5d (%5d) cap %5d  done %5d  active %5d  failed attempts %d' % (p, a, b, n, C.get(p, a - 1), len(done), len(A), fail))
elif cmd == 'submit':
    dry = '--dry' in sys.argv; C = caps(); q, act = queue(); L = jobs()
    ndis = sum(1 for l in q if l.split()[0].startswith('dis-'))
    room = min(USERCAP - len(q), DISCAP - ndis)
    st = {p: state(p) for p in PARTS}
    todo = []
    for i, (p, s) in enumerate(L, 1):
        done, att = st[p]
        if s <= min(C.get(p, 0), PARTS[p][2]) and s not in done and i not in act and att[s] < MAXATT: todo.append(i)
    msg = 'to submit %d, room %d' % (len(todo), room)
    todo = todo[:max(0, room)]
    if not todo: print(msg); sys.exit(0)
    rr = ranges(todo)
    for k in range(0, len(rr), 200):
        spec = ','.join(rr[k:k + 200])
        sc = ['sbatch', '--parsable', '--nice=100', '--array=' + spec, '--time=8:00:00', '--mem=1000MB', '--cpus-per-task=1',
              '-J', 'dis-z05x', '-o', '/dev/null', '--exclude=' + EXCL, HERE + '/run_z05x.sh']
        if dry: print(' '.join(sc)); continue
        msg += '; array %s' % sh(sc).strip()
    print(msg)
    if not dry:
        with open(Z + '/submitted.log', 'a') as f: f.write(sh(['date', '+%F %T']).strip() + ' ' + msg + '\n')
