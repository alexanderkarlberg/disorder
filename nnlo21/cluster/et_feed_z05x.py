#!/usr/bin/env python3
"""Z05x, 21 Oct: keep SLOTS tasks (running + pending) of dis-z05x-<node> on each node listed in $Z/etnodes, while the main
arrays sit pending behind other users. More than ~100-120 of the tasks starting on one et node at once never start their
script (hang), so every node gets at most SLOTS. Single-node arrays by excluding every other node (a --nodelist would ask
for all listed nodes). Lines of jobs.list are taken from the top (highest line first), skipping done seeds, vixs (cap 0), seeds
above the caps in $Z/caps and lines already running or pending in an et array. Tasks running > 20 min with no .out.tmp are
cancelled by job id (hung), logged in $Z/etfeed.log; a node with > 30 % hung in a round of >= 20 is written to $Z/etnodes.off
and gets no more tasks. Stop file: $Z/etfeed.stop.
  et_feed_z05x.py [--dry]"""
import os, re, subprocess, sys, time
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; Z = P + '/Z05x'
HERE = os.path.dirname(os.path.abspath(__file__))
SLOTS = int(os.environ.get('SLOTS', 100))
def sh(c): return subprocess.run(c, check=True, text=True, capture_output=True).stdout
def readlist(f):
    return [l.split()[0] for l in open(f) if l.strip() and not l.startswith('#')] if os.path.exists(f) else []
def ranges(v):
    v = sorted(v); out, i = [], 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[j + 1] == v[j] + 1: j += 1
        out.append(str(v[i]) if i == j else '%d-%d' % (v[i], v[j])); i = j + 1
    return out
def log(m):
    with open(Z + '/etfeed.log', 'a') as f: f.write(time.strftime('%F %T ') + m + '\n')
dry = '--dry' in sys.argv
if os.path.exists(Z + '/etfeed.stop'): sys.exit(0)
nodes = [n for n in readlist(Z + '/etnodes') if n not in readlist(Z + '/etnodes.off')]
L = [l.split() for l in open(Z + '/jobs.list')]
caps = {l.split()[0]: int(l.split()[1]) for l in open(Z + '/caps') if l.strip() and not l.startswith('#')}
q = sh(['squeue', '-u', os.environ['USER'], '-h', '-r', '-o', '%j %K %t %N %M']).splitlines()
busy = {}      # node -> number of tasks running or pending in its array
inet = set()   # lines running or pending in any et array
hung = {}      # node -> list of job ids; running count per node
run = {}
for l in q:
    f = l.split()
    if len(f) < 4 or not f[0].startswith('dis-z05x') or not f[1].isdigit(): continue
    t = int(f[1])
    # a running task counts for the node it runs on (hand-made single-node arrays too); a pending one only in an et array
    n = f[3].split(',')[0] if f[2] in ('R', 'CF') else (f[0][9:] if f[0].startswith('dis-z05x-') else None)
    if n is None or not n.startswith('et'): continue
    if f[2] in ('R', 'PD', 'CF'): busy[n] = busy.get(n, 0) + 1; inet.add(t)
    if f[2] == 'R' and len(f) >= 5:
        run[n] = run.get(n, 0) + 1
        tm = re.split('[:-]', f[4])
        mins = int(tm[-2]) + (60 * int(tm[-3]) if len(tm) >= 3 else 0) + (1440 * int(tm[-4]) if len(tm) >= 4 else 0)
        p, s = L[t - 1]
        if mins > 20 and not os.path.exists('%s/%s/%s_%s.out.tmp' % (Z, p, p, s)):
            jid = sh(['squeue', '-u', os.environ['USER'], '-h', '-r', '-n', f[0], '-t', 'R', '-O', 'JobID:30,ArrayTaskID:10']).split('\n')
            ids = [x.split()[0] for x in jid if x.split() and x.split()[1] == str(t)]
            hung.setdefault(n, []).extend(ids)
for n, ids in hung.items():
    log('%s: %d hung (of %d running) cancelled: %s' % (n, len(ids), run.get(n, 0), ' '.join(ids[:5]) + (' ...' if len(ids) > 5 else '')))
    if not dry: subprocess.run(['scancel'] + ids)
    if len(ids) >= 20 and len(ids) > 0.3 * run.get(n, 1):
        log('%s: > 30 %% hung, node switched off' % n)
        if not dry: open(Z + '/etnodes.off', 'a').write(n + '\n')
done = {}
for p in {p for p, s in L}:
    done[p] = {int(m.group(1)) for f in os.listdir(Z + '/' + p) for m in [re.match(r'%s_(\d+)\.out$' % p, f)] if m}
cand = [i + 1 for i, (p, s) in enumerate(L)
        if p != 'vixs' and int(s) <= caps.get(p, 0) and int(s) not in done[p] and (i + 1) not in inet]
cand.sort(reverse=True)
allnodes = sorted(set(sh(['sinfo', '-h', '-N', '-p', 'alma', '-o', '%N']).split()))
for n in nodes:
    room = SLOTS - busy.get(n, 0)
    if room < 10 or not cand: continue
    take, cand = cand[:room], cand[room:]
    ex = ','.join(x for x in allnodes if x != n)
    for k in range(0, len(ranges(take)), 200):
        spec = ','.join(ranges(take)[k:k + 200])
        sc = ['sbatch', '--parsable', '--nice=0', '--array=' + spec, '--time=8:00:00', '--mem=1000MB', '--cpus-per-task=1',
              '-J', 'dis-z05x-' + n, '-o', '/dev/null', '--exclude=' + ex, HERE + '/run_z05x.sh']
        if dry: print(n, len(take), spec[:60]); continue
        log('%s: +%d tasks (array %s), %d running %d busy before' % (n, len(take), sh(sc).strip(), run.get(n, 0), busy.get(n, 0)))
