#!/usr/bin/env python3
"""Submit (or resubmit) the lines of a jobs.list that are neither finished
("done" marker) nor pending/running in an array of the same job name, within
the queue caps. Array task i = line i (OFFSET 0), so this is also the
resubmission of tasks lost to NODE_FAIL/TIMEOUT. Safe to run repeatedly.

  submit.py NAME LIST TIME MEM [--cpus N] [--max N] [--nice 100] [--dry]

Caps (AK, 6 Oct): user total (all productions) < USERCAP, own dis-* < DISCAP.
"""
import os, subprocess, sys

USERCAP = int(os.environ.get('USERCAP', 19000))
DISCAP = int(os.environ.get('DISCAP', 11500))
MAXATT = int(os.environ.get('MAXATT', 4))
HERE = os.path.dirname(os.path.abspath(__file__))
BAD = '/ptmp/mpp/akarlber/cs-production/bad_nodes'


def sh(cmd):
    return subprocess.run(cmd, check=True, text=True, capture_output=True).stdout


def ranges(v):
    out, i = [], 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[j + 1] == v[j] + 1:
            j += 1
        out.append(str(v[i]) if i == j else '%d-%d' % (v[i], v[j]))
        i = j + 1
    return ','.join(out)


def main():
    a = sys.argv[1:]
    name, lst, tlim, mem = a[:4]
    a = a[4:]
    cpus, mx, nice, dry = 1, 10**9, 100, False
    while a:
        o = a.pop(0)
        if o == '--cpus': cpus = int(a.pop(0))
        elif o == '--max': mx = int(a.pop(0))
        elif o == '--nice': nice = int(a.pop(0))
        elif o == '--dry': dry = True
    assert name.startswith('dis-')
    user = os.environ['USER']
    q = sh(['squeue', '-u', user, '-h', '-r', '-o', '%j %K %t']).splitlines()
    active = set()
    ndis = 0
    for l in q:
        f = l.split()
        if f[0].startswith('dis-'):
            ndis += 1
        if f[0] == name and f[2] in ('PD', 'R', 'CF') and f[1].isdigit():
            active.add(int(f[1]))
    dirs = [l.strip() for l in open(lst) if l.strip()]
    def attempts(d):
        try:
            return sum(1 for _ in open(os.path.join(d, 'attempts')))
        except OSError:
            return 0
    todo = [i + 1 for i, d in enumerate(dirs)
            if i + 1 not in active and not os.path.exists(os.path.join(d, 'done'))]
    # a line that failed MAXATT times is left alone (deterministic failure?)
    stuck = [i for i in todo if attempts(dirs[i - 1]) >= MAXATT]
    if stuck:
        print('%s: %d lines with >= %d attempts left alone, e.g. %s' % (name, len(stuck), MAXATT, dirs[stuck[0] - 1]))
        todo = [i for i in todo if i not in set(stuck)]
    ndone = sum(os.path.exists(os.path.join(d, 'done')) for d in dirs)
    room = min(USERCAP - len(q), DISCAP - ndis, mx)
    msg = '%s: %d lines, %d done, %d active, %d to submit' % (name, len(dirs), ndone, len(active), len(todo))
    if not todo or room <= 0:
        print(msg + ('' if todo else '') + ('; queue at cap' if todo else ''))
        return
    todo = todo[:room]
    assert todo[-1] < 40000
    logs = os.path.join(os.path.dirname(os.path.abspath(lst)), 'logs')
    os.makedirs(logs, exist_ok=True)
    cmd = ['sbatch', '--parsable', '--nice=%d' % nice, '--array=' + ranges(todo), '--time=' + tlim,
           '--mem=' + mem, '--cpus-per-task=%d' % cpus, '-J', name,
           '--export=ALL,LIST=%s,OFFSET=0,OMP_NUM_THREADS=%d' % (os.path.abspath(lst), cpus),
           '-o', os.path.join(logs, '%x.%A_%a.out')]
    bad = []
    try:
        bad = [l.strip() for l in open(BAD) if l.strip()]
    except OSError:
        pass
    # EXTRA_EXCLUDE (e.g. the repeat-offender et nodes, for long jobs)
    bad += [x for x in os.environ.get('EXTRA_EXCLUDE', '').split(',') if x]
    if bad:
        cmd.append('--exclude=' + ','.join(bad))
    cmd.append(os.path.join(HERE, 'run_task.sh'))
    if dry:
        print(msg + ' (dry): ' + ' '.join(cmd)); return
    jid = sh(cmd).strip()
    msg += '; submitted %d as array %s' % (len(todo), jid)
    print(msg)
    with open(lst + '.submitted', 'a') as f:
        f.write(sh(['date', '+%F %T']).strip() + ' ' + msg + '\n')


main()
