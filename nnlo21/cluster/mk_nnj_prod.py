#!/usr/bin/env python3
"""NNLOJET production job lists of a 10 Oct set (E, Z, W) from the pilot sizing (11 Oct). NOT submitting.
  mk_nnj_prod.py          -> writes $P/runs/nnj-SET/prod/<part>/s1001.. (mkjobs), prod.list per set and
                             $P/runs/nnj-groups.txt (one line per set and part, for feeder.sh); prints the cost table.
Sizing: size_nnlojet.py with the per-set settings below; the finished pilot seeds (9001-9005) count as production jobs."""
import os, re, subprocess, glob
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'; S = '/u/akarlber/work/disorder/nnlo21/cluster'
CFG = {  # set: (env for size_nnlojet.py)
 'E': dict(OBS='dis_tauzq', NNLO_TOT_FB='5500', NNLO_REF_PB='200', FLOORFRAC='0.5'),
 'Z': dict(OBS='q2,ptavg_12,m12', NNLO_TOT_FB='800', NNLO_REF_PB='40', FLOORFRAC='0.4'),
 'W': dict(OBS='q2,ptavg_12,m12', NNLO_TOT_FB='30', NNLO_REF_PB='1.5', FLOORFRAC='0.4', FALLBACK='RRb_2=RRb_1'),
}
ncall = {}
for l in open(P + '/nnlojet/prod_ncall.txt'):
    a = l.split()
    if len(a) == 2: ncall[a[0]] = a[1]
groups = []; summary = {}
for s, cfg in CFG.items():
    N = '%s/nnj-%s' % (P, s)
    env = dict(os.environ, NLO_FLOOR='0.003', MINJOBS='2'); env.update(cfg)
    out = subprocess.run(['python3', S + '/size_nnlojet.py', N + '/prod'], env=env, capture_output=True, text=True).stdout
    open(N + '/sizing.txt', 'w').write(out)
    lst = N + '/prod.list'; open(lst, 'w').close()
    cost = 0; njobs = 0
    for l in out.splitlines():
        m = re.match(r'\s+(\S+)\s+sigma/job\s+\S+ fb\s+CPU/job\s+(\d+) s\s+->\s+(\d+) jobs', l)
        if not m: continue
        part, cpu, n = m.group(1), int(m.group(2)), int(m.group(3))
        done = len([d for d in glob.glob('%s/prod/%s/s9*' % (N, part)) if os.path.exists(d + '/done')])
        nprod = max(n - done, 0)
        if nprod == 0: continue
        subprocess.run(['python3', S + '/mkjobs.py', lst, '%s/prod/%s' % (N, part), '1001-%d' % (1000 + nprod), '--pre', 'cd {dir}',
                        '--pattern', 'Elapsed time', '--', S + '/nnlojet_prod.sh', N, part, ncall[part], '{seed}'],
                       check=True, stdout=subprocess.DEVNULL)
        pl = '%s/prod/%s/prod.list' % (N, part)
        open(pl, 'w').write(''.join('%s/prod/%s/s%d\n' % (N, part, 1000 + i) for i in range(1, nprod + 1)))
        h = 12 if cpu > 9000 else 6
        groups.append('dis-nj%s-%s %s %d:00:00 2000MB 1 et[07-40]' % (s, part, pl, h))
        cost += (nprod + done) * cpu / 3600; njobs += nprod
    summary[s] = (njobs, cost)
open(P + '/nnj-groups.txt', 'w').write('# NAME LIST TIME MEM CPUS EXCLUDE  NNLOJET production E, Z, W (11 Oct)\n' + '\n'.join(groups) + '\n')
for s, (n, c) in summary.items(): print('%s: %d production jobs, %.0f core-h (incl. done pilots)' % (s, n, c))
print('total %.0f core-h' % sum(c for n, c in summary.values()))
