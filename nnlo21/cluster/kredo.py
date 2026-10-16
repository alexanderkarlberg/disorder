#!/usr/bin/env python3
"""KREDO (UPDATE 16 Oct (b)): rerun kp / kp-x jobs with the fixed nlo31 (bin11), same commands, environment and seeds.
  kredo.py NAME OLDDIR [OLDDIR ...]  -> for each OLDDIR (a part directory with s<seed>/cmd.sh) a sibling OLDDIR + '11'
  with s<seed>/{cmd.sh,pattern} for every *done* old seed (the seeds that entered the reported results); cmd.sh is the old
  one with only the nlo31 path replaced by $P/bin11/nlo31 (checked: exactly one nlo31 path, of $P/bin*/). Appends
  'NAME-<i> <dir>/jobs.list 2:00:00 1000MB 1 EXCL' lines to $P/runs/groups-kredo.txt (idempotent)."""
import glob, os, re, sys
P = '/ptmp/mpp/akarlber/disorder-nnlo21'; NEW = P + '/bin11/nlo31'
EXCL = 'et[07-40],gt25,gt27,gt28,gt29,gt39,gt40,kt01,kt03,kt10'
G = P + '/runs/groups-kredo.txt'
name = sys.argv[1]; lines = open(G).read().split('\n') if os.path.exists(G) else []
for k, old in enumerate(sys.argv[2:]):
    old = old.rstrip('/'); new = old + '11'; os.makedirs(new, exist_ok=True); seeds = []
    for d in sorted(glob.glob(old + '/s*/'), key=lambda s: int(re.search(r's(\d+)/$', s).group(1))):
        if not os.path.exists(d + 'done'): continue
        c = open(d + 'cmd.sh').read(); m = re.findall(re.escape(P) + r'/bin\w*/nlo31', c)
        assert len(m) == 1, (d, m)
        c2 = c.replace(m[0], NEW); s = os.path.basename(d.rstrip('/')); nd = '%s/%s/' % (new, s)
        if os.path.exists(nd + 'cmd.sh'): assert open(nd + 'cmd.sh').read() == c2, nd
        else:
            os.makedirs(nd, exist_ok=True); open(nd + 'cmd.sh', 'w').write(c2); os.chmod(nd + 'cmd.sh', 0o755)
            pat = d + 'pattern'
            if os.path.exists(pat): open(nd + 'pattern', 'w').write(open(pat).read())
        seeds.append(nd.rstrip('/'))
    open(new + '/jobs.list', 'w').write('\n'.join(seeds) + '\n')
    g = 'dis-%s-%d %s/jobs.list 2:00:00 1000MB 1 %s' % (name, k, new, EXCL)
    if not any(l.split(' ')[1:2] == [new + '/jobs.list'] for l in lines if l): lines.append(g)
    print(new, len(seeds), 'jobs; old binary', m[0] if seeds else '-')
open(G, 'w').write('\n'.join(l for l in lines if l) + '\n')
