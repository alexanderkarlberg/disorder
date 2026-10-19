#!/usr/bin/env python3
# Z05: CPU per job (user+sys of /usr/bin/time -v, core-h) per part and node type, from the finished .out files.
import glob, re, sys
from collections import defaultdict
Z = '/ptmp/mpp/akarlber/disorder-nnlo21/Z05'
tot = 0
print('%-4s %6s %7s %7s %8s  %s' % ('part', 'n', 'mean', 'max', 'core-h', 'by node type: mean (n)'))
for p in ['r', 'vi', 'kp', 'b2', 'rc1', 'rc2', 'b0', 'b1', 'lo']:
    by = defaultdict(list); allc = []
    for f in glob.glob('%s/%s/%s_*.out' % (Z, p, p)):
        host = None; u = s = None
        for l in open(f, errors='replace'):
            if host is None and l.startswith('# %s ' % p):
                m = re.search(r' on ([a-z]+)\d', l); host = m.group(1) if m else '?'
            elif 'User time (seconds)' in l: u = float(l.split()[-1])
            elif 'System time (seconds)' in l: s = float(l.split()[-1])
        if u is None: continue
        c = (u + (s or 0))/3600; by[host or '?'].append(c); allc.append(c)
    tot += sum(allc)
    print('%-4s %6d %7.3f %7.2f %8.0f  %s' % (p, len(allc), sum(allc)/len(allc), max(allc), sum(allc),
          '  '.join('%s %.3f (%d)' % (h, sum(v)/len(v), len(v)) for h, v in sorted(by.items()))))
print('total %.0f core-h' % tot)
