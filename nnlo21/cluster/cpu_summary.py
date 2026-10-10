#!/usr/bin/env python3
"""Mean CPU (user+sys, core-h) per job over the finished jobs of each list: cpu_summary.py DIR... (each DIR/<name>/s*/time.log)"""
import glob, os, re, sys
for base in sys.argv[1:]:
    for d in sorted(glob.glob(base + '/*/')):
        c = []
        for s in glob.glob(d + 's*/'):
            if os.path.exists(s + 'done'):
                t = open(s + 'time.log').read()
                u = float(re.search(r'User time \(seconds\): ([\d.]+)', t).group(1)) + float(re.search(r'System time \(seconds\): ([\d.]+)', t).group(1))
                c.append(u / 3600)
        if c: print('%-28s n=%2d  CPU/job %.3f core-h' % (os.path.relpath(d), len(c), sum(c) / len(c)))
