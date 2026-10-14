#!/usr/bin/env python3
"""Completeness check of the J7/N7/Ex7 lists (14 Oct): per list, dirs with 'done', outputs present and finite.
  check7.py   (reads the groups files below; ours: run.log must hold ZCELL/CELL rows, all finite; NNLOJET: the 4 DIS*.dat finite, 1 file each)"""
import glob, os, re, sys
import numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'
G = ['ours7j/PH/groups.txt', 'ours7j/W/groups.txt', 'nnj-groups7-RRa1.txt', 'nnj-groups7-rest.txt', 'nnj-groups7-Z.txt', 'ours7/groups7.txt']
def ff(x):
    try: return float(x)
    except ValueError: return float(re.sub(r"(\d)([+-]\d{3})$", r"\1E\2", x))
tot = {}
for g in G:
    for l in open(P + '/' + g):
        if not l.strip() or l.startswith('#'): continue
        name, lst = l.split()[:2]; dirs = [d.strip() for d in open(lst) if d.strip()]
        nd = nodone = bad = 0; badl = []
        for d in dirs:
            if not os.path.exists(d + '/done'): nodone += 1; badl.append((d, 'nodone')); continue
            nd += 1
            try:
                if 'nnj-' in d:
                    fs = glob.glob(d + '/DIS*.dat'); 
                    ok = len(fs) == 4
                    for f in fs:
                        a = np.array([[ff(v) for v in x.split()] for x in open(f) if x.strip() and not x.startswith('#')])
                        ok = ok and a.size > 0 and np.isfinite(a).all()
                else:
                    r = [[ff(v) for v in x.split()[3:]] for x in open(d + '/run.log') if x.startswith((' CELL ', ' ZCELL ')) ]
                    r2 = [[ff(v) for v in x.split()[4:]] for x in open(d + '/run.log') if x.startswith(' ZCELL ')]
                    r = r or r2
                    ok = len(r) > 0 and np.isfinite(np.array(r)).all()
            except Exception as e: ok = False
            if not ok: bad += 1; badl.append((d, 'bad'))
        print('%-16s list %5d done %5d notdone %3d done-but-bad %3d' % (name, len(dirs), nd, nodone, bad), badl[:5])
