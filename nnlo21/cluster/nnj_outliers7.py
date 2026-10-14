#!/usr/bin/env python3
"""N7: per part, the most extreme seed of the NNLOJET cross value (|x - median|/1.48 MAD), and the mean with that seed removed.  nnj_outliers7.py W|Z -> results/N7_SET_outliers.txt"""
import glob, os, sys, numpy as np
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs/nnj-%s/prod' % sys.argv[1]; out = []
for pd in sorted(glob.glob(P + '/*')):
    v, s = [], []
    for d in glob.glob(pd + '/s*/'):
        if not os.path.exists(d + 'done'): continue
        f = glob.glob(d + 'DIS*.cross.s*.dat')
        if len(f) != 1: continue
        l = [x for x in open(f[0]) if not x.startswith('#') and x.strip()][0].split(); v.append(float(l[0])); s.append(int(d.rstrip('/').split('/s')[-1]))
    v = np.array(v); s = np.array(s)
    if len(v) < 20: continue
    md = np.median(v); mad = 1.4826*np.median(np.abs(v - md)); z = np.abs(v - md)/mad; i = np.argmax(z); m = np.ones(len(v), bool); m[i] = False
    out.append('%-6s N=%5d mean %10.5g +- %-8.3g median %10.5g 1.48MAD %9.4g  max z %8.1f (seed %d, value %.5g fb)  mean without it %10.5g +- %-8.3g  n(z>20) %d' % (
        os.path.basename(pd), len(v), v.mean(), v.std(ddof=1)/len(v)**.5, md, mad, z[i], s[i], v[i], v[m].mean(), v[m].std(ddof=1)/(len(v) - 1)**.5, (z > 20).sum()))
open(os.path.dirname(os.path.abspath(__file__)) + '/results/N7_%s_outliers.txt' % sys.argv[1], 'w').write('cross value [fb]\n' + '\n'.join(out) + '\n')
print('\n'.join(out))
