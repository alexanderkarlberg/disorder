#!/usr/bin/env python3
"""Create job directories <base>/s<seed>/ with cmd.sh and pattern, and append
them to a jobs.list (one directory per line; array task i runs line i).

  mkjobs.py LIST BASE SEEDS [--env K=V]... [--pre CMD]... [--pattern RE] -- command ...

SEEDS: a-b (inclusive) or a comma list. '{seed}' in the command, --pre lines
and --env values is replaced by the seed, '{dir}' by the job directory.
An existing directory with the same cmd.sh is kept (idempotent); a different
cmd.sh is an error."""
import os, sys, shlex

def main():
    a = sys.argv[1:]
    lst, base, seeds = a[0], a[1], a[2]
    a = a[3:]
    env, pre, pat = [], [], 'RESULT'
    while a and a[0] != '--':
        o = a.pop(0)
        if o == '--env': env.append(a.pop(0))
        elif o == '--pre': pre.append(a.pop(0))
        elif o == '--pattern': pat = a.pop(0)
        else: sys.exit('unknown option ' + o)
    cmd = ' '.join(shlex.quote(x) if '{' not in x else x for x in a[1:])
    if '-' in seeds and ',' not in seeds:
        lo, hi = map(int, seeds.split('-')); sl = range(lo, hi + 1)
    else:
        sl = [int(s) for s in seeds.split(',')]
    have = set(l.strip() for l in open(lst)) if os.path.exists(lst) else set()
    new = []
    for s in sl:
        d = os.path.abspath(os.path.join(base, 's%d' % s))
        sub = lambda t: t.replace('{seed}', str(s)).replace('{dir}', d)
        txt = '#!/bin/bash\nset -e\n' + ''.join('export %s\n' % sub(e) for e in env) \
              + ''.join(sub(p) + '\n' for p in pre) + 'exec ' + sub(cmd) + '\n'
        os.makedirs(d, exist_ok=True)
        f = os.path.join(d, 'cmd.sh')
        if os.path.exists(f):
            if open(f).read() != txt:
                sys.exit('different cmd.sh exists in ' + d)
        else:
            open(f, 'w').write(txt)
            open(os.path.join(d, 'pattern'), 'w').write(pat + '\n')
        if d not in have:
            new.append(d)
    with open(lst, 'a') as fh:
        for d in new:
            fh.write(d + '\n')
    print('%s: %d directories, %d new lines' % (lst, len(sl), len(new)))

main()
