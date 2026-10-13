#!/usr/bin/env python3
"""13 Oct round N7 (UPDATE 13 Oct (b)): top-up of the NNLOJET W (ZEUS dijets) production nnj-W/prod (= nnj5/W of 10 Oct).
  mk_n7.py   -> new seeds 11001- in $P/runs/nnj-W/prod/<part>/ (same directories, so combine_nnlojet.py picks them up),
                job list prod/<part>/prod7.list, groups $P/runs/nnj-groups7-RRa1.txt (RRa_1 only) and nnj-groups7-rest.txt.
Targets (total seeds) from the instructions; existing seeds are 1001- and 9001-9005 (finished pilots)."""
import os, subprocess, sys
P = '/ptmp/mpp/akarlber/disorder-nnlo21/runs'; S = '/u/akarlber/work/disorder/nnlo21/cluster'; SET = sys.argv[1] if len(sys.argv) > 1 else 'W'   # W (default) or Z (item 5, 13 Oct: RRa_3 x6, RRa_2 and RRb_3 x2)
N = P + '/nnj-' + SET
TARGET = {'RRa_1': 5500, 'RRb_1': 5100, 'RRa_2': 1850, 'RRa_3': 420, 'RRa_5': 550, 'RRb_5': 940} if SET == 'W' else {'RRa_3': 721*6, 'RRa_2': 131*2, 'RRb_3': 882*2}
HOURS = {'RRa_3': 12} if SET == 'Z' else {}
SEED0 = 11001
ncall = dict(l.split() for l in open(P + '/nnlojet/prod_ncall.txt') if len(l.split()) == 2)
grp = {}
for part, tot in TARGET.items():
    have = len([d for d in os.listdir('%s/prod/%s' % (N, part)) if d.startswith('s') and os.path.exists('%s/prod/%s/%s/done' % (N, part, d))])
    new = tot - have
    lst = '%s/prod/%s/prod7.list' % (N, part)
    subprocess.run(['python3', S + '/mkjobs.py', lst, '%s/prod/%s' % (N, part), '%d-%d' % (SEED0, SEED0 + new - 1), '--pre', 'cd {dir}',
                    '--pattern', 'Elapsed time', '--', S + '/nnlojet_prod.sh', N, part, ncall[part], '{seed}'], check=True, stdout=subprocess.DEVNULL)
    print(part, 'existing done', have, 'new', new, 'seeds %d-%d' % (SEED0, SEED0 + new - 1))
    grp[part] = 'dis-n7%s-%s %s %d:00:00 2000MB 1 et[07-40]' % ('z' if SET == 'Z' else '', part, lst, HOURS.get(part, 6))
if SET == 'Z': open(P + '/nnj-groups7-Z.txt', 'w').write('\n'.join(grp.values()) + '\n')
else:
    open(P + '/nnj-groups7-RRa1.txt', 'w').write(grp['RRa_1'] + '\n')
    open(P + '/nnj-groups7-rest.txt', 'w').write('\n'.join(v for k, v in grp.items() if k != 'RRa_1') + '\n')
