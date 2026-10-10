#!/bin/bash
# NNLOJET production pilot of a 10 Oct set: 5 seeds (9001-9005) per part, ncall as the photon ZEUS production
#   mk_nnj_pilot.sh SET      (needs the warmups; writes $P/runs/nnj-SET/prod/<part>/s9001.. and pilot.list)
P=/ptmp/mpp/akarlber/disorder-nnlo21; N=$P/runs/nnj-$1; S=/u/akarlber/work/disorder/nnlo21/cluster
sort -u $P/runs/nnlojet/prod_ncall.txt | while read part nc; do
  python3 $S/mkjobs.py $N/pilot.list $N/prod/$part 9001-9005 --pre "cd {dir}" --pattern "Elapsed time" -- $S/nnlojet_prod.sh $N $part $nc '{seed}' > /dev/null
done; wc -l $N/pilot.list
