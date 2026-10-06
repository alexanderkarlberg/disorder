#!/bin/bash
# NNLOJET production job: copy the warmup grid of PART into the job directory,
# write the runcard (production = NCALL[1], iseed = SEED) and run.
#   nnlojet_prod.sh RUNDIR PART NCALL SEED      (cwd = job directory)
set -e
N=$1; part=$2; nc=$3; seed=$4
ulimit -s unlimited; export OMP_STACKSIZE=1024M OMP_NUM_THREADS=1
w=$N/warmup/$part/s1
[ -f $w/done ] || { echo "warmup of $part not finished"; exit 3; }
for g in $w/DIS.*.y*; do case $g in *.txt|*.log) ;; *) cp -f $g . ;; esac; done
NNLOJET=$N/../../nnlojet/bin/NNLOJET /u/akarlber/work/disorder/nnlo21/cluster/nnlojet_card.py $N/template.run $part production $nc 1 $seed > job.run
exec $N/../../nnlojet/bin/NNLOJET --run job.run
