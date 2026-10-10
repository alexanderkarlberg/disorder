#!/bin/bash
# Set up the NNLOJET warmups of a validation set (10 Oct).
#   mk_nnj_set.sh SET TEMPLATE.run NNLOJET_INSTALL_DIR
# Creates $P/runs/nnj-SET/{template.run,nnlojet -> install,warmup/<part>/s1,warm-{big,small}.list}.
# Warmups: LO 1M[5]; R,V 2M[5]; VV 1M[5] (8 threads); RV, RRb 10M[6], RRa 8M[4] (5M[4] for Z RRa_2-5, E, W; cheaper) (32 threads).
set -e
SET=$1; T=$2; I=$3
P=/ptmp/mpp/akarlber/disorder-nnlo21; N=$P/runs/nnj-$SET; S=/u/akarlber/work/disorder/nnlo21/cluster
mkdir -p $N; cp $T $N/template.run
rm -f $N/nnlojet; ln -s $I $N/nnlojet
: > $N/warm-big.list; : > $N/warm-small.list
w() { # part ncall niter list
  d=$N/warmup/$1/s1; mkdir -p $d
  printf '#!/bin/bash\nset -e\nexport OMP_STACKSIZE=1024M\nulimit -s unlimited\ncd %s\nNNLOJET=%s/nnlojet/bin/NNLOJET %s/nnlojet_card.py %s/template.run %s warmup %s %s 1 > job.run\nexec %s/nnlojet/bin/NNLOJET --run job.run\n' $d $N $S $N $1 $2 $3 $N > $d/cmd.sh
  echo "Elapsed time" > $d/pattern; echo $d >> $N/$4
}
w LO 1000000 5 warm-small.list; w R 2000000 5 warm-small.list; w V 2000000 5 warm-small.list
for k in 1 2 3 4 5 6 7; do w VV_$k 1000000 5 warm-small.list; w RV_$k 10000000 6 warm-big.list; done
for k in 1 2 3 4 5; do w RRb_$k 10000000 6 warm-big.list; w RRa_$k 8000000 4 warm-big.list; done
wc -l $N/warm-*.list
