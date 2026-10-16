#!/bin/bash
# UPDATE 9 Oct: B combination with the bin4 r/rcorr sets (b0..kp unchanged: bin/ results).
# usage: combine_b4.sh [-p plain|diag|trim]; writes nnlo21/cluster/results/B_4_<mode>.{txt,json}
# BPROD=<dir>: parent of prod/B, prod/Bx (default $P; the B runs were archived 11 Oct, extract their run.log/done there).
# KP11=1 (KREDO): kp from $P/prod/B/kp11 (bin11, never archived), output in results/kp11/.
P=/ptmp/mpp/akarlber/disorder-nnlo21/runs; R=$(dirname "$0")/results; R0=$R; Q=${BPROD:-$P}; KP="$Q/prod/B/kp/s*/run.log"
if [ "$KP11" = 1 ]; then R=$R/kp11; KP="$P/prod/B/kp11/s*/run.log"; mkdir -p $R; fi
mode=${1:-plain}; opt=""
case $mode in diag) opt="--maxit r=1.79e5 --maxit rcorr=1.99e3";; trim) opt="--trimmed";; esac
python3 $(dirname "$0")/combine_b.py --nnlojet $R0/nnlojet --json $R/B_4_$mode.json $opt \
 b0="$P/pilot/zb0/s*/run.log" b1="$Q/prod/B/b1/s*/run.log" lo="$Q/prod/B/lo/s*/run.log" \
 b2="$Q/prod/B/b2/s*/run.log" vi="$Q/prod/B/vi/s*/run.log" kp="$KP" \
 r="$Q/prod/B/r-4/s*/run.log" rcorr="$Q/prod/B/rcorr-4/s*/run.log" \
 conv="$Q/prod/Bx/conv-4/s*/run.log" edge="$Q/prod/Bx/edge-r-4/s*/run.log" \
 plain-r="$Q/prod/Bx/plain-r-4/s*/run.log" plain-rcorr="$Q/prod/Bx/plain-rcorr-4/s*/run.log" > $R/B_4_$mode.txt
