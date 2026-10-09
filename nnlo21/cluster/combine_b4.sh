#!/bin/bash
# UPDATE 9 Oct: B combination with the bin4 r/rcorr sets (b0..kp unchanged: bin/ results).
# usage: combine_b4.sh [-p plain|diag|trim]; writes nnlo21/cluster/results/B_4_<mode>.{txt,json}
P=/ptmp/mpp/akarlber/disorder-nnlo21/runs; R=$(dirname "$0")/results
mode=${1:-plain}; opt=""
case $mode in diag) opt="--maxit r=1.79e5 --maxit rcorr=1.99e3";; trim) opt="--trimmed";; esac
python3 $(dirname "$0")/combine_b.py --nnlojet $R/nnlojet --json $R/B_4_$mode.json $opt \
 b0="$P/pilot/zb0/s*/run.log" b1="$P/prod/B/b1/s*/run.log" lo="$P/prod/B/lo/s*/run.log" \
 b2="$P/prod/B/b2/s*/run.log" vi="$P/prod/B/vi/s*/run.log" kp="$P/prod/B/kp/s*/run.log" \
 r="$P/prod/B/r-4/s*/run.log" rcorr="$P/prod/B/rcorr-4/s*/run.log" \
 conv="$P/prod/Bx/conv-4/s*/run.log" edge="$P/prod/Bx/edge-r-4/s*/run.log" \
 plain-r="$P/prod/Bx/plain-r-4/s*/run.log" plain-rcorr="$P/prod/Bx/plain-rcorr-4/s*/run.log" > $R/B_4_$mode.txt
