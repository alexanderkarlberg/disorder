#!/bin/bash
# Freeze an NNLOJET warmup so that its production can start:
#   freeze_warmup.sh WARMUPDIR JOBID_TASK MINITS
# If WARMUPDIR/done is missing and either >= MINITS iterations have written
# their grid or the warmup job has ended (timeout etc.) with >= 1 grid, cancel
# that one warmup task (own job, by id), wait until it has left the queue (so
# the grid file is not rewritten while productions copy it), and write done
# plus a note. NNLOJET rewrites the grid after every iteration.
w=$1; jt=$2; minits=$3
[ -f "$w/done" ] && exit 0
n=$(grep -c 'Writing grid' "$w/run.log" 2>/dev/null); n=${n:-0}
running=$(squeue -h -j "$jt" -o %T 2>/dev/null)
if [ "$n" -ge "$minits" ] || { [ -z "$running" ] && [ "$n" -ge 1 ]; }; then
    [ -n "$running" ] && scancel "$jt"
    for i in $(seq 60); do [ -z "$(squeue -h -j "$jt" 2>/dev/null)" ] && break; sleep 5; done
    [ -n "$(squeue -h -j "$jt" 2>/dev/null)" ] && { echo "$w: $jt still in queue, not frozen"; exit 1; }
    n=$(grep -c 'Writing grid' "$w/run.log")
    echo "frozen $(date) after $n iterations (job $jt), freeze_warmup.sh" > "$w/frozen"
    touch "$w/done"; echo "$w: frozen after $n iterations"
else
    echo "$w: $n iterations, job ${running:-gone}"
fi
