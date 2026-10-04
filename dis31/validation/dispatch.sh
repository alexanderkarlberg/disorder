#!/bin/bash
# Dispatcher for the NNLOJET RR reference jobs (dis31 validation, weekend
# 3-5 Oct 2026). Every 10 minutes: scan the machines and start pending jobs
# of the queue directory within AK's limits:
#  - thservs (thserv05..24): nice 10, at most 30 of my jobs, and the machine's
#    load plus new jobs at most 32 (hyperthreading slows jobs beyond that);
#  - th desktops (from ruptime): nice 19, at most half the logical cores for
#    my jobs, and never more than the free cores (nproc - load - 1); a host
#    with any of my processes stopped (state T, overheatd) gets no new jobs.
# A job directory is pending without a file 'host', running/started with it,
# done when run.log contains 'Elapsed time'. Hosts without the NFS mount of
# the queue directory are skipped (ssh -f returns 0 even if the cd fails).
# Usage: dispatch.sh <queue dir> <command run in each job dir>
Q=$1; CMD=$2
LOG=$(dirname "$0")/dispatch.log
SSH="ssh -n -o BatchMode=yes -o ConnectTimeout=5"
# only job directories count (a binary named s* in the queue directory made
# the queue never empty, 3-4 Oct); one dispatcher at a time (two started
# together overfilled thserv05/06 on 4 Oct): lock
pending() { for d in "$Q"/s*/; do d=${d%/}; [ -f "$d/host" ] || echo "$d"; done; }
exec 9> "$(dirname "$0")/dispatch.lock"
flock -n 9 || { echo "$(date) another dispatcher holds the lock, exit" >> "$(dirname "$0")/dispatch.log"; exit 1; }
while true; do
  P=($(pending))
  [ ${#P[@]} -eq 0 ] && { echo "$(date) queue empty, exit" >> "$LOG"; exit 0; }
  hosts=$(for i in $(seq -w 5 24); do echo "thserv$i"; done; ruptime 2>/dev/null | awk '$1 ~ /^th[A-Z]/ && $2 == "up" && $1 != "thA371a" {print $1}')
  for h in $hosts; do
    [ ${#P[@]} -eq 0 ] && break
    info=$(timeout 15 $SSH $h 'echo $(nproc) $(cut -d" " -f1 /proc/loadavg) $(ps -u akarlber -o comm= | grep -c "NNLOJET\|nlo31") $(ps -u akarlber -o stat=,comm= | grep "NNLOJET\|nlo31" | grep -c "^T") $([ -d '"$Q"' ] && echo mounted)' 2>/dev/null) || continue
    set -- $info; nc=$1; load=${2%.*}; mine=$3; stopped=$4; mounted=$5
    [ -z "$nc" ] && continue
    # the job directory must be visible on the host (NFS mount of thA371a)
    [ "$mounted" != mounted ] && continue
    if [[ $h == thserv* ]]; then
      nice=10; free=$((32 - load - 1)); cap=$((30 - mine))
    else
      [ "$stopped" -gt 0 ] && continue
      nice=19; free=$((nc - load - 1)); cap=$((nc/2 - mine))
    fi
    n=$(( free < cap ? free : cap ))
    [ $n -le 0 ] && continue
    for ((k = 0; k < n && ${#P[@]} > 0; k++)); do
      d=${P[0]}; P=("${P[@]:1}")
      echo "$h" > "$d/host"
      timeout 30 ssh -f -n -o BatchMode=yes $h "cd $d && nohup nice -n $nice $CMD < /dev/null > run.log 2>&1 &" 2>> "$LOG.err"
      rc=$?
      hf=$([ -f "$d/host" ] && echo ok || echo MISSING)
      echo "$(date '+%F %T') $h nice $nice rc $rc hostfile $hf $d" >> "$LOG"
      [ $rc -ne 0 ] && rm -f "$d/host"
    done
  done
  sleep 600
done
