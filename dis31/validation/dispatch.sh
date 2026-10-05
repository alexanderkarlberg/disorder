#!/bin/bash
# Dispatcher v2 (5 Oct 2026, evening). Every 10 minutes: probe all machines
# IN PARALLEL and start pending jobs of the queue directory within AK's limits:
#  - thservs (thserv05..24): nice 10, at most 30 of my jobs, and the machine's
#    load (max of the 1-min load average and the running-process count) plus
#    new jobs at most 32;
#  - order: thservs first 08:00-20:00, desktops first at night (AK);
#  - th desktops (ruptime, up, not thA371a): nice 19, at most half the logical
#    cores for my jobs, never more than the free cores; skipped if another
#    user is logged in. Display-manager greeter sessions (lightdm, gdm, sddm)
#    and my own sessions do not count as users; BOINC (the institute's
#    background task, ~3 cores on idle desktops) is disregarded: its CPU is
#    subtracted from the load (AK). Skipped if any of my processes is
#    stopped (state T, overheatd).
#  - every host must see the queue directory (NFS mount of thA371a); about 50
#    desktops hang at ssh login (their automount of /home/thA371a hangs) and
#    are skipped by the probe timeout.
# A job directory (s*/) is pending without a file 'host'. One dispatcher at a
# time (flock on dispatch.lock).
# Usage: dispatch2.sh <queue dir> <command run in each job dir>
#        DRYRUN=1 dispatch2.sh ...   (one scan, print the plan, launch nothing)
Q=$1; CMD=$2
here=$(cd "$(dirname "$0")" && pwd)
LOG=$here/dispatch.log
SSH="ssh -n -o BatchMode=yes -o ConnectTimeout=5"
pending() { for d in "$Q"/s*/; do d=${d%/}; [ -f "$d/host" ] || echo "$d"; done; }
if [ -z "$DRYRUN" ]; then
  exec 9> "$here/dispatch.lock"
  flock -n 9 || { echo "$(date) another dispatcher holds the lock, exit" >> "$LOG"; exit 1; }
fi
PROBE='echo $(nproc) $(cut -d" " -f1 /proc/loadavg) $(cut -d" " -f4 /proc/loadavg | cut -d/ -f1) \
 $(ps -u akarlber -o comm= | grep -c "NNLOJET\|nlo31\|sliced21\|tau2_nlo\|nnlo11\|disorder\|dis2j_nlo") \
 $(ps -u akarlber -o stat=,comm= | grep "NNLOJET\|nlo31\|sliced21\|tau2_nlo\|nnlo11\|disorder\|dis2j_nlo" | grep -c "^T") \
 $(who | awk "{print \$1}" | grep -v -x -e akarlber -e lightdm -e gdm -e sddm | sort -u | wc -l) \
 $(ps -eo user=,comm=,pcpu= | awk "\$1 ~ /boinc/ || \$2 ~ /boinc/ {c+=\$3} END {printf \"%d\", c/100}") \
 $([ -d '"$Q"' ] && echo mounted)'
while true; do
  P=($(pending))
  [ ${#P[@]} -eq 0 ] && { echo "$(date) queue empty, exit" >> "$LOG"; exit 0; }
  # order of filling (AK): during the day (08:00-20:00) the thservs first, so
  # that nobody logs into a desktop under heavy load; at night the desktops
  # first (their CPUs are much faster than the thservs' Xeons)
  # excluded: thA366a (overheatd stops my jobs within minutes: 1 and 5 Oct),
  # thA352a (went down while my jobs ran, 1 Oct)
  desk=$(ruptime 2>/dev/null | awk '$1 ~ /^th[A-Z]/ && $2 == "up" && $1 != "thA371a" && $1 != "thA366a" && $1 != "thA352a" {print $1}')
  serv=$(for i in $(seq -w 5 24); do echo "thserv$i"; done)
  hr=$((10#$(date +%H)))
  if [ $hr -ge 8 ] && [ $hr -lt 20 ]; then hosts="$serv $desk"; else hosts="$desk $serv"; fi
  T=$(mktemp -d); for h in $hosts; do (timeout 20 $SSH $h "$PROBE" > $T/$h 2>/dev/null) & done; wait
  for h in $hosts; do
    [ ${#P[@]} -eq 0 ] && break
    set -- $(cat $T/$h 2>/dev/null)
    nc=$1; load=${2%.*}; nrun=$3; mine=$4; stopped=$5; others=$6; boinc=$7; mounted=$8
    [ -z "$nc" ] && continue
    [ "$mounted" != mounted ] && continue
    nrun=$((nrun - 1)); [ "$nrun" -gt "$load" ] && load=$nrun
    if [[ $h == thserv* ]]; then
      nice=10; free=$((32 - load - 1)); cap=$((30 - mine))
    else
      [ "$stopped" -gt 0 ] && continue
      [ "$others" -gt 0 ] && continue
      load=$((load - boinc)); [ $load -lt 0 ] && load=0
      nice=19; free=$((nc - load - 1)); cap=$((nc/2 - mine))
    fi
    n=$(( free < cap ? free : cap ))
    [ $n -le 0 ] && continue
    if [ -n "$DRYRUN" ]; then
      echo "$h: nproc $nc load $load (boinc $boinc) mine $mine -> would start $n (nice $nice)"; P=("${P[@]:$n}"); continue
    fi
    for ((k = 0; k < n && ${#P[@]} > 0; k++)); do
      d=${P[0]}; P=("${P[@]:1}")
      echo "$h" > "$d/host"
      timeout 30 ssh -f -n -o BatchMode=yes $h "cd $d && nohup nice -n $nice $CMD < /dev/null > run.log 2>&1 &" 2>> "$LOG.err"
      rc=$?
      echo "$(date '+%F %T') $h nice $nice rc $rc $d" >> "$LOG"
      [ $rc -ne 0 ] && rm -f "$d/host"
    done
  done
  rm -rf $T
  [ -n "$DRYRUN" ] && exit 0
  sleep 600
done
