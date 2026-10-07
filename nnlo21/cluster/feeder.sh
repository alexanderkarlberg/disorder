#!/bin/bash -l
# Hourly feeder for the disorder nnlo21 production: for every group in
# $GROUPFILE (lines "NAME LIST TIME MEM [CPUS [EXCLUDE|- [NEEDS]]]"; # comments;
# NEEDS = a file that must exist first, e.g. a warmup's done marker) it submits the lines
# that are neither done nor queued (never-started lines and tasks lost to
# NODE_FAIL/TIMEOUT), within the queue caps of submit.py. Then it schedules
# itself again in an hour, unless every list is finished, the file
# $GROUPFILE.stop exists, or the deadline has passed. Run it once by hand to
# start; it runs as a tiny Slurm job (dis-feeder) afterwards.
#   GROUPFILE=/ptmp/.../groups.txt nnlo21/cluster/feeder.sh
#SBATCH --partition=alma
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --time=0:20:00
S=/u/akarlber/work/disorder/nnlo21/cluster
: "${GROUPFILE:?groups file}"
DEADLINE=${DEADLINE:-$(date -d '2026-10-14' +%s)}
log=$(dirname "$GROUPFILE")/feeder.log
mkdir -p "$(dirname "$GROUPFILE")/logs-feeder"
{
echo "=== $(date) on $(hostname)"
alldone=1
while read -r name list rest; do
    case "$name" in ''|\#*) continue;; esac
    while read -r dir; do [ -f "$dir/done" ] || { alldone=0; break; }; done < "$list"
done < "$GROUPFILE"
if [ $alldone -eq 1 ]; then echo "all lists done, feeder stops"; exit 0; fi
# schedule the next feeder first (so a node failure here does not break the
# chain); only one pending feeder at a time; never on the flaky et nodes
if [ -f "$GROUPFILE.stop" ]; then echo "stop file, feeder stops"; exit 0; fi
if [ "$(date +%s)" -gt "$DEADLINE" ]; then echo "deadline passed, feeder stops"; exit 0; fi
if [ -z "$(squeue -u "$USER" -h -n dis-feeder -t PD -o %i)" ]; then
    sbatch --parsable -J dis-feeder --begin=now+60minutes --exclude='et[07-40]' \
        --export=ALL,GROUPFILE="$GROUPFILE",DEADLINE="$DEADLINE" \
        -o "$(dirname "$GROUPFILE")/logs-feeder/feeder.%j.out" "$S/feeder.sh"
fi
while read -r name list tlim mem cpus excl needs; do
    case "$name" in ''|\#*) continue;; esac
    [ "${excl:-}" = - ] && excl=
    if [ -n "${needs:-}" ] && [ "$needs" != - ] && [ ! -e "$needs" ]; then echo "$name: waiting for $needs"; continue; fi
    EXTRA_EXCLUDE="${excl:-}" "$S/submit.py" "$name" "$list" "$tlim" "$mem" --cpus "${cpus:-1}"
done < "$GROUPFILE"
} >> "$log" 2>&1
