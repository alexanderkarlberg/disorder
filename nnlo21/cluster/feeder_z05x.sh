#!/bin/bash -l
# Hourly feeder for Z05x: z05x.py submit (seeds up to the caps in $Z/caps), then schedules itself again in an hour,
# unless $Z/feeder.stop exists or the deadline has passed. Start once by hand: sbatch .../feeder_z05x.sh
#SBATCH --partition=alma
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --time=0:20:00
#SBATCH -J dis-feeder-z05x
#SBATCH -o /ptmp/mpp/akarlber/disorder-nnlo21/Z05x/feeder.%j.out
S=/u/akarlber/work/disorder/nnlo21/cluster; Z=/ptmp/mpp/akarlber/disorder-nnlo21/Z05x
DEADLINE=${DEADLINE:-$(date -d '2026-10-25' +%s)}
{
echo "=== $(date) on $(hostname)"
if [ -f "$Z/feeder.stop" ]; then echo "stop file, feeder stops"; exit 0; fi
if [ "$(date +%s)" -gt "$DEADLINE" ]; then echo "deadline passed, feeder stops"; exit 0; fi
if [ -z "$(squeue -u "$USER" -h -n dis-feeder-z05x -t PD -o %i)" ]; then
    sbatch --parsable -J dis-feeder-z05x --begin=now+60minutes --exclude='et[07-40]' \
        --export=ALL,DEADLINE="$DEADLINE" -o "$Z/feeder.%j.out" "$S/feeder_z05x.sh"
fi
"$S/z05x.py" submit
} >> "$Z/feeder.log" 2>&1
rm -f "$Z/feeder.${SLURM_JOB_ID}.out"
