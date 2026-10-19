#!/bin/bash -l
# Hourly feeder for Z05: z05.py submit (seeds up to the caps in $Z/caps), then schedules itself again in an hour,
# unless $Z/feeder.stop exists or the deadline has passed. Start once by hand: sbatch .../feeder_z05.sh
#SBATCH --partition=alma
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500MB
#SBATCH --time=0:20:00
#SBATCH -J dis-feeder-z05
#SBATCH -o /ptmp/mpp/akarlber/disorder-nnlo21/Z05/feeder.%j.out
S=/u/akarlber/work/disorder/nnlo21/cluster; Z=/ptmp/mpp/akarlber/disorder-nnlo21/Z05
DEADLINE=${DEADLINE:-$(date -d '2026-10-25' +%s)}
{
echo "=== $(date) on $(hostname)"
if [ -f "$Z/feeder.stop" ]; then echo "stop file, feeder stops"; exit 0; fi
if [ "$(date +%s)" -gt "$DEADLINE" ]; then echo "deadline passed, feeder stops"; exit 0; fi
if [ -z "$(squeue -u "$USER" -h -n dis-feeder-z05 -t PD -o %i)" ]; then
    sbatch --parsable -J dis-feeder-z05 --begin=now+60minutes --exclude='et[07-40]' \
        --export=ALL,DEADLINE="$DEADLINE" -o "$Z/feeder.%j.out" "$S/feeder_z05.sh"
fi
"$S/z05.py" submit
} >> "$Z/feeder.log" 2>&1
rm -f "$Z/feeder.${SLURM_JOB_ID}.out"
