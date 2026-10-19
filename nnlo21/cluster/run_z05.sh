#!/bin/bash -l
# Z05 (UPDATE 18 Oct) array task: seed = task id + 60000, part from $PART; one file per job (file quota):
# the run's stdout+stderr and the /usr/bin/time -v report go to <part>_<seed>.out.tmp, renamed to
# <part>_<seed>.out when the run exits 0 with all 217 SCELL lines. An existing .out is skipped.
# The command template is $Z/<part>/cmd (SEED replaced); attempts are logged in $Z/<part>/attempts.
#SBATCH --partition=alma
#SBATCH --ntasks=1
Z=/ptmp/mpp/akarlber/disorder-nnlo21/Z05
seed=$((SLURM_ARRAY_TASK_ID + 60000))
cd "$Z/$PART" || exit 1
o=${PART}_$seed.out
[ -f "$o" ] && exit 0
# another run of this seed still writing (e.g. a cancelled task on a node stuck in COMPLETING, 19 Oct): leave it
if [ -f "$o.tmp" ] && [ $(( $(date +%s) - $(stat -c %Y "$o.tmp") )) -lt 600 ]; then exit 3; fi
echo "$seed $(hostname) ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date +%s)" >> attempts
cmd=$(sed "s/SEED/$seed/g" cmd)
echo "# $PART $seed on $(hostname), job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}, $(date)" > "$o.tmp"
echo "# $cmd" >> "$o.tmp"
/usr/bin/time -v bash -c "$cmd" >> "$o.tmp" 2>&1
status=$?
if [ $status -eq 0 ] && [ "$(grep -c '^ SCELL ' "$o.tmp")" -eq 217 ] && [ "$(grep -c RESULT "$o.tmp")" -eq 1 ]; then
    mv "$o.tmp" "$o"
else
    [ $status -eq 0 ] && status=97
fi
exit $status
