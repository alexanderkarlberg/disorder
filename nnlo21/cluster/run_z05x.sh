#!/bin/bash -l
# Z05x (UPDATE 20 Oct) array task: line SLURM_ARRAY_TASK_ID of $Z/jobs.list = 'part seed'; one file per job (file quota):
# the run's stdout+stderr and the /usr/bin/time -v report go to <part>_<seed>.out.tmp, renamed to <part>_<seed>.out when
# the run exits 0 with the right number of lines (see below). An existing .out is skipped. Command template: $Z/<part>/cmd (SEED replaced);
# attempts are logged in $Z/<part>/attempts. Same guards as run_z05.sh (a recent .out.tmp means another run is writing).
# Success: exit 0, exactly 31 LCELL lines and exactly one RESULT line (no scale cells: LABSET=zeus05 without SCALES).
#SBATCH --partition=alma
#SBATCH --ntasks=1
Z=/ptmp/mpp/akarlber/disorder-nnlo21/Z05x
read PART seed < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$Z/jobs.list")
[ -n "$PART" ] || exit 1
cd "$Z/$PART" || exit 1
o=${PART}_$seed.out
[ -f "$o" ] && exit 0
if [ -f "$o.tmp" ] && [ $(( $(date +%s) - $(stat -c %Y "$o.tmp") )) -lt 600 ]; then exit 3; fi
echo "$seed $(hostname) ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date +%s)" >> attempts
cmd=$(sed "s/SEED/$seed/g" cmd)
echo "# $PART $seed on $(hostname), job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}, $(date)" > "$o.tmp"
echo "# $cmd" >> "$o.tmp"
/usr/bin/time -v bash -c "$cmd" >> "$o.tmp" 2>&1
status=$?
if [ $status -eq 0 ] && [ "$(grep -c '^ LCELL ' "$o.tmp")" -eq 31 ] && [ "$(grep -c RESULT "$o.tmp")" -eq 1 ]; then
    mv "$o.tmp" "$o"
else
    [ $status -eq 0 ] && status=97
fi
exit $status
