#!/bin/bash -l
# Slurm array task: task i runs the job directory in line (i + OFFSET) of $LIST.
# Each job directory holds cmd.sh (the run, env included) and pattern (an
# extended regex that the run's stdout must contain for success). A finished
# directory (marker "done") is skipped, so resubmitting is always safe.
#   sbatch --array=1-200 -J dis-xxx --export=ALL,LIST=/path/jobs.list run_task.sh
#SBATCH --partition=alma
#SBATCH --ntasks=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
i=$((SLURM_ARRAY_TASK_ID + ${OFFSET:-0}))
d=$(sed -n "${i}p" "$LIST")
[ -n "$d" ] && [ -d "$d" ] || { echo "task $i: no directory in $LIST"; exit 1; }
cd "$d" || exit 1
if [ -f done ]; then echo "task $i: $d already done"; exit 0; fi
echo "$(hostname) ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date +%s)" >> attempts
echo "task $i: $d on $(hostname), $(date)"
/usr/bin/time -v bash ./cmd.sh > run.log 2> time.log
status=$?
if [ $status -eq 0 ] && grep -qE "$(cat pattern)" run.log; then
    touch done
else
    [ $status -eq 0 ] && status=97   # ran but no result line
fi
echo "task $i: exit $status, $(date)"
exit $status
