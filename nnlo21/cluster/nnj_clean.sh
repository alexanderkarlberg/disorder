#!/bin/bash
# Clean a finished NNLOJET job directory (marker "done") to save files on the
# /ptmp file quota: keep done, time.log and the result histograms (*.dat, read by
# combine_nnlojet.py); pack the logs and job files into logs.tgz; delete the
# copied warmup grid (a copy of $N/warmup/<part>/s1). Only production
# directories (path contains /prod/): warmup directories hold the grids.
#   nnj_clean.sh DIR...
for d in "$@"; do
    case "$(cd "$d" && pwd)" in */prod/*) ;; *) continue;; esac
    [ -f "$d/done" ] && [ -f "$d/job.run" ] || continue
    ( cd "$d" || exit
      fs=$(ls run.log attempts cmd.sh pattern job.run DIS*.s*.log 2>/dev/null)
      tar -czf logs.tgz $fs && rm -f $fs
      for g in DIS*.y*; do case $g in *.dat) ;; *) rm -f "$g" ;; esac; done )
done
