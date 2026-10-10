#!/bin/bash
# submit the NNLOJET pilot jobs of every part whose warmup is done (10 Oct sets), until all pilots are done
P=/ptmp/mpp/akarlber/disorder-nnlo21; S=/u/akarlber/work/disorder/nnlo21/cluster
declare -A sub
while true; do
  left=0
  for s in E Z W; do N=$P/runs/nnj-$s
    for w in $N/warmup/*; do part=$(basename $w)
      l=$N/prod/$part/pilot.list; [ -f $l ] || { grep "/prod/$part/" $N/pilot.list > $l; }
      n=$(grep -c . $l); d=$(for x in $(cat $l); do [ -f $x/done ] && echo; done | wc -l)
      [ $d -lt $n ] && left=1
      if [ -f $w/s1/done ] && [ -z "${sub[$s$part]}" ]; then
        $S/submit.py dis-nj${s}p-$part $l 6:00:00 2000MB --cpus 1 >/dev/null 2>&1; sub[$s$part]=1; fi
    done
  done
  [ $left -eq 0 ] && break
  sleep 300
done
