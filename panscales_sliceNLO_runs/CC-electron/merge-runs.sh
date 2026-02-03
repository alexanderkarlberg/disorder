#!/bin/bash

dir=$(basename "$PWD") # CUrrent directory

FODEST=/home/karlberg/PanScales/panscales-vbf/logbook/2024-01-17-slice-NLO/paper-runs-NLO-checks/DIS/disorder_refs/$dir
PHENODEST=/home/karlberg/PanScales/panscales-vbf/logbook/2024-01-17-slice-NLO/paper-runs-pheno/DIS/disorder_refs/$dir

# Check if the destinations exist, otherwise create them
[ -d "$FODEST" ] || mkdir -p "$FODEST"
[ -d "$PHENODEST" ] || mkdir -p "$PHENODEST"

# Fixed order runs
mergedata 1 disorder_lo_seed000*
mv fort.12 ${FODEST}/lo.top

mergedata 1 disorder_nlo_seed000*
mv fort.12 ${FODEST}/nlocoeff.top

# Pheno runs
mergedata 1 phenodisorder_lo_seed000*μR_1.0_μF_1.0.dat
mv fort.12 ${PHENODEST}/lo.top

mergedata 1 phenodisorder_nlo_seed000*μR_1.0_μF_1.0.dat
mv fort.12 ${PHENODEST}/nlo.top

mergedata 1 phenodisorder_nnlo_seed000*μR_1.0_μF_1.0.dat
mv fort.12 ${PHENODEST}/nnlo.top

mergedata 1 phenodisorder_n3lo_seed000*μR_1.0_μF_1.0.dat
mv fort.12 ${PHENODEST}/n3lo.top

for mur in 1.0 2.0 0.5
do
    for muf in 1.0 2.0 0.5
    do
	f=phenodisorder_nlo_seed0001_pdfmem000_μR_${mur}_μF_${muf}.dat
	if [ ! -f "$f" ]; then
            continue
	fi
	mergedata 1 phenodisorder_nlo_seed000*_μR_${mur}_μF_${muf}.dat
	mv fort.12 nlo_μR_${mur}_μF_${muf}.top
    done
done

mergedata 4 nlo_μR_*_μF_*.top
mv fort.12 ${PHENODEST}/nlo-max.top

mergedata 5 nlo_μR_*_μF_*.top
mv fort.12 ${PHENODEST}/nlo-min.top







