#!/bin/bash


CMDLINE="../../aux/multirun.pl 8 ../disorder -Qmin 5 -pdf NNPDF40MC_nlo_as_01180 -ncall1 10000000 -ncall2 10000000 -CC -noNC -neutrino -positron -prefix pheno -scaleuncert "

echo "y" | $CMDLINE " -lo"
wait
echo "y" | $CMDLINE " -nlo"
wait
echo "y" | $CMDLINE " -nnlo"
wait
echo "y" | $CMDLINE " -n3lo"
