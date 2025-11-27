#!/bin/bash


CMDLINE="../../aux/multirun.pl 8 ../disorder -Qmin 5 -pdf NNPDF40MC_nlo_as_01180 -ncall1 10000000 -ncall2 10000000000 -scale-choice 0 -noNC -CC -neutrino -positron "

echo "y" | $CMDLINE " -order-min 1 -order-max 1"

wait

echo "y" | $CMDLINE " -order-min 2 -order-max 2"
