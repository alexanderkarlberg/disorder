#!/bin/bash
# Build the NLO 3+1 integrator nlo31 standalone (dis31 sources + LHAPDF).
#   dis31/build_nlo31.sh <build directory> [lhapdf-config]
# Needs gfortran (>= 9) and LHAPDF 6 with its Fortran interface.
set -e
D=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$1"; W=$(cd "$1" && pwd)
LC=${2:-lhapdf-config}
cd "$W"
FF="gfortran -O2 -ffixed-line-length-132 -std=legacy -I$D/mcfm/Inc"
F9="gfortran -O2 -ffree-line-length-none -I$D/mcfm/Inc"
O=""
for f in $D/mcfm/*.f $D/mcfm/loop/*.f; do
  o=$(basename $(dirname $f))_$(basename $f .f).o
  $FF -c $f -o $o; O="$O $o"
done
for m in psmc born31 me41 dip41 virt31 iop31 nlo31; do
  $F9 -c $D/$m.f90 -o $m.o; O="$O $m.o"
done
L=$($LC --libdir)
gfortran -O2 -o nlo31 $O -L$L -lLHAPDF -Wl,-rpath,$L -lstdc++
echo "built $W/nlo31"
