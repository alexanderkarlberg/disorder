#!/bin/bash
# Build nnlo21/tests/harness_hard21 (two-loop DIS 2+1 hard function tests)
# and harness_born21 (2+1 Born and one-loop helicity classes with photon + Z
# against DISENT's MATTHR/VIRTHR; run with disorder's flags, e.g.
# ./harness_born21 -toyQ0 2 -includeZ -positron).
#   nnlo21/build_harness.sh <disorder build dir> <dir of libnnlojet_core.so> [work dir]
# The disorder build must be configured with -DDISORDER_SLICING=ON and have
# built tau2_nlo (its link line is reused for DISENT and the slicing modules).
# libnnlojet_core.so: NNLOJET v1.0.2 (GPL-3.0-or-later), for the helicity
# coefficients of the 16 kinematic regions (makecoef, makecoeftaylorU).
set -e
D=$(cd "$(dirname "$0")" && pwd); R=$(dirname "$D")
B=$(cd "$1" && pwd); L=$(cd "$2" && pwd)
mkdir -p "${3:-$B/nnlo21}"; W=$(cd "${3:-$B/nnlo21}" && pwd)
cd "$W"
FF="gfortran -O2 -ffixed-line-length-132 -std=legacy -I$R/dis31/mcfm/Inc"
$FF -c $R/dis31/mcfm/spinoru.f -o spinoru.o
$FF -c $D/mcfm/ampqqbgll.f -o ampqqbgll.o
gfortran -O2 -ffree-line-length-none -I$D -c $D/hard21.f90 -o hard21.o
gfortran -O2 -ffree-line-length-none -I$B/modules -I. -c $D/tests/harness_hard21.f90 -o harness_hard21.o
gfortran -O2 -ffree-line-length-none -c $R/dis31/ew31.f90 -o ew31.o
gfortran -O2 -ffree-line-length-none -I$(hoppet-config --prefix)/include/hoppet -I$B/modules -I. -c $D/tests/harness_born21.f90 -o harness_born21.o
sed "s#-o tau2_nlo#-o $W/harness_hard21#; s#CMakeFiles/tau2_nlo.dir/tau2_nlo.f90.o#$W/harness_hard21.o $W/hard21.o $W/spinoru.o $W/ampqqbgll.o#; s#\$# -L$L -Wl,-rpath,$L -lnnlojet_core#" \
  $B/slicing/CMakeFiles/tau2_nlo.dir/link.txt > link.sh
(cd $B/slicing && bash $W/link.sh)
sed "s#harness_hard21#harness_born21#g; s#$W/hard21.o#$W/hard21.o $W/ew31.o#" link.sh > link_born21.sh
(cd $B/slicing && bash $W/link_born21.sh)
echo "built $W/harness_hard21 and $W/harness_born21"
