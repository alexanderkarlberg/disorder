#!/bin/bash
# build the matrix-element harnesses against a disorder build directory
# (objects of disorder_core and the test support), e.g.
#   dis31/tests/build.sh /path/to/build
# run with: ./harness_gluon -toyQ0 2 ; ./harness_quark -toyQ0 2 ; ./harness_me31 -toyQ0 2
set -e
B=$(cd "$1" && pwd); D=$(cd "$(dirname "$0")/.." && pwd); W=${2:-$PWD}
cd "$W"
for f in z2jetsq storecsz subqcd spinoru ampqqb_qqb aqqb_zbb; do
  gfortran -O2 -ffixed-line-length-132 -I$D/mcfm/Inc -c $D/mcfm/$f.f -o $f.o
done
O=$(ls $B/CMakeFiles/disorder_core.dir/src/*.o | grep -v "mod_dsigma\|mod_analysis")
gfortran -O2 -ffree-line-length-none -c $D/me31.f90 -o me31.o
for h in harness_gluon harness_quark harness_me31 time_me31; do
  gfortran -O2 -ffree-line-length-none -I$(hoppet-config --prefix)/include/hoppet -I$B/modules -c $D/tests/$h.f90 -o $h.o
  gfortran -O2 $h.o me31.o z2jetsq.o storecsz.o subqcd.o spinoru.o ampqqb_qqb.o aqqb_zbb.o $O \
    $B/CMakeFiles/disorder_core.dir/analysis/pwhg_bookhist-multi.f.o \
    $B/tests/CMakeFiles/disorder_test_support.dir/__/analysis/simple_analysis.f.o \
    $(hoppet-config --libs) $(lhapdf-config --libs) -o $h
done
