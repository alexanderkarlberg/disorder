#!/bin/bash
# build the matrix-element harnesses against a disorder build directory
# (objects of disorder_core and the test support), e.g.
#   dis31/tests/build.sh /path/to/build [workdir]
# run with: ./harness_gluon -toyQ0 2 ; ./harness_quark -toyQ0 2 ;
#   ./harness_me31 -toyQ0 2 ; ./harness_lim41 ; ./harness_ccsum ; ./harness_born31 ; ./harness_dip41 ; ./harness_virt31 ; ./harness_iop31 ; ./dump_fd31 && python3 fd31.py
set -e
B=$(cd "$1" && pwd); D=$(cd "$(dirname "$0")/.." && pwd); W=${2:-$PWD}
cd "$W"
M="z2jetsq storecsz subqcd spinoru ampqqb_qqb aqqb_zbb xzqqggg amp_qqggg msq_gqqQQg makemb_photon nagyqqqqg subqcdn spinork checkndotp"
for f in $M; do
  gfortran -O2 -ffixed-line-length-132 -I$D/mcfm/Inc -c $D/mcfm/$f.f -o $f.o
done
MO=$(for f in $M; do echo -n "$f.o "; done)
O=$(ls $B/CMakeFiles/disorder_core.dir/src/*.o | grep -v "mod_dsigma\|mod_analysis")
gfortran -O2 -ffree-line-length-none -c $D/ew31.f90 -o ew31.o
gfortran -O2 -ffree-line-length-none -c $D/me31.f90 -o me31.o
gfortran -O2 -ffree-line-length-none -c $D/me41.f90 -o me41.o
gfortran -O2 -ffree-line-length-none -c $D/born31.f90 -o born31.o
gfortran -O2 -ffree-line-length-none -c $D/dip41.f90 -o dip41.o
for h in harness_gluon harness_quark harness_me31 time_me31; do
  gfortran -O2 -ffree-line-length-none -I$(hoppet-config --prefix)/include/hoppet -I$B/modules -c $D/tests/$h.f90 -o $h.o
  gfortran -O2 $h.o me31.o ew31.o $MO $O \
    $B/CMakeFiles/disorder_core.dir/analysis/pwhg_bookhist-multi.f.o \
    $B/tests/CMakeFiles/disorder_test_support.dir/__/analysis/simple_analysis.f.o \
    $(hoppet-config --libs) $(lhapdf-config --libs) -o $h
done
# standalone (MCFM routines only)
for h in harness_lim41 dump_fd31 harness_born31 harness_dip41 harness_ccsum; do
  gfortran -O2 -ffree-line-length-none -c $D/tests/$h.f90 -o $h.o
  gfortran -O2 $h.o me31.o ew31.o me41.o born31.o dip41.o $MO -o $h
done
# one loop (BDK amplitudes from MCFM, dis31/mcfm/loop)
L=""
for f in $D/mcfm/loop/*.f; do
  gfortran -O2 -ffixed-line-length-132 -std=legacy -I$D/mcfm/Inc -c $f -o loop_$(basename $f .f).o
  L="$L loop_$(basename $f .f).o"
done
gfortran -O2 -ffree-line-length-none -I$D/mcfm/Inc -c $D/virt31.f90 -o virt31.o
gfortran -O2 -ffree-line-length-none -c $D/tests/harness_virt31.f90 -o harness_virt31.o
gfortran -O2 harness_virt31.o me31.o ew31.o born31.o virt31.o $MO $L -o harness_virt31
gfortran -O2 -ffree-line-length-none -c $D/iop31.f90 -o iop31.o
gfortran -O2 -ffree-line-length-none -c $D/tests/harness_iop31.f90 -o harness_iop31.o
gfortran -O2 harness_iop31.o ew31.o born31.o virt31.o iop31.o $MO $L -o harness_iop31
cp $D/tests/fd31.py .
