#!/bin/bash
# Build sliced21 (below-cut part of tau_2-sliced NNLO DIS 2+1) and nlo31
# (above-cut part, mode 1) into one build directory.
#   nnlo21/build_sliced21.sh /path/to/MCFM-10.3 <dir of libnnlojet_core.so> <builddir>
set -e
M=$1; L=$2; B=$3
here=$(cd "$(dirname "$0")" && pwd); R=$(dirname "$here")
"$here/lp21/build.sh" "$M" "$L" "$B" > /dev/null
B=$(cd "$B" && pwd); L=$(cd "$L" && pwd)
mkdir -p "$B/dis31"; cd "$B/dis31"
"$R/dis31/build_nlo31.sh" . > /dev/null      # nlo31 itself (above-cut part, mode 1)
D=$R/dis31
sed -n '1,/^end module nlo31_mod/p' $D/nlo31.f90 > nlo31_mod.f90
gfortran -O2 -ffree-line-length-none -I$D/mcfm/Inc -c nlo31_mod.f90
cd "$B"
FF="gfortran -O2 -ffixed-line-length-132 -ffree-line-length-none -I$M/src/Inc -I$M -I. -Idis31"
$FF -c "$here/sliced21.f90"
OBJ=$(ls *.o | grep -v -e '^test_lp21.o$' -e '^spinoru_dis.o$' -e '^lnrat.o$')   # dis31 has its own lnrat, spinoru
DOBJ=$(ls dis31/*.o | grep -v -e 'dis31/nlo31.o$')
gfortran -O2 -o sliced21 $OBJ $DOBJ $(lhapdf-config --libs) -Wl,-rpath,$(lhapdf-config --libdir) \
  -L$L -Wl,-rpath,$L -lnnlojet_core -lstdc++
echo "built $B/sliced21 and $B/dis31/nlo31"
