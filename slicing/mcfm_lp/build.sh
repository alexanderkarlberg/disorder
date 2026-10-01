#!/bin/bash
# Build dis_tau1_lp (NNLO leading-power tau_1 cumulant for DIS 1+1) from the
# SCET pieces of an MCFM 10.3 source tree (not copied into this repository).
# Usage: ./build.sh /path/to/MCFM-10.3 [builddir]
set -e
M=${1:?path to MCFM-10.3}; B=${2:-build}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$B"; cd "$B"
for f in SCET/ap.f SCET/auxfunctions.f SCET/I1qq.f SCET/I2qq.f SCET/I1gg.f SCET/I2gg.f SCET/plus.f \
         SCET/xbeam1bis.f SCET/xbeam2bis.f SCET/qqcoeff.f SCET/hardqq.f SCET0j/softqqbis.f SCET0j/assemble.f \
         SCET1j/jet.f90 Need/lnrat.f; do cp "$M/src/$f" .; done
for f in Li2.f Li3.f WGPLG.f; do cp "$M/lib/SpecialFunctions/$f" .; done
FF="gfortran -O2 -ffixed-line-length-132 -I$M/src/Inc -I$M -I."
gfortran -c "$here/stubs_mod.f90"
for f in *.f jet.f90; do $FF -c "$f"; done
$FF -c "$here/dis_tau1_lp.f90"
gfortran -o dis_tau1_lp *.o $(lhapdf-config --libs) -Wl,-rpath,$(lhapdf-config --libdir)
echo "built $B/dis_tau1_lp"
