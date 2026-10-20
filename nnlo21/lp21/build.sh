#!/bin/bash
# Build test_lp21 (the tau_2 cumulant of DIS 2+1 at NNLO from MCFM 10.3's
# SCET pieces) from an MCFM 10.3 source tree (not copied into the repository).
#   nnlo21/lp21/build.sh /path/to/MCFM-10.3 <dir of libnnlojet_core.so> [builddir]
# Change to MCFM's soft1.f90 (sed below): the I_ij,m array has six entries
# (preIijm(5) -> preIijm(6)); lp21 passes all six (MCFM sets I13x2 = I23x1).
# Change to MCFM's plus.f (20 Oct audit): the subtracted plus-distribution
# integrand was dropped for 1 - z < 1e-5 (O(1e-3) of the gluon-channel Born at
# O(alpha_s^2)). 1e-8: within 1e-4 of the exact plus (1e-13, which equals an
# independent implementation to 6.5e-10); below ~1e-8 the subtraction is
# rounding noise and the adaptive z integration stalls at large xi.
set -e
M=$(cd "${1:?path to MCFM-10.3}" && pwd); L=$(cd "${2:?dir of libnnlojet_core.so}" && pwd)
B=${3:-build}
here=$(cd "$(dirname "$0")" && pwd); R=$(dirname "$(dirname "$here")")
mkdir -p "$B"; cd "$B"
for f in SCET/ap.f SCET/auxfunctions.f SCET/I1qq.f SCET/I2qq.f SCET/I1gg.f SCET/I2gg.f SCET/plus.f \
         SCET/xbeam1bis.f SCET/xbeam2bis.f SCET1j/jet.f90 SCET1j/soft1.f90 SCET1j/assemblejet.f90 Need/lnrat.f; do
  cp "$M/src/$f" .
done
for f in Li2.f Li3.f WGPLG.f dgauss.f xspenz.f xcdil.f; do cp "$M/lib/SpecialFunctions/$f" .; done
sed -i 's/preIijm(5)/preIijm(6)/g' soft1.f90
sed -i 's/abs(one-z) < 1.e-5_dp/abs(one-z) < 1.e-8_dp/' plus.f
grep -q '1.e-8_dp' plus.f
FF="gfortran -O2 -ffixed-line-length-132 -ffree-line-length-none -I$M/src/Inc -I$M -I."
gfortran -c "$R/slicing/mcfm_lp/stubs_mod.f90"
for f in *.f jet.f90 soft1.f90 assemblejet.f90; do $FF -c "$f"; done
gfortran -O2 -ffree-line-length-none -c "$R/slicing/mod_slicing_scet.f90"
gfortran -O2 -ffixed-line-length-132 -std=legacy -I$R/dis31/mcfm/Inc -c $R/dis31/mcfm/spinoru.f -o spinoru_dis.o
gfortran -O2 -ffixed-line-length-132 -std=legacy -I$R/dis31/mcfm/Inc -c $R/nnlo21/mcfm/ampqqbgll.f -o ampqqbgll_dis.o
gfortran -O2 -ffree-line-length-none -I$R/nnlo21 -c $R/nnlo21/hard21.f90
$FF -c "$here/lp21.f90"
$FF -c "$here/test_lp21.f90"
gfortran -o test_lp21 $(ls *.o | grep -v "^sliced21.o$") $(lhapdf-config --libs) -Wl,-rpath,$(lhapdf-config --libdir) $(hoppet-config --libs) \
  -L$L -Wl,-rpath,$L -lnnlojet_core
echo "built $B/test_lp21"
