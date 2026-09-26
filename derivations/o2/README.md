# DISENT's O(αs²) ingredients for γ/Z and W exchange

Derivation and validation of the Z/W (parity-violating, and CC-specific)
extensions of DISENT's O(αs²) matrix elements (`src/libdisent.f`: MATFOR,
CONTHR → CONTHR3, VIRTHR, SUBFOR), following the same approach as for the
three-parton tree (`../matthr/`). Results and the logic are in
`docs/notebook.md` (2026-09-26).

## Photon re-derivation (point by point against DISENT)

* `harness/harness.f90`: evaluates DISENT's photon functions (MATFOR
  structures, CONTHR, ERTV/LEIV) and the new routines at phase-space points
  (`harness/kin.py`, `kinee.py`, `kinlim.py`). Build with e.g.
  `gfortran -std=legacy -o harness harness.f90 ../../../src/libdisent.f
  ../../../src/disent_virt3.f ../../../src/disent_o2_trees.f`.
* `gen_o2.py` (with `gen2.py`): FORM programs for the tree structures in an
  all-outgoing notation with crossing by momentum substitution, resolved in
  lepton helicity × quark chirality (projectors). `formeval.py`, `evalprog.py`
  evaluate FORM output numerically; `cmp_rest.py`, `gen_qgg.py`/`cmp_qgg.py`
  compare with DISENT.
* `spinors.py`, `num4q.py`: independent explicit-Dirac-spinor evaluation of
  the four-quark amplitudes (checks FORM's D1, D2, E and the CC classes).
* `oneloop/bdk.py`: Python port of the Bern–Dixon–Kosower one-loop
  amplitudes for 0 → q̄ q l l̄ g from MCFM; `fit_photon.py`, `fit_ee.py`,
  `fit_ee_avg.py` compare with DISENT's ERTV/LEIV.

Photon results (every point, arbitrary CF, CA): q→qgg = 1/64 [CF X1 +
(CF−CA/2) X2]; g→qq̄g 1/32; D1 (boson on the incoming line 1→2, g→(3,4))
1/32; D2 (boson on the (3,4) pair) 1/32; E −1/16 Re; CONTHR
∓CF·4π²(4π/137)²/Q⁴ (×1/VV); one loop: with all four helicities summed,
BDK = DISENT with normalisation −1/(8π²) and tree terms
−7/2+π²/2−(L13²+L23²)/2 (leading), 7/2+L12²/2 (subleading colour).

## Z/W pieces

* `make_trees.py` regenerates `src/disent_o2_trees.f` (FORM → optimised
  Fortran via `gen_prod.py`, `wrap.py`): helicity differences (PV) and sums
  of q→qgg (FQGG3), D1 (FD13), E (FE3), the CC classes EXX/EXY (FEXX3,
  FEXY3) and the spin-correlated three-parton ME (FCONTH3). The sums
  reproduce DISENT's photon functions to 1e-14.
* `src/disent_virt3.f`: Fortran port of the BDK amplitudes; its photon
  combination reproduces ERTV/LEIV to 1e-13, its PV output agrees with the
  Python port.
* The consistency of MATFOR (incl. PV) with SUBFOR (incl. CONTHR3) in all
  single-unresolved limits is `tests/test_subtraction.f90` (also
  `harness` modes `subph/subnc/subcc` with `kinlim.py`).
