# DIS 3+1 at NLO and 4+1 at LO with Catani–Seymour dipoles: plan

Started 2 Oct 2026 (branch `2026-10-dis31`, off `2026-10-tau-slicing`).

## Why

NNLO DIS 2+1 by τ₂ slicing needs, above the cut, DIS 3+1 at NLO (the third
parton resolved): 4+1 trees (double real) and 3+1 one-loop (real-virtual),
with single-unresolved subtraction only. Catani–Seymour dipoles, as in
DISENT for 2+1. Below the cut: the NNLO factorisation theorem (hard,
beam, jet, soft functions; to come). The same NNLO 2+1 serves N3LO 1+1 via
P2B in disorder and the (3,0)/(0,3) pieces of N3LO VBF.

## Ingredients

1. Trees, photon exchange first (Z, W later), all crossings with one
   incoming parton:
   - 3+1: γ* q → q g g, γ* q → q Q Q̄, γ* q → q q q̄ (identical),
     γ* g → q q̄ g; with colour-correlated ⟨T_i·T_j⟩ and gluon
     spin-correlated Borns for the dipoles.
   - 4+1: γ* q → q g g g, γ* q → q Q Q̄ g (and identical), γ* g → q q̄ g g,
     γ* g → q q̄ Q Q̄.
   Source: MCFM 10.3's Z+2 jet / Z+3 parton amplitudes (`src/Z2jet`: trees,
   real, and the colour structures its own CS dipoles use), crossed to DIS
   (the lepton pair with one incoming electron, q² < 0), as for
   `src/disent_virt3.f`. Efficiency matters at every step: spinors once
   per point, charges applied after the amplitudes; if MCFM's trees are the
   bottleneck, derive our own helicity-summed trees with FORM (e.g. a
   five-parton MATFVE in DISENT's style).
2. One loop, 3+1: MCFM's Z+2 jet virtual (BDK amplitudes, `qqb_z2jet_v`),
   crossed; renormalised in MS-bar, finite part in the CS convention.
3. Subtraction: CS dipoles for one incoming parton (FF, FI, IF) with
   four-parton Borns, and the integrated I, K, P operators; generalising
   DISENT's three-parton implementation (`src/libdisent.f`).
4. Phase space: 3+1 Born points (fixed x, Q² first, as the slicing
   drivers), 4+1 from 3+1 by inverse CS maps (multichannel over dipoles);
   lab and Breit frames as in disorder.

## Validation

1. Matrix elements point by point against NNLOJET v1.0.2's public core
   library (GPL; `libnnlojet_core.so`): B2g0Z, Bt2g0Z (q q̄ g g), C0g0Z,
   D0g0Z (four quarks), B3g0Z, Bt3g0Z, Btt3g0Z (q q̄ g g g), C1g0Z, D1g0Z,
   and at one loop B2g1Z, Bt2g1Z, Btt2g1Z, C0g1Z, D0g1Z (colour-ordered
   pieces, combined with their colour factors). No code from the private
   dis-1jet program.
2. Pointwise against Feynman diagrams where the symmetrised comparisons
   are blind (the charge-odd e_q e_Q terms of the four-quark channels:
   `dis31/tests/fd31.py`), and 4+1 against 3+1 in all single-collinear
   limits (`dis31/tests/harness_lim41.f90`).
3. Dipole limits: the 4+1 real against the sum of dipoles in all single
   collinear and soft limits (as `tests/test_subtraction.f90` for DISENT).
4. Poles of the one loop against the I operator (the finite part's μ_R
   dependence must be the renormalisation-group one).
5. Integrated NLO 3+1 cross sections and distributions against an
   independent code (NLOJet++'s dis3jet if obtainable; otherwise
   cross-checks of cutoff independence and of the IR limits).
6. Above-cut part of NNLO 2+1: τ₂ > τ_cut on all events, combined with the
   slicing below the cut, against DISENT-based NNLO 1+1 projections where
   possible.

## Steps

1. Harness: crossed MCFM trees (3+1, 4+1) vs NNLOJET, photon exchange.
2. Colour- and spin-correlated 3+1 Borns; dipoles; limit tests.
3. LO 4+1 integration (with jet cuts), checks of the cross section.
4. One loop 3+1 (crossed BDK) vs NNLOJET; I, K, P; NLO 3+1.
5. Validation 4; then the τ₂-sliced NNLO 2+1.

## Next: NNLO 2+1 below the cut (survey, 3 Oct)

The O(α_s²) terms of σ(τ₂ < τ_cut) = H × B × J × J × S for one beam and two
jets. MCFM 10.3 (`src/SCET1j`, Z+jet by 1-jettiness, the same three
coloured directions crossed) gives:
- beam functions at NNLO (`SCET/xbeam1bis.f`, `I2qq.f`, `I2gg.f`) and jet
  functions (`SCET1j/jet.f90`, GSTW (A.10)–(A.12)): usable as they are (Q_B,
  Q_J = 2E as in our NLO slicing);
- soft function (`SCET1j/soft1.f90`, geometric measure): the logarithmic
  and abelian two-loop terms are analytic in general y_ij; the non-abelian
  two-loop constant is a fit (Campbell, Ellis, Mondini, Williams,
  arXiv:1711.09984) in y31, y23 with y12 = 1 (back-to-back beams): not
  usable for DIS as it is;
- two-loop hard function (`Zampqqbgsq`, Becher–Lorentzen–Schwartz
  coefficients): timelike regions only (q² > 0), no spacelike continuation:
  not usable.
Replacements:
- hard: Gehrmann, Glover, arXiv:0904.2665, two-loop helicity amplitudes for
  (2+1)-jet DIS, with FORM/Fortran files in the arXiv sources and the
  variable transformations for the spacelike regions; check against
  NNLOJET's DIS double-virtual pointwise;
- soft: Bell, Dehnadi, Mohrmann, Rahn, arXiv:2312.11626 (SoftSERVE, NNLO
  N-jettiness soft function, numerical grids for 1- and 2-jettiness as
  ancillary files): to check whether its 1-jettiness grids (three Wilson
  lines) cover our geometry and normalisations (Breit frame, Q_i = 2E_i),
  possibly after a boost; otherwise compute with the SoftSERVE method or
  Buonocore et al., arXiv:2604.13167. Check: the τ_cut independence of the
  sliced NNLO 2+1 and its projection against DISENT-based NNLO 1+1.

Soft function, follow-up (2312.11626, from the abstract page and HTML): with
Q_i = 2ω_i in a given frame the soft function depends only on
n_ij = 1 − n̂_i·n̂_j. For three directions that is three numbers; the
1-jettiness grids (about 30,000 points, Laplace space, α_s/4π) have the
two beams back to back (n_12 = 2), a two-parameter slice. Our Breit-frame
definition has all three n_ij general, so it is not covered. Option: define
τ₂ event by event in the frame where the incoming parton and one jet are
back to back (as MCFM's `tauboost` uses the Z+jet rest frame), with
Q_i = 2E_i there. Then the three directions are "beam, back-to-back jet,
other jet", i.e. the hadronic 1-jettiness geometry, and the grids (or the
CEMW fit, where valid) apply, with the colour assignment by which direction
is the gluon. The NLO slicing (beam/jet Q_i, the I_ijm) must then move to
the same frame. To settle by reading the paper and the grids before
building.
