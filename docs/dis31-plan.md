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
   `src/disent_virt3.f`.
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
2. Dipole limits: the 4+1 real against the sum of dipoles in all single
   collinear and soft limits (as `tests/test_subtraction.f90` for DISENT).
3. Poles of the one loop against the I operator (the finite part's μ_R
   dependence must be the renormalisation-group one).
4. Integrated NLO 3+1 cross sections and distributions against an
   independent code (NLOJet++'s dis3jet if obtainable; otherwise
   cross-checks of cutoff independence and of the IR limits).
5. Above-cut part of NNLO 2+1: τ₂ > τ_cut on all events, combined with the
   slicing below the cut, against DISENT-based NNLO 1+1 projections where
   possible.

## Steps

1. Harness: crossed MCFM trees (3+1, 4+1) vs NNLOJET, photon exchange.
2. Colour- and spin-correlated 3+1 Borns; dipoles; limit tests.
3. LO 4+1 integration (with jet cuts), checks of the cross section.
4. One loop 3+1 (crossed BDK) vs NNLOJET; I, K, P; NLO 3+1.
5. Validation 4; then the τ₂-sliced NNLO 2+1.
