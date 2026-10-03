# NNLO DIS 2+1 by τ₂ slicing (branch 2026-10-nnlo21)

Plan: `docs/nnlo21-plan.md`. Above the cut: NLO 3+1 (`dis31/`). Below the
cut: H × B × J × J ⊗ S at O(α_s²) in the jets'-rest-frame geometric measure
(`slicing/mod_tau2_run.f90`, `-measure cm`).

## Hard function: `hard21.f90`

Two-loop γ* → q q̄ g helicity coefficients (Garland et al.; DIS
continuation Gehrmann–Remiddi, Gehrmann–Glover arXiv:0904.2665) from
NNLOJET v1.0.2's 16 kinematic regions (`makecoef`, linked from
`libnnlojet_core.so`, GPL-3.0-or-later), contracted with the lepton current
by MCFM's `ampqqbgll` (`mcfm/`, GPL), converted from Catani's scheme to
SCET (MS-bar) at μ = Q.
- `scheme_conversion.py` derives the conversion (sympy) and checks it
  against MCFM 10.3's timelike `schemeconvC0`, `schemeconv2lM0` (equal);
  writes `x2zero.inc`.
- NNLOJET's region coefficients equal the Catani-scheme finite remainders
  of 0904.2665's Fortran files (arXiv sources; all eight DIS regions there,
  30 coefficients, ≤ 1e-6) after a one-loop convention shift (see the header
  of `hard21.f90`). 0904.2665's files cover only one of the two q ↔ q̄
  partner regions of the quark channel, which the helicity sum needs;
  NNLOJET has both.

Tests (`tests/harness_hard21.f90`, build with `build_harness.sh`):
- tree: |M0|² ∝ DISENT's MATTHR (ratio × Q⁴ constant to all digits,
  quark and gluon);
- one loop: H⁽¹⁾/H⁽⁰⁾ equals the independent DISENT-based hard function of
  the NLO slicing (hard_fact + non-factorising LEIV/ERTV) to 3e-11 (quark)
  and 1e-9 (gluon) at 40 random points;
- two loop: continuity along angular scans across the region boundaries;
  one step of 0.2% at 2 p₁·p₂ = Q² (NNLOJET's displacement of v by 1e-3 at
  v → 1), negligible after integration; to refine if needed.
