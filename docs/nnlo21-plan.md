# NNLO DIS 2+1 by τ₂ slicing: plan

Started 3 Oct 2026 (branch `2026-10-nnlo21`, off `2026-10-dis31`, with the
DISENT fix 5262826 cherry-picked). AK: "If it all looks right I suggest you
move on to the 2+1 NNLO".

σ_NNLO(2+1) = σ(τ₂ > τ_cut) + σ(τ₂ < τ_cut) + O(τ_cut):
- above the cut: NLO 3+1 (`dis31/nlo31`, validated against NNLOJET on
  3 Oct) with τ₂ > τ_cut, plus the 2+1 observable;
- below the cut: the factorisation theorem H × B × J × J ⊗ S at O(α_s²),
  times the 2+1 Born.

## The measure: geometric in the jets' rest frame

τ₂ = T₂/Q with, for each partition π of the outgoing partons into beam /
jet 1 / jet 2, u = (P_J1 + P_J2)/m and

    T_π = Σ_beam P·p_k / P·u + Σ_jets [u·P_J − |P_J|_u],   T₂ = min_π T_π

(`measure = 2`, `-measure cm`, in `slicing/mod_tau2_run.f90`). In every
singular limit the frame is the Born's partonic CM frame, where the two
jets are back to back. Q_i = 2E_i there.

Why this measure:
- The NNLO soft function for three directions with two of them back to
  back depends on a single angle. Bell, Dehnadi, Mohrmann, Rahn
  (arXiv:2312.11626) give it as numerical grids (1-jettiness:
  `1JN.dipoles.renormalised.csv`, colour-stripped dipoles c_ij^(1),
  c_ij^(2,CA), c_ij^(2,nf), 24 angles). MCFM's CEMW fit
  (arXiv:1711.09984, `SCET1j/soft1.f90`) is an independent computation
  of the same function.
- **Change of route** (against the notebook of 2 Oct, which named the
  invariant measure with our own implementation of arXiv:2604.13167):
  - With the invariant measure, the three normalisations differ (beam
    against jets), so the soft function depends on two variables and is
    not tabulated anywhere.
  - The jets'-frame geometric measure needs only the published 1D
    function.
  - 2604.13167 stays the fallback, and a cross-check if the grids turn out
    too coarse.

Checks of the measure:
- unit (scratch `jetframe/jf_test.f90`): boost invariance to 1e-15; soft
  limit equal to min_i n_i·k in the Born CM frame;
- NLO 2+1 slicing against DISENT (`tau2_nlo -measure cm`, x = 0.01,
  Q² = 400): in progress.

## Ingredients below the cut

1. **Hard function**, two loops, γ* q → q g and γ* g → q q̄ (spacelike
   q²).
   - Source: Gehrmann, Glover, arXiv:0904.2665, finite remainders (Catani
     subtraction) in eight kinematic regions, Fortran in the arXiv sources
     (`~/work/disorder-comparisons/nnlo21/gg0904/fortran`, not copied into
     the repository). Needs `hplog.for` and `tdhpl.f` (Gehrmann–Remiddi;
     tdhpl is in MCFM).
   - Convert Catani's scheme to the SCET (MS-bar) hard function, as
     MCFM's Z+jet does.
   - Check against our one loop (virt31 crossed, DISENT's VIRTHR), and the
     two loop against NNLOJET's DIS double virtual pointwise.
2. **Beam functions** at NNLO (quark, gluon): MCFM `SCET/xbeam1bis.f`,
   `I2qq.f`, `I2gg.f` (Gaunt, Stahlhofen, Tackmann). Q_B = 2E_a in the
   jets' frame.
3. **Jet functions** at NNLO: MCFM `SCET1j/jet.f90` (GSTW (A.10)–(A.12)).
4. **Soft function** at NNLO:
   - the logarithmic terms are fixed by RG consistency (analytic, MCFM
     `soft1.f90`);
   - the constant c^(2) comes from the 2312.11626 grids, interpolated in
     the angle, with the asymptotics of their section 4 near the edges,
     and is checked against the CEMW fit.
   - DIS 2+1 has three coloured legs, so the tripoles vanish by colour
     conservation.
5. **Assembly**: the cumulant to O(α_s²) in L = ln(τ_cut/μ) (MCFM
   `SCET1j/assemblejet.f90` structure), μ_R, μ_F logs.

## Validation

1. Each ingredient pointwise as above; the NLO cumulant of the new
   measure against DISENT (step 1 of the NLO 2+1 test).
2. τ_cut independence of NNLO 2+1 = below + above.
3. Against NNLOJET's DIS dijet NNLO (epLJJ, photon exchange): the inclusive
   dijet cross section and distributions at HERA kinematics.
4. The projection to 1+1 at N3LO later (P2B with disorder).
