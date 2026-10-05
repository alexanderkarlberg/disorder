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

## Steps (3 Oct, evening)

1. Done: measure (NLO converges); two-loop hard function (`nnlo21/hard21.f90`,
   one loop equal to the DISENT-based one to 1e-9).
2. Below-cut cumulant per 2+1 Born (`nnlo21/lp21`), from MCFM 10.3's
   SCET1j pieces (built from an MCFM tree, not copied):
   - beam: `xbeam1bis`/`xbeam2bis`, z-integrated and tabulated in ξ at
     μ = Q, Q_B = μ; the Q_B = 2E_a dependence as an exact shift of the
     cumulant polynomial in ln(τ_cut/μ) + ln(Q_B/μ);
   - jets: `jetq`/`jetg` at Q_J = 2E_J;
   - soft: `soft_ab_*` + `soft_nab_*` (quark channel = MCFM's qgq with
     1 = jet q, 2 = jet g, 3 = beam q; gluon channel = qag with
     1 = jet q, 2 = jet q̄, 3 = beam g), y_ij in the jets' frame;
   - hard: `hard21` (the n_f,γ term per quark charge);
   - assembly: `assemblejet`, with jet 2 in place of beam b.
   Check: the O(α_s) part equals the NLO cumulant of `mod_tau2_run`
   (measure 2) pointwise.
3. Driver: 2+1 Born integration in the HERA set-up (dijets) with the
   below-cut weight; above the cut, `nlo31` with τ₂ > τ_cut in place of
   ≥ 3 jets and a 2+1 observable.
4. Validation: τ_cut independence; NNLOJET epLJJ dijet NNLO (photon).

## Slicing-adapted phase space for the above-cut integrals (5 Oct)

Cause of the small-τ_cut drift: VEGAS adaptation bias in nlo31's flat
sequential-decay phase space (near-2+1 region of relative measure ~1e-5).
Design (nlo31 mode 1, option `psmc`):
- 3+1 = 2+1 Born (log η̃, isotropic two-body) × one Catani–Seymour emission:
  FF(i,j;k), FI(i,j;a), IF(a,i;k) for all labels, 12 channels, plus the flat
  generator; emission variables (y or 1−x, z̃ or u, φ) log/logistic down to
  1e-10.
- 4+1 = 3+1 (the above) × a second emission: 30 channels, plus the flat
  4+1 generator.
- Weight 1/Σ_c α_c g_c(Φ) with every channel density from the exact inverse
  CS maps (CS phase-space factorisation, hep-ph/9605323 section 5).
- Validation: phase-space volume against the flat generator; lo slice
  against DISENT with few iterations (early iterations unbiased); r at large
  τ_cut against the existing runs; then the τ_cut test.

## Order of the next steps (agreed with AK, 5 Oct)

1. Photon exchange: τ_cut test with psmc (`runs/tcut5`), then the ZEUS-like
   dijet NNLO against NNLOJET (`validation/nnlojet_epLJJ_zeus2j.run`).
2. NC and CC (Z, γZ, W) in all new pieces, which are photon-only today: trees
   (`me31`, `me41`, `born31`: helicity couplings, q and q̄ lines separately,
   W flavour structures), one loop (`virt31`: plus boson-on-loop pieces for
   Z), two-loop hard function (`hard21`: couplings on the non-singlet part;
   Z singlet terms Σv_q and axial; Gehrmann–Tancredi 1112.1531), the 2+1 Born
   in `sliced21` (DISENT's MATTHR has NC/CC). Dipoles, I/K/P, soft, beam,
   jet functions and psmc are unchanged. Validation: pointwise against
   NNLOJET's DIS with Z and DISWM/DISWP, as for the photon. VBF needs this.
3. Then efficiency (profiling, hoppet convolutions), then P2B at N3LO.
