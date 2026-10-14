# Fiducial power corrections in slicing: literature notes

Notes taken 6 Oct 2026 (AK: "make notes about these papers so we don't
forget"), when the ZEUS-like dijet NNLO by τ₂ slicing showed no plateau
(docs/notebook.md, 5–6 Oct). Read from the arXiv PDFs; only the parts listed
here were read in detail.

## The problem in one paragraph

Slicing evaluates the observable at Born kinematics below the cut. Real events
below the cut differ from the Born by soft/collinear radiation (jet masses,
initial-state recoil ~ √τ Q). Cuts on the final state (lepton or jet p_T, η,
masses, isolation) respond linearly to these shifts, so the slicing residual
gets **fiducial power corrections ∝ √x_cut** (x = τ or q_T²/Q²) on top of the
hadronic ones ∝ x_cut. The origin is the breaking of the Born's azimuthal
symmetry by the cuts.

## Ebert, Tackmann, JHEP 03 (2020) 158 — "Impact of isolation and fiducial cuts on q_T and N-jettiness subtractions"

- First analytic study of the cut-induced power corrections. For Drell–Yan with
  symmetric lepton p_T cuts: p = 1/2, proportional to p_T^min/Q. Worse for q_T
  than for N-jettiness (up to an order of magnitude, depending on cuts and
  cutoffs).
- Proposed computing them numerically by projection (P2B) — the basis of the
  next two papers. (Ref. [45] in Campbell–Neumann–Vita.)

## Campbell, Neumann, Vita, arXiv:2408.05265, JHEP 05 (2025) 172 — "Projection-to-Born-improved subtractions at NNLO" (MCFM)

- Classification (their eqs. 2.10–2.12): hadronic power corrections
  x_cut ln^{2l−1} x_cut (l = perturbative order); fiducial ones
  √x_cut ln^{2l−2} x_cut (one log fewer).
- **P2B-improved slicing (their eq. 2.14):**
  σ(O) = σ_sub(x_cut) O(Born) + ∫_{x>x_cut} dσ_{h+j}(O)
         + ∫_{x<x_cut} dσ_{h+j}^{full}(O − Õ),
  Õ = the observable on the projected Born; dσ_{h+j} = the process with one
  more jet at one order lower (for us: the 3+1 at LO resp. NLO). The remaining
  slicing residual is that of Õ, which has no fiducial cuts dependence: only
  hadronic power corrections ∝ x_cut.
- Projection for colour singlets (following ref. [14] = Chen et al., N3LO H):
  rescale the incoming partons p̃_{a,b} = ξ_{a,b} p_{a,b}, ξ fixed by preserving
  q² and the boson rapidity; decay products boosted (their eqs. 2.17–2.19).
- Numbers (Z → ℓℓ with symmetric p_T cuts): at NLO the residual goes from
  √τ to linear τ. At NNLO the unimproved 0-jettiness slicing shows **no
  convergence at all** (our ZEUS situation); with P2B per-mille truncation
  errors; 0.5% on σ_NNLO reachable at τ_cut ≈ 2e-3 (unimproved would need
  τ_cut ≪ 1e-4); adding the NLP LL hadronic corrections: τ_cut ≈ 1e-2 suffices.
  Warning: without NLP terms the diagonal channel showed an erroneous plateau
  near 1e-3 — a local extremum can mislead.
- Photon isolation: P2B (and q_T recoil) does not remove the dominant isolation
  power corrections in the fragmentation channel; they devise an extra method
  (section 4).

## Alioli, Billis, Broggio, Stagnitto, arXiv:2504.11357, JHEP 01 (2026) 065 — "NNLO predictions with nonlocal subtractions and fiducial power corrections in GENEVA"

- NNLO in GENEVA with N-jettiness (T₀, T₁) nonlocal subtraction for colour
  singlet **and colour singlet + jet (Z+jet)**, combined with P2B to include the
  fiducial power corrections (FPCs) below the IR cutoff T_δ; validated against
  NNLOJET. The closest published analogue to our DIS 2+1 case (jets in the
  final state, 1-jettiness).
- Structure (their eqs. 2.3, 3.9, 4.16): the P2B identity
  O^NNLO(Φ_N) = dσ_N^NNLO O(Φ_N) + ∫ dσ_{N+1}^NLO [O(Φ_{N+X}) − O(Φ_N)], with
  the first term by slicing/subtraction in T_N and the second exact (NLO_{N+1}
  with FKS, RV and RR terms each with O(event) − O(projected)).
- Three mappings matter: (i) the splitting map of the NLO_{N+1} calculation,
  (ii) the projection used by the subtraction counterterm's splitting
  functions, (iii) the P2B projections Φ_{N+1,N+2} → Φ_N for the observable.
  (ii) and (iii) must preserve their generation-level restrictions (e.g. the
  defining cut on the Born jet); otherwise large weights and spikes. They
  project the "two closest partons" (with flavour- and map-validity checks).
- Takeaway for us: P2B-FPC with jets works in practice with N-jettiness;
  the choice of projection matters for convergence and must respect the cuts
  defining the Born.

## Related

- NLP (next-to-leading-power) hadronic corrections for N-jettiness: Moult,
  Rothen, Stewart, Tackmann, Zhu; Ebert, Moult, Stewart, Tackmann, Vita, Zhu
  (e.g. arXiv:1802.00456); V+1 jet NLP: arXiv:1907.12213. A different remedy
  (inclusive power corrections), complementary to P2B.
- N3LO 0-jettiness power corrections with fiducial cuts, disentangled with
  P2B-improved slicing: arXiv:2401.03017.
- Our own P2B: Cacciari, Dreyer, Karlberg, Salam, Zanderighi, PRL 115 (2015).

## What is new for disorder/nnlo21

- DIS 2+1 with jet cuts (Breit and lab frame): a projection 3+1 / 4+1 → 2+1
  that preserves x, Q², y, is IR safe, respects the Born-defining jet cuts
  where needed, and agrees with the Born of the τ₂ cumulant up to power
  corrections. First test at NLO (ZEUS-like dijets, where our slicing NLO
  converges like √τ to NNLOJET): with the P2B term it should be flat from
  much larger τ_cut. Then NNLO (vi, kp, r with O − Õ below the cut).

## Idea (AK, 6 Oct; with K. Melnikov): numerical P2B beyond two legs

P2B needs the cross section differential in the projected Born variables,
σ_proj(Φ_B) = ∫ dσ δ(Φ_B − proj(Φ)), to the full order. For DIS 1+1 these are
the structure functions (analytic). For 2+1 with an IR-safe projection
(e.g. `project21`), σ_proj is itself an IR-safe observable, so it can be
computed numerically and tabulated once, like "2+1 structure functions":
- for photon exchange, the φ dependence (lepton plane against hadron plane)
  is exactly 1, cos φ, cos 2φ, so a few functions of (x, Q², x_p, z_p);
- PDFs enter through a convolution grid (APPLgrid/fastNLO style);
- they can be filled by subtraction (DISENT, nlo31) or by slicing at fixed
  Born points, where Õ is constant: only hadronic power corrections, so
  τ → 0 extrapolation or NLP terms can be used.

Then any observable is ∫ σ_proj Õ + ∫ dσ_{N+1}(O − Õ). Proof of concept one
order lower: NLO 2+1 = σ_proj^NLO grid + LO 3+1 (O − Õ), which needs no
subtraction. Check against NNLOJET in the ZEUS set-up.

Refinement (AK, 6 Oct, same day): not a table. Start from an ordinary
subtraction and add/subtract R·Õ(proj Φ_R):
σ(O) = ∫dΦ_B Õ [B + V + I + Δ(Φ_B)] + ∫dΦ_R R (O − Õ), with
Δ(Φ_B) = ∫[R δ(Φ_B − proj Φ_R) − Σ_k D_k δ(Φ_B − Φ̃_k)].
- R(O − Õ) needs no subtraction.
- At NNLO, RR/RV (O − Õ) need only the NLO_{N+1} subtraction (nlo31).
- The double-unresolved structure sits only in the Born-local,
  observable-independent Δ^NNLO. Its counterterms need only be correct after
  integrating the radiation at fixed Φ_B: azimuthally averaged, nonlocal, or
  slicing at fixed Φ_B (hadronic power corrections only).
- Born-first phase space through the inverse of the projection, so no
  tabulation.
- Crux: one global projection that factorises the N+2 phase space onto Φ_B in
  all double-unresolved limits.
- Testable at NLO first (AK). Deferred until the slicing is finished.

## 14 Oct 2026: the ZEUS-jet drift, literature against our data

Constant slope in ln τ_cut, jets only (mostly p̄_T 8–15), boson independent, inclusive flat;
NNLOJET W agrees near τ_cut 1e-3, ours 4–5σ above at ≤ 2e-4. Candidates (arXiv numbers checked):
- Fiducial √τ·L^k from the symmetric E_T cut (Ebert–Tackmann 1911.08486; CNV 2408.05265, who warn
  of a false plateau near 1e-3; Frixione–Ridolfi hep-ph/9707345; Banfi–Dasgupta hep-ph/0312108;
  NNLOJET dijets 1703.05977). √τ·L² has a near-constant slope over 5e-3…1e-5. Against it: our ZEUS
  runs use P2B, which removes these terms; P2B − plain (11/12 Oct) is small in the hard bins, but
  p̄_T 8–15 was not examined. Diagnostic: the same events with other cut designs (p̄_T only;
  staggered 10/8 GeV; Salam–Slade 2106.08329).
- VEGAS bias of 1/σ²-weighted iterations, ∝ 1/N_ev per run, not reduced by averaging seeds
  (Lepage 2009.05112; MCFM-9 1909.09117; NNLOJET merging 2503.22804). Tested by J7v.
- Hadronic τ·L³ power corrections (1612.00450, 1612.02911, 1802.00456; measure/frame 1807.10764;
  H+jet 1906.01020): slope falls ≈ 30× over our range, not the drift.
- kT-jet vs jettiness partition (2510.03097, 2502.09327): power suppressed.
Convergence: P2B with jet cuts (2408.05265, GENEVA Z+jet 2504.11357); τ_cut normalised to the Born
hard scale, evaluated on the projected Born (1906.01020, 2504.11357); differential T_N subtraction
between τ_δ and T_off (1505.04794, eq. 2.19); unbiased run combination (frozen grid, equal weights),
adapting at the smallest τ_cut; fits τ(c3 L³ + c2 L²), with √τ terms as a check.
