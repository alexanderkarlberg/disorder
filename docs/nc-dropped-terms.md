# Photon + Z in NNLO DIS 2+1: terms that are dropped, and why

Notes started 8 Oct 2026 (AK: "we also need some notes about the potential
impact of the various dropped terms ... whether these are terms that exist
and we drop or they do not"). For the paper. Status and numbers as of the
date; the notebook (`docs/notebook.md`, 8 Oct entries) has the details of
each check.

Conventions: c_f(h, l) is the coupling of a quark line of flavour f to the
exchanged boson for quark helicity h and lepton helicity l (`dis31/ew31.f90`):
photon Q_f q_l, plus z_f(h) z_l(l) Q²/(Q² + M_Z²). Vector and axial parts
of the Z coupling of quark f: v_f = (z_f(L) + z_f(R))/2,
a_f = (z_f(L) − z_f(R))/2 (up to normalisation). With five light flavours
and no top, Σ_f v_f ≠ 0 and Σ_f a_f = a_b ≠ 0 (the b has no partner), while
Σ_f e_f = 1/3.

The guiding principle ("option 1", AK, 8 Oct): keep every term in which the
boson couples to the quark line that is connected to the incoming parton or
to an outgoing jet in the same way as in the structure functions we use
(HOPPET; disorder's inclusive NC), and drop the terms in which the boson
couples to a different quark line or to a closed quark loop whenever their
axial part would need the top quark to be consistent. This is also what
disorder's DISENT (MATFOR, VIRTHR) does with Z exchange, and what HOPPET does
for the structure functions (`InitC3NNLO` sets the singlet parts of C3 to
zero).

## 1. Classification

| # | term | order (first) | exists? | kept / dropped | effect on our observables |
|---|------|---------------|---------|----------------|---------------------------|
| 1 | Boson on two different open quark lines, interference (3+1 tree: q → q Q Q̄ g...; 4+1 tree) — vector × vector part | α_s² (3+1 tree), i.e. NLO 2+1, N3LO 1+1 | yes, point by point | dropped with Z | photon-like e_q e_Q term: odd under exchange of the pair Q ↔ Q̄ momenta, integrates to zero for any observable that does not tag quark flavour/charge (all ours). Exactly zero after integration. |
| 2 | Same, parts with an axial coupling on a line | α_s² | yes | dropped | does not vanish after integration because Σ_Q a_Q ≠ 0 over the pair flavours (no top). Computed: ≈ 2–6·10⁻⁶ of LO for dσ/dx dQ² at (0.01, 400), (0.1, 5000), (0.3, 20000) (`dis31/tests/axsum31`); ≤ 10⁻³ of the NNLO coefficient. e⁺: sign flips at low Q². |
| 3 | Identical quarks: the boson on the "pair" line in the same pairing class (|D|², |E|² of me31) | α_s² | yes | dropped (the Q = q member of rows 1–2) | included in the axsum31 numbers above. The direct × exchange interference (D·E*) is kept: there the couplings combine to c_q² (not a different line). |
| 4 | One loop 3+1 (q q̄ g g + V): boson on a closed quark loop, vector part (BDK's A6^v; loop with V and three gluons, colour d^abc, (N − 4/N)) | α_s³ (NNLO 2+1) | **yes**, ∝ Σ_f v_f (photon Σ e_q = 1/3) | **kept** since 8 Oct (before: left out, also for the photon, with the wrong claim that Furry's theorem removes it — it does for two gluons on the loop, not three). Validated against MCFM's `xzqqgg_v` | photon: −0.03 % of vi in NLO 3+1 (HERA 3-jet set-up), 3·10⁻⁶ of vi at (0.01, 400); NNLOJET's DIS process does not contain it either |
| 5 | Same, axial part (triangle anomaly) | α_s³ | yes, but ∝ Σ_f a_f over the loop flavours: cancels in each massless isodoublet, survives through the top–bottom mass splitting (in MCFM: `toploops`) | dropped (we have n_f = 5 massless flavours and no top) | **computed (8 Oct):** MCFM's BDK axial amplitudes (a64ax, a65ax; top in the heavy-mass approximation) with e⁻ couplings in nlo31's HERA 3-jet set-up: 0.0041 ± 0.0002 pb = 2.7·10⁻⁴ of vi (q q̄ g g channels, 15.1 pb) ≈ 4·10⁻⁵ of LO (93.6 pb); independent of m_t (173 vs 500 GeV: 4.05 vs 4.27·10⁻³ pb), as expected for this non-decoupling term. Same status as in MCFM's Z + 2 jets without top loops and in NNLOJET's DIS process. |
| 6 | One loop four quarks (q q̄ Q Q̄ + V): interference of the boson on the two open lines | α_s³ | yes | dropped (as row 1–2 at one loop) | vector part odd under pair exchange (zero after integration); axial part same order of size as row 2 relative to its order: expected ~10⁻⁶ of the α_s³ coefficient. |
| 7 | One loop four quarks, boson on a closed loop (loop with V and two gluons) | α_s³ | vector: no (Furry, two gluons); axial: as row 5 | dropped (MCFM's qqb_z2jet_v sets a63z to zero too) | as row 5 |
| 8 | One-loop interference, reflection-odd part (ε-tensor × absorptive part) | α_s³ | yes, with Z only (needs c_LL ≠ c_RR); zero for the photon | kept in virt31 point by point; nlo31's flavour sums keep only the reflection-even part | zero after integration for observables symmetric under the reflection of the hadronic final state through the lepton plane (y → −y), i.e. all unpolarised, azimuthally symmetric observables we compute. virt31 (MCFM's crossing) and NNLOJET have it with opposite signs point by point: NNLOJET's one-loop Z and W functions equal ours with all helicities flipped (9 Oct; the 8 Oct "differ" was an incomplete comparison). Which sign is physical is not settled. |
| 9 | 2+1 Born, gluon channel: parity-odd part | α_s (Born) | yes, point by point | DISENT's MATTHR and sliced21 drop it | odd under q ↔ q̄ exchange: zero after integration for flavour-blind observables. Not an approximation for us. |
| 10 | 2+1 one-loop hard function: boson on a closed loop (V + q q̄ g, loop with V and two gluons) | α_s² (NLO 2+1) | vector: no (Furry); axial: yes, anomaly, ∝ Σ_f a_f (top–bottom) | dropped | (α_s/2π)·c × the 2+1 Born with c = −5·10⁻⁵ (x = 0.01, Q² = 400), −3·10⁻⁴ (0.05, 1000), −4·10⁻³ (0.1, 5000), −7·10⁻³ (0.2, 10⁴), e⁻, τ_zQ ∈ [0.05, 0.5) (section 3): ≤ 1.4·10⁻⁴ of the 2+1 Born. Also absent from DISENT's NC VIRTHR (validated reference for NLO 2+1) and from NNLOJET's DIS process. |
| 11 | 2+1 two-loop hard function: N_F,V term (boson on a closed loop, d^abc-type; hard21's G with N_F,γ = Σ e_q/e_q) | α_s³ (NNLO 2+1) | yes, vector part ∝ Σ_f v_f | **kept** (8 Oct, stage 3): vector couplings of the loop, effective ratio per helicity class (`sliced21` `class_weights`) | — |
| 12 | Same, axial part | α_s³ | yes, ∝ Σ_f a_f (top–bottom) | dropped | as row 10 at one order higher; estimate pending |
| 13 | Top-quark loops in general (gluon self-energies etc.) | α_s³ | yes | dropped (n_f = 5, decoupled top), standard | standard; consistent with the PDFs and α_s |

## 2. Consistency with the structure functions (P2B, N3LO 1+1)

The N3LO 1+1 result built from NNLO 2+1 (P2B or slicing) needs the 2+1 piece
and the inclusive structure functions to contain the same terms. HOPPET
(local checkout `~/work/hoppet`, d891044; `src/structure_functions.f90`,
`structure_function_general_full`; `src/coefficient_functions_holder.f90`;
MVV parametrisations in `src/param-coefs/`), read 8 Oct:

- **Flavour weights.** For every order F2 and FL are built as
  Σ_q (q + q̄) C ⊗ ... with the per-flavour weights e_q² (photon),
  v_q² + a_q² (Z) and e_q v_q (γZ); F3 as Σ_q (q − q̄) C3 with a_q v_q and
  e_q a_q. These are boson-on-the-incoming-line weights, i.e. option 1.
- **O(α_s²) (NNLO 1+1).** C2, CL: non-singlet, pure-singlet and gluon
  pieces, no fl11 terms (they start at O(α_s³)). C3: the singlet parts are
  set to zero (`InitC3NNLO`: "singlet piece should have no impact on Z
  case"). So the axial boson-on-another-line terms (our rows 2, 3) are
  absent, and the vector ones (row 1) integrate to zero anyway: consistent.
- **O(α_s³) (N3LO 1+1).** C2 and CL have the fl11 terms (`C2N3LO_fl11`,
  `CLN3LO_fl11`: MVV's terms with the boson on a closed quark loop / on
  another line, d^abc-like). MVV write them for photon exchange as a
  non-singlet piece with the per-flavour weight e_q² × `FL` = e_q² · 3⟨e⟩
  (nf = 5: 0.2), a pure-singlet piece with ⟨e⟩²/⟨e²⟩ − 3⟨e⟩ (`FLS` − `FL`)
  and a gluon piece with ⟨e⟩²/⟨e²⟩ (`FLG`), all with the same functions
  (`xc2ns3p.f`, `xc2sg3p.f`, `xclns3p.f`, `xclsg3p.f`).
  `structure_function_general_full` then weights the sum like every other
  term: e_q² (photon), v_q² + a_q² (Z), 2 e_q v_q (γZ), and leaves it out for
  W.
- **Photon: exact (checked 9 Oct).** The diagram coupling is e_q Σ_f e_f
  (quarks) and (Σ_f e_f)² (gluon). MVV's decomposition reproduces it for
  every flavour, because 3e_q² − e_q = 2/3 for both up- and down-type
  charges, so the non-singlet mismatch is flavour-independent and the
  pure-singlet piece absorbs it (exact rational check for nf = 3, 4, 5;
  numerically, HOPPET's photon fl11 equals the exact couplings to 1e-9).
  Our rows 1, 4, 11 (e_q Σ e_f per flavour) are therefore consistent with
  HOPPET for the photon.
- **Z and γZ: approximate.** The diagram couplings for F2 and FL are
  g_V,q Σ_f g_V,f (Z) and e_q Σ_f g_V,f + g_V,q Σ_f e_f (γZ) for quarks,
  (Σ g_V)² and 2 Σe Σg_V for the gluon (the loop couples through its
  vector part; the axial part of the open line times the vector loop is
  parity-odd and goes into F3, the axial loop is the anomaly term of rows
  5, 10, 12). v_q² + a_q² is not linear in g_V,q, so no identity like the
  photon's helps: HOPPET's Z and γZ fl11 terms are not the diagram
  couplings, and for γZ they even have the opposite sign. Size (program
  `~/work/disorder-comparisons/hoppet-fl11/fl11z.f90`, outside the
  repository; NNPDF30_nlo_as_0118, sin²θ_W = 1 − M_W²/M_Z², μ = Q):

  | x, Q² | F2^γZ fl11 coefficient HOPPET / exact | difference × (α_s/2π)³ / LO | F2^Z, same | FL: γZ, Z |
  |---|---|---|---|---|
  | 0.01, 400 | −0.37 / +1.66 (rest of N3LO: −159) | −3.2·10⁻⁵ | 7.0·10⁻⁶ | −1.1·10⁻⁵, 2.8·10⁻⁶ |
  | 0.1, 5000 | −0.15 / +0.80 (−2.7) | −2.0·10⁻⁵ | 1.2·10⁻⁶ | −6.8·10⁻⁶, 3.6·10⁻⁷ |
  | 0.3, 20000 | −0.04 / +0.28 (−14.5) | −1.4·10⁻⁵ | −1.2·10⁻⁶ | −4.3·10⁻⁶, −4.5·10⁻⁷ |
  | 0.01, 10000 | −0.46 / +1.98 (−132) | −1.3·10⁻⁵ | 3.3·10⁻⁶ | −4.7·10⁻⁶, 1.3·10⁻⁶ |

  I.e. ≲ 3·10⁻⁵ of the γZ structure function, ≲ 10⁻⁵ of the Z one; ≈ 1–2%
  of the N3LO coefficient of F2^γZ, ≲ 0.3% for F2^Z. In the e±p cross
  section F2^γZ comes with g_V^e (≈ −0.05) and the propagator factor, so the
  effect there is below 10⁻⁶. Negligible, but a real inconsistency: for a
  P2B/slicing N3LO 1+1 with Z, our NNLO 2+1 closed-loop terms (rows 4, 11,
  diagram couplings) and HOPPET's fl11 must use the same couplings. The
  clean fix is in HOPPET (weight the non-singlet fl11 piece with the
  diagram coupling instead of e_q² · 3⟨e⟩, and the gluon piece with the
  loop couplings); to be raised with the HOPPET authors. APFEL++ (checked
  9 Oct, github vbertone/apfelxx 27deaec) does the same: its NC builders
  switch the fl11 pieces on (`C23nsp{nf}` etc., off for CC) and weight
  everything with one per-flavour effective charge, `ElectroWeakCharges` =
  e_q² − 2 e_q v_q v_e P_Z + (v_e² + a_e²)(v_q² + a_q²) P_Z². So the
  HOPPET–APFEL++ agreement does not test this. HOPPET's weighting dates from
  2016 (1c110a8, unchanged since); the 2023 fixes after V. Bertone's
  comparison (1058001, 8c4ac2a: δ(1−x) and fl11 parts of the coefficient
  functions themselves) are unrelated.
  yadism (NNPDF, 0.13.11 = github 8bbfdc8; EKO only evolves) has a dedicated
  fl11 coupling, ⟨Q_b⟩ Q_b′ (loop average × line), which is the right
  structure, but `partonic_coupling_fl11` takes photon/Z from the mode, not
  from the position: the γZ orderings come out as ⟨Z⟩Z_q and ⟨γ⟩γ_q, and ZZ
  adds ⟨a⟩a_q (axial loop with the vector-loop function). Checked by running
  yadism's own routine; error up to ≈ 6·10⁻⁵ of F2^γZ. Its F3 fl02 weights
  are the loop sums (correct), as are APFEL++'s (valence channel × ΣCh).
- F3 at O(α_s³): HOPPET's C3 uses MVV's non-singlet plus/minus functions
  without the fl02 (d^abc d_abc) term for the NS± combinations and the
  valence function with it (`C%NS_V = cfN3LO_F3NS_val`, `V = 1`) for the
  total valence Σ(q − q̄), as MVV prescribe (the fl02 term has both bosons
  on a closed loop and is C-odd, so it multiplies the total valence). For
  the record: a 9 Oct version of this paragraph claimed C3 has no fl02 term
  at all; that was wrong (a missed line), the 8 Oct statement was right.
  **Couplings of fl02 with Z and γZ: correct (checked 9 Oct).** fl02 has
  both bosons on the closed loop (MVV 0812.4168, Fig. 1), so its coupling
  is the flavour trace over the loop, Σ_f 2 g_V,f g_A,f (Z), Σ_f 2 e_f g_A,f
  (γZ), (n_f/2) per W charge, times the flavour-blind total valence.
  HOPPET adds the valence piece equally to every flavour column, so its
  per-flavour weights 2 v_q a_q, 2 e_q a_q sum to exactly that loop sum:
  numerically equal to 7 digits (`fl02z.f90`, x = 0.01/0.1/0.3). Only fl11
  needs the patch.

## 3. Estimates still to do

Done: row 5 (above) and row 10 (10 Oct): MCFM 10.3's Z + jet axial
amplitude (`A53` in `virt5`/`A5NLO`, large-m_t expansion as for row 5, with
qqb_z1jet_v's couplings), crossed to DIS (all-outgoing momenta, lnrat
continuation), against the 2+1 Born from the same amplitudes (MCFM's
virtual and tree normalisations related through the leading-colour double
pole), PDF-weighted over the Born phase space at fixed (x, Q²) in τ_zQ ∈
[0.05, 0.5) (`disorder-comparisons/axial21/axsum21.f90`). Coefficients of
α_s/2π relative to the Born: see row 10; growing with Q² (top term ∝
Q²/m_t² and the Z propagator). Row 12 (the same at two loops) is expected
to have the same relative size one order higher: negligible. The program is a scratch harness (MCFM 10.3's
`xzqqgg_v`, `fax`, `faxsl` with virt31's crossing and nlo31's phase space;
outside the repository).

- Rows 5, 10, 12 (axial closed-loop terms): NNLOJET contains the one- and
  two-loop V → q q̄ g pure-singlet amplitudes with vector and axial
  couplings (`src/process/Z/B1gNZax.f`, based on arXiv:2211.13596 and
  2306.10170, with the top mass through ln(m_t²/s); and the one-loop Z + 2
  parton axial pieces `*Zax*`), but its DIS process does not call them,
  and its pure-singlet coefficients (`helcoeffPS`) exist only for the four
  timelike regions (s45 > 0; DIS kinematics stops with "kinematical region
  not implemented"). An estimate in DIS kinematics therefore needs the
  analytic continuation of those coefficients to the spacelike regions (as
  hard21 has for the main coefficients) — a separate piece of work. For row
  5 (one loop 3+1) MCFM's BDK axial amplitudes (`a64ax`, `a65ax`, `fax`,
  `faxsl`; analytic continuation through lnrat, valid in any region) could
  be evaluated in DIS kinematics directly.
- Row 6 at one loop (axial part): could be computed as axsum31 with virt31
  (`virt31_keepint`), if needed.
- Row 8: a third reference (or an analytic check of the crossing of the
  ε-tensor terms) to settle which code is right point by point. It does not
  affect any result of ours.

## 4. Summary for the paper (draft)

All terms we drop with Z exchange are terms in which the boson couples to a
quark line other than the one carrying the incoming parton or a jet's
flavour, or to a closed quark loop. The vector closed-loop terms that do
not vanish (one loop q q̄ g g, two-loop N_F,V in the 2+1 hard function) are
kept. The vector parts of the dropped terms vanish identically (Furry) or
after integration over flavour-blind final states (odd under the exchange
of the quark and antiquark of a pair); their axial parts survive
only because the five-flavour theory has no top partner for the b quark
(Σ a_q ≠ 0). The tree-level ones are 2–6·10⁻⁶ of the LO cross section at
HERA kinematics; the one-loop axial closed-loop term of the 3+1 virtual is
3·10⁻⁴ of the virtual (4·10⁻⁵ of LO) in a HERA 3-jet set-up; the 2+1 ones
(one and two loop) are not yet computed in DIS kinematics. The same terms are absent
from the inclusive coefficient functions we combine with (HOPPET) and from
DISENT and NNLOJET's DIS processes. One mismatch remains for an N3LO 1+1
with Z: HOPPET's inclusive fl11 terms (MVV's photon closed-loop/other-line
terms) carry the ordinary Z and γZ flavour weights instead of the diagram
couplings, which is exact for the photon but not for Z; the difference is
≲ 3·10⁻⁵ of the γZ structure function and below 10⁻⁶ in the e±p cross
section at HERA (section 2).
