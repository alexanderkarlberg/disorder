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
| 5 | Same, axial part (triangle anomaly) | α_s³ | yes, but ∝ Σ_f a_f over the loop flavours: cancels in each massless isodoublet, survives through the top–bottom mass splitting (in MCFM: `toploops`) | dropped (we have n_f = 5 massless flavours and no top) | estimate pending (section 3). Same status as in MCFM's Z + 2 jets without top loops and in NNLOJET's DIS process. |
| 6 | One loop four quarks (q q̄ Q Q̄ + V): interference of the boson on the two open lines | α_s³ | yes | dropped (as row 1–2 at one loop) | vector part odd under pair exchange (zero after integration); axial part same order of size as row 2 relative to its order: expected ~10⁻⁶ of the α_s³ coefficient. |
| 7 | One loop four quarks, boson on a closed loop (loop with V and two gluons) | α_s³ | vector: no (Furry, two gluons); axial: as row 5 | dropped (MCFM's qqb_z2jet_v sets a63z to zero too) | as row 5 |
| 8 | One-loop interference, reflection-odd part (ε-tensor × absorptive part) | α_s³ | yes, with Z only (needs c_LL ≠ c_RR); zero for the photon | kept in virt31 point by point; nlo31's flavour sums keep only the reflection-even part | zero after integration for observables symmetric under the reflection of the hadronic final state through the lepton plane (y → −y), i.e. all unpolarised, azimuthally symmetric observables we compute. **Open:** virt31 (MCFM's crossing) and NNLOJET differ in this part point by point (8 Oct). |
| 9 | 2+1 Born, gluon channel: parity-odd part | α_s (Born) | yes, point by point | DISENT's MATTHR and sliced21 drop it | odd under q ↔ q̄ exchange: zero after integration for flavour-blind observables. Not an approximation for us. |
| 10 | 2+1 one-loop hard function: boson on a closed loop (V + q q̄ g, loop with V and two gluons) | α_s² (NLO 2+1) | vector: no (Furry); axial: yes, anomaly, ∝ Σ_f a_f (top–bottom) | dropped | estimate pending (section 3). Also absent from DISENT's NC VIRTHR (validated reference for NLO 2+1) and from NNLOJET's DIS process. |
| 11 | 2+1 two-loop hard function: N_F,V term (boson on a closed loop, d^abc-type; hard21's G with N_F,γ = Σ e_q/e_q) | α_s³ (NNLO 2+1) | yes, vector part ∝ Σ_f v_f | **kept** (8 Oct, stage 3): vector couplings of the loop, effective ratio per helicity class (`sliced21` `class_weights`) | — |
| 12 | Same, axial part | α_s³ | yes, ∝ Σ_f a_f (top–bottom) | dropped | as row 10 at one order higher; estimate pending |
| 13 | Top-quark loops in general (gluon self-energies etc.) | α_s³ | yes | dropped (n_f = 5, decoupled top), standard | standard; consistent with the PDFs and α_s |

## 2. Consistency with the structure functions (P2B, N3LO 1+1)

The N3LO 1+1 result built from NNLO 2+1 (P2B or slicing) needs the 2+1 piece
and the inclusive structure functions to contain the same terms. HOPPET's NC
coefficient functions with Z: the C3 pure-singlet parts are set to zero
(`InitC3NNLO`), i.e. the axial "boson on a different line/loop" terms are
absent there too. To be documented precisely for the paper: which singlet
terms HOPPET keeps for F2/FL with Z couplings (vector Σ v_q terms, the
analogues of rows 1 and 11) and at which order (the "fl11"-type terms at
O(α_s²) and O(α_s³)); our rows 1 (vector, zero after integration) and 11
(kept) must match them. Open item.

## 3. Estimates still to do

- Rows 5, 10, 12 (axial closed-loop terms): NNLOJET contains the one- and
  two-loop V → q q̄ g pure-singlet amplitudes with vector and axial
  couplings (`src/process/Z/B1gNZax.f`, based on arXiv:2211.13596 and
  2306.10170; and the one-loop Z + 2 parton axial pieces `*Zax*`), but its
  DIS process does not call them. Plan: evaluate the one-loop 2+1 axial term
  (row 10, the lowest order one) relative to the Born at the fixed point
  (0.1, 5000), integrated over the 2+1 phase space with the τ_zQ bins, as
  the size of the dropped term in the NLO 2+1 coefficient.
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
HERA kinematics; the loop ones are [to be filled]. The same terms are absent
from the inclusive coefficient functions we combine with (HOPPET) and from
DISENT and NNLOJET's DIS processes.
