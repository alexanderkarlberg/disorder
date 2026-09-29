# disorder development notebook

Running log of work on disorder: what was done, what was checked and how,
what was found, and what was corrected. Newest entries at the bottom.
Corrections to earlier entries are made by adding a note, not by editing
the original text away.

---

## 2026-09-25 — Unit tests, ctest-driven validation, CI keepalive

Branch `2026-09-ci-unit-tests`. Goal: unit tests for phase-space
generation, boosts, the matrix element and the disent interface; run the
validation through ctest; stop CI being switched off after 60 days of
inactivity.

### Baseline

- Full validation (`validate_or_generate.sh validate`) on a clean clone of
  `main` (67fa107) **failed on this machine** (thA371a) before any change:
  every run using the exclusive analysis aborted with `std::length_error` in
  `fastjet::ClusterSequence::_fill_initial_history`.
- Cause (local environment, not the code): stale FastJet 3.5.1 headers in
  `~/.local/include/fastjet`, configured *without* thread safety and with no
  library next to them. LHAPDF's include dir (`~/.local/include`) comes first
  on the compile line, so they shadow `/usr/include/fastjet`, whose library
  (`/usr/lib64/libfastjet.so`) is built *with* thread safety → ABI mismatch.
  Prepending the FastJet include dir does not help (GCC drops `-I/usr/include`
  as a system dir; CMake puts includes before `CMAKE_CXX_FLAGS`).
- Workaround used for all local builds: a `g++` wrapper that adds
  `-I<dir with a symlink to /usr/include/fastjet>` first
  (`-DCMAKE_CXX_COMPILER=<wrapper>`). With it the baseline passes:
  276 comparisons, 0 failures, 3m52s wall. `~/.local` left untouched;
  suggested to move that directory aside.
- Timings: toy-PDF N3LO runs 5–20 s each; the MSHT20an3lo
  `-pdfuncert -alphasuncert` run ~310 s CPU, i.e. almost all of the time.
- `ref_runs_quick/` was byte-identical to the corresponding `ref_runs/` files
  (154 files); `quick` = full matrix minus that one run.

### Refactoring (no change to results)

- `user`, `dis_cuts`, `disent_muf`, `DOT` moved from `program disorder`
  (internal procedures, not callable from tests) to
  `src/mod_disent_interface.f90`. The momentum reordering, the x/y/Q²
  reconstruction and the mapping of disent's 3 μF points onto the 7-point
  variation are small helpers/parameters (`disent_to_momenta`,
  `disent_kinematics`, `disent_imuf`/`disent_imur`) used by `user`.
- `read_PDF` and the HOPPET/structure-function start-up moved to
  `src/mod_pdf_setup.f90` (`setup_structure_functions`).
- CMake: physics code built once as object library `disorder_core`, linked by
  `disorder` and by the tests.
- Check: all `.dat` outputs of the validation matrix (cross sections and
  histograms) **bitwise identical** to `ref_runs/` after the refactor.

### Scale choice 4

- Code (`muR_muF`) used μ = Q·sqrt((1−x)/x); header/help text and the unused
  `muRlcl`/`muFlcl` said Q(1−x)/x. Confirmed by AK: μ² = Q²(1−x)/x, i.e. the
  code is right. Fixed the texts, deleted the unused functions. Test added.
- The 1.0.0 manual lists the options as 2: Q²(1−y), 3: Q²(1−x)/x, whereas the
  code uses 3 and 4 (2 is Q with scale variations). Manual not edited.

### Unit tests (`tests/`, ctest label `unit`, ~1 s total)

- `test_boosts`: q from physical HERA kinematics (random x, y, azimuth; both
  `isPlus`). Round trips, invariants, q → (0,0,0,−Q), x·P → (Q/2)(1,0,0,1),
  z-mirror symmetry. First version used random q with q⁰<0 and tiny
  q⁰−q³ (unphysical); rounding error measured to grow like A²·1e-16 with A
  the largest matrix element of the transformation → tolerance
  1e-15·(1+A)².
- `test_phase_space`: kinematics at a grid of points; ∫Jacobian over the unit
  square = allowed area in (x,Q²) (1e-4), = Q² range at fixed x, = x range at
  fixed Q² (fixed-Q² integrand has steps → 1D midpoint rule with 2e5 points).
- `test_matrix_element`: run with disorder's flags, 13 processes at LO +
  2 beyond LO. At LO compared with parton-model expressions written out in
  the test: NC γ/Z/γZ from H1 1206.7007 eqs. 1–7; CC and ν-NC from V−A
  (ε_L = T3 − e sin²θ, ε_R = −e sin²θ; no b→t; charged lepton averaged over
  2 helicities, neutrino 1). All agree to 1e-5. Reduced cross sections
  checked against 1206.7007. Beyond LO (γ only): combination
  Y+F2 − y²FL of HOPPET's structure functions.
- `test_disent_interface`: `DISENTFULL` (2000 events, up to O(αs²)) with a
  checking callback: conservation, masslessness, cuts, Q² = xys, lab-frame
  beams, lab lepton = projected-Born lepton. Plus `dis_cuts`,
  `disent_muf`, scale-point mapping (disent's MUF² = 1, 4, 1/4 in
  `KPFUNS_SCL_VAR`).
- `guard_*`: invalid flag combinations refused (message match, since `stop`
  exits with 0).
- Mutation check (20 deliberate bugs in boosts, phase space, disent
  interface, EW couplings, F3/FL signs, propagators, scale choices): all
  caught. Note: the first round hit the unused `eval_matrix_element` /
  `compute_sigmas` (first occurrence in the file) — dead code, left alone.
  The F_L sign was only caught after adding the NLO case (F_L = 0 at LO).

### Validation through ctest

- `validation/configurations.txt` (labels | prefix | args) = the old
  `cmdline_full` matrix; the MSHT20 PDF-uncertainty run labelled `slow`.
- `validation/run_validation.py` runs one configuration and compares all
  output files. Findings while making it robust:
  - Old logs were written by GNU parallel = all stdout, then all stderr. The
    runner now captures both and writes them in that order.
  - FastJet and LHAPDF banners come from C++ streams; their position relative
    to the Fortran output depends on buffering → dropped.
  - The echoed command line contains the executable path → runner calls
    disorder through a symlink `./disorder`, and the path is normalised.
  - Numbers can be glued to text (`=245788719.9890 %`) → numbers extracted
    anywhere in a line.
  - VEGAS χ²/iteration after iteration 1 is 0 up to rounding (`0.29E-08`
    vs `0.0`) → ignored.
- Tolerance calibration: build with `-O2 -march=native` (FMA), no LTO, vs
  refs: xsct ≤ 3.8e-9, histograms ≤ 5e-8, except one P2B bin at 9.9e-7
  (value −1.8 with MC error 35). Default rtol = 1e-5. Checked that a 1e-5
  change in a cross section and a changed histogram bin are flagged.
- Final check from a clean copy of the working tree:
  `validate_or_generate.sh validate` → 45/45 tests passed, 3m48s wall.
  (Two earlier attempts at the clean copy were broken by the copy itself:
  `git ls-files` escapes the non-ASCII `μ` in names unless `-z` is used.)

### CI

- Workflow configures with `-DNEEDS_FASTJET=ON
  -DANALYSIS=exclusive_lab_frame_analysis.f` and runs `ctest` (previously the
  `ctest` step was a no-op and the script ran separately). GNU parallel no
  longer installed. `workflow_dispatch` added.
- `keepalive` job (scheduled runs only): `gh api --method PUT
  repos/<repo>/actions/workflows/cmake-single-platform.yml/enable` with
  `actions: write`. Re-enabling resets GitHub's 60-day inactivity timer (the
  mechanism used by the common keepalive actions). **Not yet verified on
  GitHub.**

### Corrections

- Suggested earlier a test "P2B total = inclusive total". Wrong idea: in P2B
  every disent event is filled with a projected-Born counter-event of
  opposite weight, so disent adds exactly zero to the total, which comes from
  the structure functions. The identity holds by construction. Replaced by
  the `DISENTFULL` checking-callback test.
- Reported (from memory) that the reduced CC cross section was 2× the
  H1/ZEUS definition. **Wrong.** H1 1206.7007 eq. 11 (cited by code and
  manual) defines σ̃_CC = 4πx/G_F² ((M_W²+Q²)/M_W²)² d²σ/dxdQ²; the code,
  the manual (eq. reducedsigma) and the published table follow it. The 2πx
  normalisation is that of the HERA combination (1506.06042), half of this.
  Test corrected.

### Open items / observations

- `ref_runs/inclusive_cc_Qmin_1_x_0.01_xsct` records a QCD scale uncertainty
  (+) of 245,788,720 % (and 78,657 % / 105,438 % on the reduced ones). Not
  from these changes; the configuration integrates Q down to 1 GeV with a toy
  PDF starting at Q0 = 2 GeV. Worth a look.
- Unused `eval_matrix_element` and `compute_sigmas` remain in
  `mod_matrix_element.f90`.
- The manual still describes validation with GNU parallel; changes are
  documented in `docs/release-notes.md` instead.

## 2026-09-25/26 — Cross-check against NNLOJET 1.0.2 and POWHEG-BOX-RES (DIS)

Full log, scripts, results and plots: `~/cernbox/disorder-comparisons/`
(`NOTES.md`, `results/`, `plots/`, `bug_reports/`). disorder used from the
command line only; no change to `src/`.

Setup: 27.5 × 920 GeV, 150 < Q² < 15000 GeV², 0.1 < y < 0.9,
NNPDF40MC_nlo_as_01180 (NNLO set for O(αs²)), α = 1/137, M_Z = 91.1876,
M_W = 80.398, zero widths, identity CKM, μR = μF = Q. The observables are in
`analysis/cmp_obs_core.f`, shared by the disorder analysis
`analysis/cmp_nnlojet_powheg.f` and a POWHEG analysis: σ, Q², x, y; the
leading lab-frame anti-kt (R = 1) jet (pT, η); and the Breit-frame,
E-normalised current-hemisphere shapes τ_zE, B_zE, ρ_E (NNLOJET definitions).
These files are not committed yet.

### Findings

- LO: all processes (NC γ/Z/γZ e∓, ν, ν̄; CC e∓, ν, ν̄) agree for all
  histograms. χ²/n ≲ 1. Totals agree to 1e-5–3e-4, within 1.4σ.
- O(αs):
  - Inclusive: disorder, NNLOJET and POWHEG agree for every process each code
    supports, with totals within 1.1σ.
  - Lab jets: NNLOJET and POWHEG agree, and the photon case also agrees with
    disorder P2B.
  - Event shapes: all three codes agree. The exceptions are NC e⁺, where
    NNLOJET returns the e⁻ result (its limitation), and one 4σ bin.
- O(αs²), photon only (disorder P2B `-nnlo` vs NNLOJET epLJJ NLO):
  - ρ_E agrees.
  - τ_zE and B_zE disagree in their tails, by up to 7%.
  - Cause: without a current-hemisphere energy cut, τ_zE and B_zE are not
    infrared safe at this order. A soft gluon alone in the current hemisphere
    gives a finite value. Varying DISENT's `-cutoff` from 1e-6 to 1e-10 moves
    the τ tail from −25% to +30%, while ρ_E stays unchanged.
  - So this is not a code bug. Future comparisons need E_cur > εQ, or the
    Q-normalised shapes.
- No bug in disorder found.
- One bug in the new comparison analysis, now fixed: at O(αs²), DISENT passes
  exactly-zero momenta in its collinear counterterms, which made `cmp_antikt`
  loop forever.
- Issues in NNLOJET and POWHEG are drafted as bug reports in the directory
  above.

### Speed (CPU time for 1e-4 on the total)

- disorder, inclusive: ~7 s.
- disorder, P2B NLO: ~120 s.
- POWHEG: LO 60–75 s, NLO 1800–4000 s.
- NNLOJET: LO 300–400 s.
- O(αs²) shapes:
  - disorder: 24 CPU-h for all shapes at 0.3–1% per bin.
  - NNLOJET: ~80 core-h per (process, shape) at 0.2–0.6%.
  - A full NC/CC set at 0.3% needs O(10⁴) core-h with NNLOJET.

## 2026-09-26 — MATTHR for γ/Z and W exchange (p2b at NLO for all processes)

Goal (AK): re-derive DISENT's photon three-parton matrix element, find the
form in which Z and W exchange fit into it while keeping photon exchange
bit-identical, then implement it for all NC/CC processes and validate
against NNLOJET and POWHEG. Only `MATTHR` (and the guards) are touched: at
NLO p2b it is the only DISENT matrix element that reaches the user routine.

### Derivation (FORM 4.3, `derivations/matthr/`)

- |M|² of l q → l q g, l q̄ → l q̄ g and l g → l q q̄ for each lepton
  helicity and quark chirality, with p3 and p1·p2 eliminated; sympy
  confirms exact identities A(l,h) = 32 (−q²) PAIR/(s13 s23)
  (s12 s13 for the gluon), where PAIR is (k·p1)² + (k'·p2)² for equal
  helicities and (k'·p1)² + (k·p2)² otherwise (swapped for antiquarks).
  For the gluon the pairs are (k·p3)² + (k'·p2)² and (k'·p3)² + (k·p2)².
- The helicity sum reproduces DISENT's `QQ`/`GQ`; with spin/colour
  averages, QQ = |M|²/(αs/2π) exactly.
- Hence M(i) = C2(i) QQ + C3(i) QQ3, where QQ3 has the numerator
  same − opposite and C2, C3 are the per-parton couplings of F2/x and F3
  in photon units, the same combination as at LO (Y₊ C2 + Y₋ C3).
  Antiquarks flip C3. In the gluon channel the parity-violating part is
  odd under 2 ↔ 3 and integrates to zero over DISENT's symmetric phase
  space (checked that GENTHR/GENDEC generate z and 1−z and the azimuth
  symmetrically), so it is dropped and the coupling is (C2(i)+C2(−i))/2.

### Implementation

- `parton_couplings` (`src/mod_matrix_element.f90`) mirrors the coupling
  and propagator logic of `eval_matrix_element_new` (γ, Z, γ/Z, Z-only,
  interference-only, CC, e±, ν/ν̄, HOPPET's quark couplings and complete
  generations for W exchange). `disent_couplings`
  (`src/mod_disent_interface.f90`, external, called from F77) adds the
  gluon couplings.
- `MATTHR`: `M(I) = C2(I)*QQ + C3(I)*QQ3`. For photon exchange C2 = EQ²,
  C3 = 0 and CG = EQ², so the arithmetic is the original one.
- Guards: Z/CC with p2b allowed up to NLO; beyond that still refused
  (VIRTHR, CONTHR, MATFOR, ... are photon-only, and they also call MATTHR).

### Checks

- `test_matrix_element` (13 processes): LO cross section from
  `parton_couplings` = quark-parton-model formulas (1e-5).
- `test_matthr` (13 processes): photon bitwise identical to the original
  expression; decomposition; QQ3/QQ → Y₋/Y₊ in the initial- and
  final-state collinear limits.
- Full ctest (59 tests) passes; the p2b validation outputs (photon,
  O(αs²), which also use MATTHR via VIRTHR/COLFOR/SUBFOR) are bitwise
  identical to `ref_runs` apart from the header stamp.
- The photon O(αs) p2b run of 2026-09-25 (1e8 events, seed 1), repeated
  with the new binary: histograms, grids and cross section byte-identical.
- Physics validation (`~/cernbox/disorder-comparisons`, stages
  `p2b_jets`, `p2b_shapes`): disorder `-nlo -p2b`, 4 × 1e8 events per
  process, against the O(αs) NNLOJET and POWHEG runs of 2026-09-25/26.
  All 10 processes, all histograms (Q², x, y, lab jet pT and η, τ_zE,
  B_zE, ρ_E) agree within statistics. Totals: χ² 397.8/360 (NNLOJET shapes),
  799.3/720 (POWHEG shapes), with no process or observable standing out.
  NNLOJET is not available for ν beams, nor for NC e⁺ (it computes e⁻).
- Negative control: the same CC e⁻ run with QQ3 dropped deviates by up to
  11% (τ_zE χ² ≈ 1.9e4/24, lab η_j 3.2e3/10 vs NNLOJET), so the
  comparison is sensitive to the new term.
- Photon ρ_E lowest bin: with seeds 1–4, disorder was 0.3% above both
  references (also in nc_gZ_e±, which share the seeds and the low-Q² photon
  events). Independent seeds 5–12 agree with NNLOJET/POWHEG, so this was a
  correlated ~2.5σ fluctuation in the lowest τ/B/ρ bins. Near the singular
  region the P2B weights have heavy tails, so the errors of the first bins
  may be slightly optimistic.
- Speed: old and new binary, pinned to separate cores, photon, 2e7 events:
  75.5 s and 75.0 s, with identical histograms.
- Full ctest (59 tests) passes from a clean clone of the branch.

## 2026-09-26 (cont.) — O(αs²) ingredients: photon re-derivation

Goal (AK): re-derive DISENT's photon O(αs²) ingredients in its own form,
validated point by point, then extend to Z/W. Tools in `derivations/o2/`:
a Fortran harness (`harness/`) that evaluates DISENT's photon functions at
DIS phase-space points, FORM programs generated by `gen_o2.py`, a Python
evaluator, an explicit Dirac-spinor check (`spinors.py`, `num4q.py`), and
a Python port of the BDK one-loop amplitudes from MCFM (`oneloop/`).

Photon results (exact, every point, arbitrary CF, CA):

| DISENT structure | FORM (helicity sum) × factor |
|---|---|
| q→qgg, HF(CF A + (CF−CA/2) B + CA C) | 1/64 [CF X1 + (CF−CA/2) X2] |
| g→qq̄g | 1/32 [same] |
| D1 = ERTD/LEID(4,−1,3,2)+(3,2,4,−1): boson on the incoming line 1→2, g→(3,4) | 1/32 |
| D2: boson on the (3,4) pair | 1/32 |
| E (identical-quark interference, all attachments) | −1/16 Re |
| CONTHR (quark / gluon Born) | ∓(CF/4)(4π/137)² 16π² /Q⁴ |

- ERTD(I,J,K,L) has the gluon producing (I,K) and the boson on (J,L).
- My first E comparison failed because I had not taken the real part.
  The independent explicit-spinor calculation reproduces DISENT's E
  exactly. This was my error, not DISENT's.

One loop:
- BDK/MCFM helicity amplitudes against 2 LEIV − Q²/2 ERTV.
- In e⁺e⁻ kinematics they agree exactly, up to tree-proportional terms,
  in both colour structures.
- In DIS kinematics, with all four helicity combinations summed, they also
  agree exactly: normalisation −1/(8π²). The tree terms are
  −7/2 + π²/2 − ½(L13² + L23²) (leading colour) and 7/2 + ½L12²
  (subleading colour), with L_ij = log(2|p_i·p_j|/Q²).
- **Correction:** I first summed only LL + LR, assuming RR = LL by parity,
  and found a mismatch in DIS that looked like missing π² terms in
  DISENT's continuation. That was wrong.
  - At one loop, absorptive parts give T-odd terms ∝ π ε(k,p1,p2,p3),
    which flip sign between LL and RR. They cancel only in the full sum.
  - MCFM's box function Ls₋₁ was checked against QCDLoop in all sign
    regions and agrees.
  - DISENT's virtual is correct.
- For Z/W, the T-odd terms multiply new coupling combinations. They are
  odd under reflection of the event and integrate to zero for
  reflection-symmetric observables over DISENT's symmetric phase space, so
  they are dropped.
- The parity-violating finite part is (LL + RR) − (LR + RL), with the same
  normalisation and tree terms. At tree level this combination reproduces
  MATTHR's QQ3 exactly.

### Implementation (Z/W at O(αs²))

- **MATFOR** takes the couplings from `disent_couplings4` (NC and CC
  separately) and computes
  `M(I) = C2N·Q + C3N·QNC3 + C2C·QCC + C3C·QCC3 + Σ_J CG(J)·QQ`.
  - QNC3 is the parity-violating part of Q: q→qgg, D1 and E.
  - For W exchange, E (identical quarks) is replaced by its attachment
    classes:
    - EXX: W on the incoming line in both amplitudes, with the pair
      flavour being the partner of the incoming quark; weighted with HF.
    - EXY: W on the incoming line interfering with W on the pair, e.g.
      d u ū from an incoming u; not identical particles, so no HF.
  - Dropped:
    - the parity-violating parts of the gluon-initiated and QQ (D2)
      terms, which are odd under q ↔ q̄ of the final pair;
    - the "different-line" interference, which is odd for vector
      couplings and for axial couplings ∝ Σ a_q′. That sum vanishes for
      complete generations; it is neglected as in the structure
      functions.
- **VIRTHR:** adds `C3·QQ3` with `QQ3 = −(4π/137)²·4·CF(CF·NX3 + CA·NY3)`
  from `VIRT3PV` (`src/disent_virt3.f`). Entries with |I| > NF keep their
  original values, since the scale-variation weights are normalised by V(I).
- **CONTHR3:** the parity-violating CONTHR,
  `−CF·4π²(4π/137)²·C3PV/(VV·Q⁴)`.
- **SUBFOR:**
  - s3, s9 use `C2·CONTHR + C3·CONTHR3`;
  - s7, s11 use CG;
  - s5 averages the quark- and antiquark-initiated Borns. Its local
    parity-violating part would otherwise not match the dropped one of
    the gluon-initiated ME, nor COLFOR, which sums quarks and antiquarks.
- **Photon:** all original expressions keep their operation order, and
  the new terms are only computed and added when non-zero.
- **Guards:** Z and CC with `-p2b` are allowed at NNLO.

### Checks

- Full ctest passes. The p2b validation outputs (photon, O(αs²)) are
  byte-identical to `ref_runs` apart from the header stamp.
- The generated tree functions reproduce DISENT's photon q→qgg, D1 and E
  to 1e-14, and CONTHR to 1e-13.
- The CC classes agree with the explicit-spinor calculation.
- The one-loop port reproduces ERTV/LEIV to 1e-13.
- New `test_subtraction` (13 processes): MATFOR/ΣSUBFOR → 1 in the
  final-state collinear (3∥4, 2∥3, 2∥4), soft (4, 2) and initial-state
  collinear (4∥1, 2∥1) limits, at λ = 1e-7 to 3e-3.
  - Correction to my first attempt: a direct comparison fails even for
    the photon. DISENT uses fixed labels for some structures (D2 is
    singular only at p2∥p1, while its dipoles sit at p3∥p1, p4∥p1) and
    relies on its permutation-symmetric phase space. Dipoles whose Born
    contains the unresolved parton are singular in these limits too; in
    a calculation the observable removes them.
  - After symmetrising over labels and leaving out unresolved-Born
    dipoles, the photon and all Z/W cases converge like √λ.
  - Negative control: dropping QNC3 gives ratios of 0.83.
- Speed: O(αs²) CC and Z-only runs are about 25% slower than photon.
- What the limit test can and cannot see (question from AK): removing the
  old SUBFOR fix (perm 4, s5: quark/gluon orientation of the Born)
  leaves all limit ratios unchanged. In the collinear limit that
  orientation error is equivalent to swapping partons 2 and 3, which the
  symmetrisation averages over. The bug only changes the finite part of
  the local dipole away from the limit, where it no longer matches its
  integrated counterpart. Such finite-part errors are only visible
  against an independent code (NNLOJET distributions).
- Local vs integrated consistency of the new spin-correlated dipole: the
  azimuthal average of CONTHR3 around the gluon equals half of MATTHR's
  QQ3 exactly, as CONTHR does for QQ. COLFOR/VIRTHR, which see only the
  average, therefore match the local dipoles.

### Validation against NNLOJET at O(αs²) (2026-09-26 evening)

Setup: disorder `-nnlo -p2b`, 30 × 5e7 events per process, compared with
NNLOJET 1.0.2 epLJJ/epNJJ/epNbJJ at NLO (O(αs²)).
- NNLOJET runs: 30 jobs per part; LO/V/R job lengths 30/30/100 min
  (the R jobs took about 3.3 h).
- Both use NNPDF40MC_nnlo_as_01180 and μ = Q.
- τ_zE and B_zE are taken with the current-hemisphere energy cut
  E_cur > Q/10: NNLOJET `dis_eventshapes = 0.1` in the PROCESS block,
  our `*_Ec` histograms. Without the cut these observables are not
  infrared safe at this order (see 2026-09-25). ρ_E needs no cut.

| process | observable | χ²/24 |
|---|---|---|
| NC γ e⁻ | τ_zE (cut) | 21.8 |
| NC γ e⁻ | B_zE (cut) | 32.6 |
| NC Z e⁻ | τ_zE (cut) | 27.3 |
| NC Z e⁻ | ρ_E | 36.8 |
| CC e⁻ | τ_zE (cut) | 35.2 |
| CC e⁻ | B_zE (cut) | 17.0 |
| CC e⁻ | ρ_E | 24.1 |
| CC e⁺ | τ_zE (cut) | 20.5 |
| CC e⁺ | ρ_E | 18.9 |

Total χ² = 234.2/216 (+0.9σ); no |pull| above 2.9. The photon τ/B
disagreement of 2026-09-25 is gone with the IR-safe definitions.

Further checks:
- Negative control, CC e⁻, same 30 seeds, with the one-loop
  parity-violating finite part (VIRT3PV) switched off:
  - The distributions shift by up to +5% (ρ, τ) and +16% (B) in the
    lowest bins; the correlated error of the difference is below 0.001%.
  - χ² against NNLOJET becomes 1568 (ρ), 2623 (τ) and 10854 (B) for 24
    bins.
  - So the comparison is highly sensitive to the new one-loop term. This
    also supports dropping the T-odd parts.
- Cutoff independence, CC e⁻: DISENT cutoff 1e-6 vs 1e-10 (5 × 2e7 each).
  - ρ, τ_Ec and B_Ec agree (χ² 28.6, 27.6, 15.3 for 24 bins).
  - The uncut τ and B change by 15–60% (χ² ≈ 390), as expected for
    observables that are not IR safe.
- Fluctuation check of the largest χ², NC Z e⁻ ρ_E (36.8/24, p ≈ 0.05;
  2026-09-27 night). Both codes were rerun with 30 new seeds (101–130),
  NNLOJET on the VEGAS grids of the original run.

  | comparison | χ²/24 |
  |---|---|
  | new NNLOJET vs new disorder | 20.3 |
  | new vs old disorder | 20.0 |
  | new vs old NNLOJET | 36.7 |
  | combined NNLOJET vs combined disorder | 27.0 |

  - The excess came from the NNLOJET sample, not from disorder: the two
    NNLOJET samples differ from each other at the same level (pulls of
    +2.7, +2.4 and −2.3 at ρ ≈ 0.13–0.22), while the two disorder samples
    agree.
  - NNLOJET's per-bin errors from the real-emission part may be somewhat
    small there, or this is a fluctuation; this check cannot tell which.
  - With the combined samples the comparison is 27.0/24.
  - Runs: `~/cernbox/disorder-comparisons/runs/topup_ncZ_rho`, comparison
    `tools/topup_compare.py`.

Not directly tested against an external code, because NNLOJET has no
e⁺ NC or ν beams: NC e⁺, NC ν/ν̄, CC ν/ν̄. Their O(αs²) matrix elements
use the same structures with different coupling values:
- e⁺ flips the sign of C3;
- ν/ν̄ change C2/C3.

These couplings are validated against the quark-parton model and, at
O(αs), against POWHEG. The CC runs cover both signs of C3 and the
CC-specific interference classes. `test_subtraction`/`test_matthr` run for
all of these flag sets.

## 2026-09-29 — p2b NC γ/Z and CC entries in the validation matrix

The validation matrix had p2b runs for photon exchange only, so a change
to DISENT's Z/W couplings would not have been caught by CI. Two entries
were added to `validation/configurations.txt` (NNLO, toy PDF, Q = 30,
x = 0.1, 7-point scale variation):
- `p2b_nc_includeZ_Q_30_x_0.1_positron_`: NC γ/Z, e⁺, so both the
  parity-conserving and parity-violating parts (with the e⁺ sign of C3);
- `p2b_cc_Q_30_x_0.1_`: CC e⁻.

Checks:
- The first attempt used Q = 50, whose Q² lies outside the
  `exclusive_lab_frame_analysis.f` window (25 < Q² < 1000). All 51 bins
  of those references were zero, so they tested nothing of DISENT.
  They were discarded before committing; the new entries fill 42 of 51
  bins, like the photon entry.
- Totals equal the inclusive NNLO results exactly (NC γ/Z e⁺ 1.225914,
  CC e⁻ 0.047969 pb/GeV² at Q = 30; also at Q = 50).
- The Z enters the jet histograms: the cut cross section is 0.44% above
  that of a photon-only e⁺ run with the same seed, which is about 400
  times the 1e-5 tolerance. The inclusive total is 0.46% above.
- Regenerating the references and re-validating them gives identical
  output (deviation 0). Full ctest: 71/71 pass (228 s). Each new run
  takes about 6 s.

The references are regression values. The physics validation of these
processes is the comparison with NNLOJET/POWHEG in the entries above.

README updated (same day):
- third-party code: the Z/W extension of disent (FORM-generated trees)
  and the MCFM 10.3 port in `src/disent_virt3.f` (GPLv3 or later; MCFM's
  copyright and SPDX notice added to that file's header);
- citation policy: Bern–Dixon–Kosower (hep-ph/9708239) for p2b at NNLO
  with Z or W exchange;
- a short description of the two modes;
- a Tests section (ctest), which the ctest commit had not added.

## 2026-09-29 — Pointwise cross-check against an independent DIS+jet implementation

AK has a separate, private implementation of the NLO 2+1-jet matrix
elements for all beams and NC/CC: a POWHEG-BOX process with VBFNLO-style
helicity amplitudes (Born e q → e q g, one loop, real e q → e q g g,
e q → e q Q Q̄, and crossings). Its code is not part of this repository, and
neither is the harness that calls it. The harness is linked against a
clean snapshot of this branch (4113237) and runs with disorder's own
flags. The other code is used unmodified, with the Z/W
widths set to zero as in disorder.

Checks, at 20 random points in DISENT's layout (`derivations/o2/harness/kin.py`,
Q = 12–100 GeV). Each check covers 11 processes: NC γ e±, γZ e±, Z e⁻, ν,
ν̄, and CC e±, ν, ν̄. Each covers every incoming parton (q, q̄, g).

- **Born, MATTHR.**
  - Independent reference: BDK/MCFM tree helicity amplitudes (V3SPIN/V3HEL
    of `src/disent_virt3.f`), crossed to each channel. Each lepton/quark-line
    chirality gets an SM coupling computed from sin²θ_W. This uses neither
    code's couplings.
  - The other code agrees with the reference to 1e-12 in every process and
    channel, once its lepton line is crossed explicitly (see the list at the
    end).
  - MATTHR agrees with both to 1e-13. For g-initiated channels this holds for
    the flavour sum symmetrised in q ↔ q̄: the dropped PV part is odd, and the
    unsymmetrised amplitudes agree with the reference.
  - For ν beams disorder is twice the other code, as it should be: one
    helicity is averaged over, not two.
- **One loop, VIRTHR's non-factorising part including VIRT3PV.**
  - The reference V/B is MCFM's assembly: N(A51 + A52/N²) + N/6 (the UV
    counterterm), converted from DRED to CDR by −C_F − C_A/6, exactly as
    POWHEG-BOX's Zj process does.
  - The other code's box-line virtual has the same loop functions to 1e-11
    (after one bug of its own; see the list at the end).
  - DISENT: (C2·QQ_nf + C3·QQ3)/MATTHR, minus its photon value at the same
    point, equals (1/3) × the same difference of the reference V/B, to 1e-10
    in all 11 processes. This requires averaging the reference over the
    reflected point (p_y → −p_y). Without it the relation fails by O(10).
    This confirms that what disorder drops is exactly the reflection-odd
    (T-odd, absorptive) part.
- **Real, MATFOR.** Compared after averaging over the 6 permutations of the
  outgoing partons 2, 3, 4, since MATFOR is not the pointwise |M|² for
  identical partons, only equal to it for label-symmetric observables.
  Processes: NC γ, γZ, Z, ν; CC e⁻, ν. For e⁺ and ν̄ the other code's real
  cannot be crossed.
  - Exact (≤ 7e-14):
    - q → q g g, including its PV part (QG3);
    - g → q q̄ g;
    - q → q Q Q̄ per pair flavour, as e_q²·TR·D1 + e_Q²·TR·D2 for photon
      exchange, and for Z/W summed over complete generations, including the
      PV part of D1 (FD13) and CG(Q)·D2;
    - the CC classes with the W on the incoming line and on the pair.
  - With DISENT's identical-quark interferences added to the other code (E
    for NC; EXX, EXY for CC), MATFOR agrees to 2e-14 (photon, CC, and NC Z
    with the b pair removed). The other code lacks these terms. They are up
    to ~1% (NC) and 2.4% (CC) of the matrix element at these points, and
    were checked earlier against explicit Dirac spinors (`num4q.py`).
  - The one remaining difference is the neglected axial interference of the
    boson on the two quark lines. It is ∝ a_q·a_Q and cancels over u+d and
    c+s, but not for the b pair at nf = 5. It is up to 0.4% (γZ e⁻) and 2%
    (pure Z, ν) of the four-parton matrix element at these points, and
    changes sign from point to point. As documented in MATFOR, disorder
    neglects it deliberately, as in the structure functions.

So every Z/W ingredient of disorder's p2b matrix elements up to O(αs²)
(MATTHR, VIRTHR incl. VIRT3PV, MATFOR incl. QNC3/QCC/QCC3) now also agrees
point by point with an independent implementation. The only exceptions are
the approximations documented in the code.

Issues found in the other code were reported to AK separately.

## 2026-09-29 — NaN in p2b CC at NNLO with scale variations; validation matrix for all p2b modes

Correction to the entry "p2b NC γ/Z and CC entries in the validation
matrix" above. The CC entry (`p2b_cc_Q_30_x_0.1_`) passed, but its
reference log had 1.5 million lines: 198,782 events were reported as NaN
and dropped by `analysis`. The reference had been generated with the
same NaN, so the comparison could not see it. I only noticed from the
size of the log before pushing.

Cause:
- With `-scaleuncert`, VIRTHR (and VIRTWO, COLTHR, COLFOR) normalise the
  three μ_F weights of each incoming parton to the central one.
- With W exchange some partons have no contribution at all (d, ū, s, c̄,
  b, b̄ for e⁻ CC): V(I) = 0 at all three scales, giving 0/0.
- In VIRTHR the q → g collinear term QG is zero by construction in MS-bar
  (KPFUNS is called with −X); that convolution is done in COLFOR. This
  was checked: all three weights vanish for these partons.
- So no contribution was lost apart from the NaN events themselves. The
  analysis drops those, so the central values of such runs were biased.
- Only NNLO runs with `-scaleuncert` are affected. The NNLOJET comparisons
  above ran without scale variations and are unaffected: their logs have
  no NaN.

Fix: `SCLNRM` (`src/libdisent.f`) leaves the ratio at 1 where the central
weight vanishes, since it multiplies a zero weight. If a varied weight is
non-zero there it warns once; that has not happened. For non-zero weights
the arithmetic is unchanged: all 19 existing fast validation entries
still pass against their old references, except the CC one.

Validation matrix:
- `run_validation.py` fails on any NaN or infinity in the output, also
  when generating references.
- New p2b entries (toy PDF, Q = 30, x = 0.1, fixed kinematics):
  - NC γ/Z at NLO;
  - Z only;
  - interference only (e⁺);
  - NC ν and ν̄;
  - CC at NLO;
  - CC without scale variations;
  - CC e⁺, ν and ν̄;
  - NC+CC.
  With the two earlier entries, p2b is covered for every NC/CC mode and
  beam, at NLO and NNLO, and on both scale-variation code paths.
- Checks of the new references:
  - every histogram is populated (42 of 51 bins);
  - every total and scale uncertainty equals the inclusive result;
  - the CC central histograms are bit-identical with and without
    `-scaleuncert`.
- Full ctest: 82/82.

The welcome line ("Welcome to disorder v. …") is now ignored by the
comparison, so that version changes do not require new reference logs.
Version set to 2.2.0 (welcome message and CMake project).

## 2026-09-29 — NNLOJET at O(αs²): NC γ/Z e⁻ with the interference

The table of 2026-09-26 has NC γ and NC Z separately. The full NC process
(γ + Z + γZ interference) was not in it. It is compared here for all three
shapes.

Setup (`~/cernbox/disorder-comparisons/runs/shapes_nlo_gz`, set up by
`tools/setup_shapes_nlo_gz.py`):
- disorder `-nnlo -p2b -NC -includeZ` at commit c1eff06 of this branch,
  with the comparison analysis: 30 seeds × 5e7 events.
- NNLOJET 1.0.2 epLJJ at NLO: 30 jobs per part, LO/V/R 30/30/100 min.
- HERA e⁻p (27.5 × 920 GeV), 1000 < Q² < 15000 GeV², 0.1 < y < 0.9,
  NNPDF40MC_nnlo_as_01180, μ_R = μ_F = Q, M_Z = 91.1876 GeV, M_W = 80.398
  GeV, zero widths.
- τ_zE and B_zE are taken with the current-hemisphere energy cut
  E_cur > Q/10, as before.
- The later commits change only the scale-variation ratio for zero weights
  (not used here), the help text and the version.

| observable | χ²/24, disorder γ/Z vs NNLOJET | max \|r−1\| | photon-only disorder vs NNLOJET γ/Z |
|---|---|---|---|
| τ_zE (cut) | 20.8 | 5.1e-3 | 8724 |
| B_zE (cut) | 24.2 | 4.7e-3 | 18726 |
| ρ_E | 16.8 | 1.9e-2 | 3402 |

- Total χ² = 61.8/72; no |pull| above 2.4.
- The control (photon-only disorder, 10 seeds) is 7–9% below NNLOJET, with
  pulls up to 50. So the comparison resolves the Z and γZ contributions
  with a large margin.
- Comparison: `tools/gz_compare.py`; bin-by-bin tables in
  `results/shapes_nlo_gz.txt`.

With this, disorder's p2b O(αs²) event shapes agree with NNLOJET for:
- NC γ, NC Z and NC γ/Z e⁻;
- CC e⁻ and e⁺.
