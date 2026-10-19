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

## 2026-09-30 — CI failure on main after the 2.2.0 merge

- CI failed on main at 9ee483c and 51aa33c (Test step; ctest exit code 8).
  The last green run was the scheduled one on 67fa107.
- Reproduced from a clean clone of main (51aa33c), configured as in CI
  (Release, FastJet, exclusive lab-frame analysis):
  - all unit and guard tests pass;
  - 31 of the 32 validation tests fail, each only on its `.log`, and each on
    the same line of the welcome banner: `Written by Alexander Karlberg
    (2023-2024)` in the references vs `(2023-2026)` since cc5b9f2 ("Updated
    banner"). Numbers are compared anywhere in a line, so 2024 vs 2026 is a
    relative deviation of 9.9e-4.
  - All result and histogram files agree to 0.0e+00.
- Cause: the years in the banner are volatile in the same way as the
  version on the line above it, which was already dropped.
- Fix: `run_validation.py` also drops the "Written by" line of the banner.
  The references are left as generated.
- Why it was not caught: CI runs only on pushes and PRs to main (and on the
  schedule), not on branches. The full ctest was not rerun after the banner
  commit went onto the branch before the merge.

## 2026-10-01 — DISENT dropped the O(αs²) part of 0.5% of its events (branch 2026-10-disent-xcut)

Found with the τ₂ slicing of NLO DIS 2+1 (branch 2026-10-tau-slicing, see its
notebook entries of the same day). AK: fix it on a new branch.

- **Mechanism.** `VIRTHR` samples the collinear momentum fraction X of the K
  and P terms (`GENCOL(2,…)`) and then did `IF (1-X.LT.CUTOFF) RETURN 1`.
  The alternate return ends the event, but the three-parton Born (NA = 1)
  has already been given to `USER`. So for those events the virtual, the
  collinear term, the real emission and its counter-events (all NA = 2) were
  never computed.
- **Rate.** `GENCOL` draws X = 1 − (1 − X_min) R^npow2 half of the time, so
  P(1 − X < CUTOFF) = CUTOFF^(1/npow2)/2 = 0.5% for npow2 = 4 and CUTOFF = 1e-8
  (disorder's defaults). Measured in the slicing runs: 0.49%.
- **Consequence.** DISENT's O(αs²) result was low by about 0.5% of the
  O(αs²) coefficient, a cutoff effect that falls only like CUTOFF^(1/4)
  (1.6% at 1e-6, 0.16% at 1e-10). In P2B the O(αs²) 2+1 events and their
  projections drop together, so the totals are unaffected and distributions
  lose about 0.5% of their O(αs²) 2+1 part. The cutoff test of 2026-09-26
  (1e-6 vs 1e-10) was not precise enough to see it.
- **Not affected:** `VIRTWO` has the same check on its own X (`GENCOL(1,…)`,
  npow1 = 2), but that X generates the 2+1 phase space and the abort comes
  before the Born; it only cuts the x_p → 1 corner (probability 5e-5).
- **Fix** (`src/libdisent.f`):
  - `VIRTHR` no longer aborts: for 1 − X < CUTOFF it sets XJAC = 0 and stores
    it (new entry `SETCOL` of `GENCOL`), and evaluates the X-sampled K and P
    terms at a safe X. All of them are proportional to XJAC, while the δ(1−x)
    terms and the virtual are kept. X can be exactly 1 (probability 5e-5), where
    the plus terms would give 0·∞.
  - `COLFOR` gets XJAC = 0, so its weight vanishes; the main loop no longer
    hands that zero-weight configuration to `USER` (its parton 4 can have zero
    momentum).
  - `GENFOR`, initial-state-spectator branch (which reuses X): `RETURN 1`
    when XJAC = 0. That dipole configuration is singular, like the ones its
    existing z cuts drop.
  - Remaining cutoff dependence: O(CUTOFF ln CUTOFF), independent of npow2,
    at no extra CPU cost.
- **Unit test** (`tests/test_disent_interface.f90`): every event whose Born
  is handed over must also get its O(αs²) virtual, and all weights must be
  finite. Negative control: with the old `libdisent.f` the new check fails in
  both settings (33 of 35 checks), with the fix all 35 pass.
- **Full ctest** (FastJet + exclusive analysis, g++ wrapper for the stale
  `~/.local/include/fastjet`): 69 of 82 pass, including all unit and guard
  tests, all inclusive validation runs and the two P2B NLO runs (no VIRTHR).
  The 13 failures are exactly the P2B NNLO runs, and only their histogram files
  (`disorder_*.dat`); every `xsct` total is unchanged. Previously aborted
  events now draw their remaining random numbers, so the event stream changes.
  Against the references χ²/n = 1.05 over 546 bins (central scale, both
  errors), max |pull| 2.6. The references are not regenerated; to be decided
  after review.

## 2026-10-01 — N-jettiness slicing for DIS: NLO 1+1 with tau_1^b (branch 2026-10-tau-slicing)

Context: the route to fully differential N3LO DIS (and VBF line by line) is
P2B with an NNLO DIS 2+1 calculation. AK and I agreed to try N-jettiness
slicing first (NLO 3+1 with dipoles above the cut, the SCET singular
cumulant below), keeping P2B + EFT matching (arXiv:2609.36007) as the later
upgrade. First step: the SCET ingredients at NLO, in the simplest setting.

`slicing/tau1b_nlo.py`: NLO DIS at fixed (x, Q^2, y), photon exchange,
mu_R = mu_F = Q, deterministic integration (scipy quad, rel. 1e-9 to 1e-11):
- tau_1^b of Kang, Lee, Stewart (arXiv:1303.6952): q_B = xP, q_J = q + xP
  (DIS thrust in the Breit frame). For two final partons
  tau = min(z, a(1-z)) + min(1-z, a z), a = (1-x_p)/x_p.
- Below the cut: the O(alpha_s) cumulant of KLS Eqs. (173)-(174) (hard
  function, quark jet function, hemisphere soft function, quark beam
  function with I_qq and I_qg). Checked by hand that its delta(1-z) terms
  are the sum of H, J, S and B (-9 - 2pi^2/3 - 6 ln tau - 4 ln^2 tau at
  mu = Q, in units of alpha_s CF/4pi, with the ln tau P_qq term).
- Above the cut: the O(alpha_s) real emission in (x_p, z_p), F2 and FL
  kernels, over the z intervals with tau > tau_cut (solved exactly; tau is
  piecewise linear in z).
- Exact NLO: the MS-bar coefficient functions C_q, C_g, C_Lq, C_Lg; they
  agree with hoppet's StrFctNLO (F2 to 1.5e-5, FL to 1e-6).

Two mistakes found on the way, both mine:
- The plus-distribution identity (1+z^2)[ln(1-z)/(1-z)]_+ =
  [(1+z^2) ln(1-z)/(1-z)]_+ + (7/4) delta(1-z): first coded with -7/4.
- The gluon F2 real kernel: first written from memory with 16 x_p(1-x_p),
  which gave a tau-independent offset of exactly 8 TR x_p(1-x_p) per flavour.
  The correct constant is 8 x_p(1-x_p) (no polynomial term in the
  transverse part; FL = 8 TR x_p(1-x_p)). This was fixed by requiring the
  tau -> 0 limit with the standard C_g and the beam-function constant
  2z(1-z), so it still needs an independent pointwise check against DISENT's
  MATTHR. The quark channel converged with the kernels as first written.

Result (sum - exact, relative to the O(alpha_s) correction):

| x, Q, y | tau_cut 1e-2 | 1e-3 | 1e-4 | 1e-5 | 1e-6 |
|---|---|---|---|---|---|
| 0.01, 20, 0.5 | 5.0e-2 | 6.1e-3 | 7.2e-4 | 8.3e-5 | 9.4e-6 |
| 0.001, 10, 0.9 | 2.8e-2 | 4.2e-3 | 5.6e-4 | 7.0e-5 | 8.4e-6 |
| 0.1, 50, 0.1 | -4.2e-2 | -1.0e-2 | -1.6e-3 | -2.2e-4 | -2.8e-5 |
| 0.4, 100, 0.3 | 1.6e-1 | 2.8e-2 | 4.0e-3 | 5.2e-4 | 6.4e-5 |

O(tau_cut ln tau_cut) convergence at all four points (y from 0.1 to 0.9,
so the FL separation is tested); FL alone converges linearly (FL is
power suppressed). Quark-only and gluon-only runs converge separately. The
same table relative to the Born is below 1.5e-6 at tau_cut = 1e-6.

Next: NLO 2+1 with tau_2 slicing against DISENT (one-loop hard function of
gamma* q -> q g, three-direction one-loop soft function, beam and jet
functions for both channels), and the pointwise MATTHR check of the real
kernels.

## 2026-10-01 — tau_2 slicing for DIS 2+1 at NLO: set-up (branch 2026-10-tau-slicing)

AK: go ahead with the 2+1 at NLO. Code: `slicing/mod_slicing_scet.f90`,
`slicing/tau2_nlo.f90`, design and formulas in `slicing/README.md`
(GSTW conventions, arXiv:1505.04794; geometric measure, Q_i = 2 E_i, Breit
frame). Ingredients and how they were checked:
- jet and beam functions (quark and gluon; gluon beam coefficients from
  arXiv:1405.1044) and the three-direction soft function (JSTW
  arXiv:1102.4344 non-hemisphere integrals I0, I1, as 1D integrals with a
  complex dilogarithm);
- hard function from DISENT: VIRTHR's constant (factorising virtual + CS
  I operator) + the finite part of the CS I operator + (pi^2/12) sum C_i,
  plus DISENT's non-factorising virtual. The conversion was checked on the
  two-parton case: VIRTWO's CF(2 - pi^2) gives the DIS form factor
  CF(-8 + pi^2/6). For three partons the single logs cancel analytically.
- `test_scet`: complex dilog vs mpmath; I0, I1 vs mpmath at 30 digits
  (1e-11; a plain scipy 2D quadrature, first used as the reference, was off
  by 3e-5 and failed where the log point is inside the region); the full
  1+1 cumulant (beam, jet, two-direction soft, hard) vs the validated
  Python tau_1^a cumulant (5e-7).
- `tau2_nlo`: in one DISENT run, the reference (all O(alpha_s^2) pieces),
  the real events above tau_cut (counter-events and collinear terms have
  tau_2 = 0 exactly, checked in the run), and the Born events reweighted
  with the cumulant. Observable: tau_zQ (= tau_1^b) bins above 0.05.
  Speed about 0.4 ms per DISENT event including the reweighting.
- First production (x = 0.01, Q^2 = 400 GeV^2, s = 101200 GeV^2,
  NNPDF30_nlo_as_0118; 26 x 2M events, `slicing-runs/nlo21-x0.01-Q400`):
  (sum - DISENT)/DISENT for tau_zQ in [0.05, 1) went -0.47, -0.36, -0.28,
  -0.25, -0.24 for tau_cut = 1e-2 ... 5e-4 (errors 1%): flattening, but not
  to zero. Quark-only and gluon-only runs (`-pdfmask`) were both off. With
  the O(alpha_s) Born rates per bin, the offset per Born (alpha_s/2pi units)
  was about -1 for tau_zQ < 0.5, but -27 to -44 in tau_zQ in [0.5, 1), with a
  slope in ln tau_cut (single-log mismatch there, while the ln^2 coefficients
  agree: below -11.3, above +11.0 per Born).
- Cause: my observable, not the slicing. A Born with an empty current
  hemisphere has tau_zQ = 1 exactly (a finite region of the 2+1 phase
  space), and a soft gluon into the current hemisphere moves it to 1 - eps.
  A bin with upper edge tau_zQ = 1 is therefore not IR safe; this is why the
  earlier shape validations cut on E_cur > Q/10. The first runs' [0.5, 1)
  and [0.05, 1) bins are invalid. In the IR-safe bins (tau_zQ < 0.5) the
  offset per Born was -4.4, -2.6, -1.4, -1.1, -1.2 (+- 0.25) for tau_cut =
  1e-2 ... 5e-4.
- On the way, the hard function was cross-checked against the CS paper's
  explicit e+e- -> 3 jets one-loop matrix element and I operator
  (CS (D.16)-(D.18)): DISENT's QQ is their V + I constant with the DIS
  analytic continuation (-CF pi^2: only the outgoing pair is timelike), my
  I0 reproduces their I-operator finite part term by term, and the same
  conversion gives the timelike form factor CF(-8 + 7 pi^2/6).
- Rerun with IR-safe bins (tau_zQ in [0.05, 0.5), five bins), tau_cut 2e-2
  ... 1e-4, 106 x 2M events (`slicing-runs/nlo21v2-x0.01-Q400`).
- IR-safe rerun (106 x 2M events, `nlo21v2-x0.01-Q400`): the offset per
  Born (alpha_s/2pi) in tau_zQ in [0.05, 0.5) is -6.8, -4.3, -2.5, -1.3,
  -0.94, -0.79 +- 0.06, -0.91 +- 0.09, -0.87 +- 0.11 for tau_cut = 2e-2 ...
  1e-4: a plateau at about -0.85, i.e. -7 to -8% of the NLO coefficient,
  roughly uniform over the bins. Not power corrections.
- By channel (30 x 2M each, `-pdfmask`): quark only -0.55 +- 0.1, gluon
  only -1.27 +- 0.12; x = 0.1, Q^2 = 1000 (quark dominated) -0.56 +- 0.1;
  x = 0.001, Q^2 = 50 (gluon dominated) -1.1 +- 0.2. So about -0.42 C_a
  per Born, C_a the Casimir of the incoming parton (close to -pi^2/24 C_a;
  ratio g/q 2.3 +- 0.5, also compatible with 2 = number of quark jets).
- Checked and excluded since: the 1+1 test with the geometric-measure
  tau_1 (axes by minimisation, the measure of the 2+1 code; KLS (173)
  cumulant) converges to the exact NLO at three points (1e-5 relative at
  tau_cut = 1e-6), so beam, jet, two-direction soft and hard pieces are
  right in 1+1; the three-direction soft I-terms equal a direct angular
  integral of the eikonal times ln(m_hemi/m_true) (4.780559 both); the
  soft integrals agree with mpmath at extreme ratios; the hard function
  matches CS (D.16)-(D.18) and the incoming-parton I operator of CS (8.25)
  (DISENT's QQ contains only the symmetric operator; I^(1), I^(2) are in
  its sampled K+P terms).
- Open. Next diagnostic: the measure in a frame boosted along z (the beam
  piece is boost invariant, jet and soft are not); runs with Y = +1 and
  Y = -1 (the latter on idle th desktops, nice 19, AK's rule).

## 2026-10-01 (afternoon) — tau_2 slicing 2+1: boosted frames, more checks, colour decomposition

Continuing the open -0.42 C_a offset per Born (alpha_s/2pi) of the NLO 2+1
slicing against DISENT (x = 0.01, Q^2 = 400 GeV^2, tau_zQ in [0.05, 0.5)).

- Measure in a frame boosted along z (`-boostY`), per-Born offset at
  tau_cut = 5e-4 / 2e-4 / 1e-4:
  Y = 0 (106 runs) -0.78 +- 0.06 / -0.89 +- 0.08 / -0.88 +- 0.10;
  Y = +1 (30 runs) -1.20 +- 0.09 / -1.30 +- 0.13 / -1.28 +- 0.16;
  Y = -1 (28 runs) -0.86 +- 0.13 / -0.64 +- 0.18 / -0.57 +- 0.23.
  Y = +1 is about 3 sigma more negative, Y = -1 slightly less. Not
  conclusive: the power corrections at Y = +1 are about those of Y = 0 at
  twice the cut, so slow convergence cannot be excluded. Settling it at
  tau_cut = 1e-5 would need O(10^3) runs (the error there is about twice
  that at 1e-4), so I moved to cheaper diagnostics.
- 1+1 with the geometric measure in a boosted frame (`tau1b_nlo.py --geo
  --boostY`, lambda_B = e^Y tau, lambda_J = e^-Y tau; the cumulant changes
  by 2 C_F Y^2 q + Y (P x f)): converges at Y = +1 and -1 (remainder
  3e-6 and -2e-5 of the Born at tau_cut = 1e-5). So the energy
  dependence of beam vs jet logs is right, which matters in 2+1 where
  E_a = Q/(2 x_p).
- Soft I-terms, independent check with a fresh derivation of the
  normalisation from the eikonal (dP = (alpha_s/2pi) sum_{i/=j} (-T_i.T_j)
  dE/E dOmega/2pi W_ij): the third-direction correction is
  sum_{i/=j} (-T_i.T_j) K_ij with K_ij = int dOmega/(2pi) W_ij ln(m_ij/m_true).
  The code's (2 sum_m I_ij,m + 2 sum_m I_ji,m)/2 equals -K_ij at four
  configurations to 1e-4 (quadrature-limited). The individual I_ij are
  not symmetric in i <-> j; only the dipole sum enters.
- Hard function, independent check: with the BDK/MCFM remainder of
  2026-09 (tree terms N(-7/2 + pi^2/2 - (L13^2 + L23^2)/2) +
  (1/N)(7/2 + L12^2/2), plus N/6 and the DRED -> CDR shift -C_F - C_A/6),
  MCFM's finite part is N(-4 + pi^2/2 - (l13^2 + l23^2)/2) + (1/N)(4 + l12^2/2).
  That is exactly QQ + I0_CS (worked out analytically: the single logs and
  n_f cancel). So the quark-channel hard function is confirmed in DIS
  kinematics, independent of DISENT's QQ and of my I-operator conversion.
- **Correction of an assumption I was about to make:** I considered DISENT's
  NLO 2+1 itself as the culprit, and looked at NLOJET++ for an independent
  check. That is not needed: DISENT's O(alpha_s^2) event shapes
  (P2B -nnlo) were validated against NNLOJET on 2026-09-26 (tau_zE with
  E_cur > Q/10, B_zE, rho_E; chi^2 234/216, errors 0.2-0.6%). An error of
  2.4% of the 2+1 Born would have shown there. The known gluon-channel
  DISENT bug (Borsa et al., arXiv:2010.07354 App. A) is fixed in our copy
  (src/libdisent.f:1760). So the reference is trusted and the problem is
  on the slicing side.
- Bookkeeping checked against DISENT's event loop: COLFOR puts parton 4
  along the beam (tau_2 = 0, observable unchanged); counter-events have
  tau_2 = 0; the Born and the virtual share the phase-space weight.
- Colour decomposition (AK suggested leading colour first): colour factors
  are now run-time options (`-CF -CA -TR`, passed to DISENT as well;
  diagnostics only). C_F = 3/2 is refused by DISENT's range check, so
  leading colour is C_F = 4/3, C_A = 8/3. For C_A = 0 the I-operator term
  T_g.T_k gamma_g/C_g is a removable 0/0 (limit -gamma_g/2), now coded;
  for T_R = 0 the gluon Born vanishes (guarded). Quark channel
  (`-pdfmask 1`), 30 x 2M each: leading colour, leading colour with
  T_R = 0, C_A = 0, T_R = 0, C_A = T_R = 0, and QCD
  (`slicing-runs/col-*`), on thserv23/24 (nice 10) and 11 idle th desktops
  (nobody logged in, nice 19, half the cores).
- Cut list extended to 3e-5 and 1e-5 (nt = 10).

### Later on 2026-10-01: corrections, and DISENT's aborted events

- **Correction (job lists):** in the first colour-decomposition batch, the
  QCD lines of the job lists ended in a blank, and `xargs -L 1` treats a
  trailing blank as a line continuation. So 22 of the 30 "QCD" runs swallowed
  the next leading-colour line and actually ran with C_A = 8/3 (with the QCD
  seeds), and those 22 leading-colour jobs never ran. Regrouped by the colour
  factors printed in each run.log: 30 leading-colour runs and 8 QCD runs. The
  numbers I first quoted for "QCD" were mostly leading colour.
- Regrouped fit, quark channel, o = (4/3) a + b C_A + c T_R n_f per Born at
  tau_cut = 2e-4 / 1e-4: a = -0.47 +- 0.08 / -0.47 +- 0.10 (abelian C_F^2,
  -0.63 per Born), b = -0.01 +- 0.04 / -0.04 +- 0.06, c = -0.11 +- 0.05 /
  -0.07 +- 0.06. Mostly abelian.
- **Correction (plateau):** with cuts down to 1e-5 the offset is not
  constant; all configurations drift further negative below 1e-4 (about -0.8
  per Born between 1e-4 and 1e-5).
- Forward-jet power corrections: refuted. An IR-safe cut (`-smin`) on both
  jets of the minimising 2-jettiness partition, s_1J > 0.003 or 0.03
  (removing 2% or 10% of the Born rate), leaves the offset unchanged
  (-0.80 +- 0.11 and -0.88 +- 0.12 at 5e-4).
- Hard function, complete check: the full one-loop V/B used by tau2_nlo
  (QQ + I0_CS + DISENT's non-factorising part / Born, normalisation
  included) equals MCFM's BDK assembly (N A51 + A52/N + N/6 - C_F - C_A/6,
  averaged over the reflected point) at 40 DIS points to 2e-13.
- **Cause found (candidate):** DISENT's VIRTHR samples the collinear X and
  executes `IF (1-X.LT.CUTOFF) RETURN 1`, which aborts the rest of the
  event (alternate return to label 1000). By then the three-parton Born
  (NA = 1) has already been passed to USER. GENCOL draws
  X = 1 - (1 - X_min) R^npow2 half the time, so
  P(1 - X < CUTOFF) = (1/2) CUTOFF^(1/npow2) = 0.5% for npow2 = 4,
  CUTOFF = 1e-8 (the defaults of disorder as well). Measured: 335 of 67943
  Born events in the bins (0.49%) have no O(alpha_s^2) part.
  - In tau2_nlo those events still received the below-cut weight (about
    -320 per Born at 1e-4, growing like ln^2 tau_cut), but no
    above-cut, virtual or reference. tau2_nlo now drops their Born and
    below-cut weights (`-keepall` restores the old behaviour) and prints the
    count. Paired reruns (same seeds) of Y = 0, Y = +1 and the abelian case
    are running.
  - Note on the size: p times the below-cut per Born is -1.6 at 1e-4, more
    than the observed -0.88, so this alone may over-correct; the paired
    reruns will tell.
  - Consequence for DISENT itself: its O(alpha_s^2) result is missing the
    whole NA = 2 part (virtual, collinear, real, counter-events) of a
    fraction p of the events, i.e. it is low by about p = 0.5% of the
    O(alpha_s^2) coefficient, a cutoff effect scaling like CUTOFF^(1/npow2).
    In disorder's P2B mode the same applies to the O(alpha_s^2) 2+1
    events (and to their projections, so totals are unaffected but
    distributions are). To report to AK; src/ not touched.
- Paired reruns with the fix (same seeds as the old runs): new - old per Born
  = +0.33, +0.75, +0.98, +1.31, +1.60 at tau_cut = 5e-3, 1e-3, 5e-4, 2e-4,
  1e-4, i.e. exactly the below-cut weight of the aborted events (0.49%).
  **Correction:** the abort fully explains the earlier "plateau": with the
  fix the offset is positive and still moves with tau_cut. QCD, all
  channels (89 runs): +0.22, +0.49, +0.81, +1.12, +1.14 at 5e-4, 2e-4,
  1e-4, 3e-5, 1e-5.
- Decomposition with the fix (quark channel unless stated), slope of the
  offset in ln(1/tau_cut) between 2e-4 and 1e-5: QCD all 0.32 +- 0.07,
  quark 0.34 +- 0.11, gluon 0.45 +- 0.13, abelian (90 runs) 0.17 +- 0.04,
  C_A = 0 (90 runs) 0.23 +- 0.05, T_R = 0 0.10 +- 0.15, leading colour
  0.18 +- 0.11. The abelian part is significant; n_f and C_A parts are
  consistent with zero.
- 1+1 test of the soft L ln s_hat terms: the geometric measure with
  rescaled directions, min(rho_B n_B.k, rho_J n_J.k), is a dipole with
  s_hat = rho_B rho_J (Lorentz covariance); `tau1b_nlo.py --geo --rho`.
  Converges for rho = (2,2), (0.5,0.5), (3,0.7) (remainder 3e-7 to 1e-6 of
  the Born at tau_cut = 1e-6, linear in tau_cut). So the dipole's L ln s_hat
  and constant terms are right; with the soft I-terms, hard function,
  beam and jet already checked, every abelian ingredient is verified.
- The high-statistics sets flatten below 3e-5 (QCD 1.12 -> 1.14, abelian
  0.56 -> 0.50), and the pure-log fit for QCD is poor (chi^2 6.3/2). Two
  readings are open: a single-log mismatch, or convergence to a constant
  (about +1.2 per Born in QCD, +0.55 abelian) with slow power corrections
  (a sqrt(tau_cut) form fits). Not resolved.

### 1 Oct, evening: DISENT fix on its own branch, invariant measure, literature

- DISENT's abort fixed on branch `2026-10-disent-xcut` (worktree
  `~/work/disorder-disent-xcut`, see its notebook): the fixed DISENT is
  cutoff independent within 0.3% (1e-10 vs 1e-6: −0.41 ± 0.27%), the unfixed
  one shows the predicted CUTOFF^(1/4) dependence (+1.18 ± 0.26%, expected
  +1.44%). New slicing runs use the fixed DISENT (`tau2_nlo_fixed2`), where
  the driver's own drop of aborted events is a no-op (0 aborted events).
- Frame-independent measure (`-measure inv`, AK agreed after reading
  arXiv:2408.05265): Q_i = Q for all regions,
  τ₂ = min(min_j 2η_j P·p_j, min_{j<k} 2p_j·p_k)/Q², η_j = x(1 + s_kl/Q²),
  cumulant at λ_B = λ_J = τ_cut, soft function of the Born invariants
  ŝ_ij = 2p_i·p_j/Q² (`soft_from_s`). Same output as before for the
  geometric default with the same seeds.
- Results, 90 × 2M each, offset per Born at τ_cut = 2e-4 / 1e-4 / 3e-5 / 1e-5:
  - invariant, QCD: −0.17 ± 0.08 / −0.13 ± 0.09 / −0.06 ± 0.14 / +0.11 ± 0.24,
    slope 0.07 ± 0.07: converges (from below, as power corrections), consistent
    with zero, i.e. within about 1% of the NLO coefficient (9.9 per Born);
  - invariant, abelian quark channel: −0.42 ± 0.08 / −0.24 ± 0.09 / −0.30 ± 0.14
    / −0.30 ± 0.18, slope 0.08 ± 0.05; similar in all τ_zQ bins; open (1.7σ at
    1e-5);
  - geometric (Breit frame), QCD: +0.48 / +0.80 / +1.09 / +1.14 ± 0.24, still
    rising (slope 0.31 ± 0.07): unexplained; the invariant measure is now the
    working choice. Decomposition of the invariant measure by channel and
    colour is running.
- Literature (AK asked): no DIS 2+1 at NNLO by slicing exists. NNLO DIS dijets
  only from NNLOJET (antenna; 1606.03991, 1703.05977), N3LO DIS jets by P2B on
  NNLOJET (1803.09973). Slicing with jets: H+1 jet in MCFM (1906.01020: NLO
  converges to the dipole result like ε ln ε, 0.15%, with the geometric
  measure in the rest frame of the Born system; the hadronic-frame definition
  converges much more slowly), V+jet NLP (1907.12213: √(T₁/Q) terms cancel;
  frame of the Born system recommended), DIS single-inclusive jet at NNLO
  with τ₁ᵃ (1607.04921). P2B-improved slicing (2408.05265) removes fiducial
  √τ corrections; same decomposition as disorder's P2B.
- k_T-like slicing (y₂₃; Buonocore, Grazzini, Guadagni 2508.19226,
  2512.03954): NNLO only for two-jet final states (quark jet function, two-
  direction soft function, E-scheme factorisation-breaking term); "to compute
  NNLO corrections to e+e− → 3 jets, the gluon jet function is needed". Not
  usable for DIS 2+1 at NNLO yet.
- τ₂ at NNLO: all ingredients exist: two-loop γ*qg amplitudes (Garland et al.
  hep-ph/0112081, hep-ph/0206067, crossed by Gehrmann–Remiddi), jet functions
  (Becher–Neubert hep-ph/0603140, Becher–Bell 1008.1936), beam functions
  (Gaunt–Stahlhofen–Tackmann 1401.5478, 1405.1044), the three-direction NNLO
  soft function (Campbell et al. 1711.09984, generic kinematics; Bell et al.
  2312.11626). Crossing (AK asked): Bell et al. §2.4: the process dependence
  (incoming vs outgoing) sits only in the real-virtual tripoles, weighted by
  λ_AB, and "in processes with three hard partons, the sum of the tripole
  contributions vanishes because of colour conservation"; the dipole part has
  a universal cos(πε). So the pp → V+jet soft function applies to DIS 2+1
  (qq̄g: one-dimensional colour space, tripoles vanish identically).
- Invariant-measure decomposition (fixed DISENT, 90 × 2M each; settings
  checked in every run.log), offset per Born at τ_cut = 1e-4 / 3e-5 / 1e-5,
  slope over 2e-4 … 1e-5:
  - QCD, all channels: −0.13 ± 0.09 / −0.06 ± 0.14 / +0.11 ± 0.24, 0.07 ± 0.07
  - QCD, quark channel: +0.16 ± 0.10 / +0.20 ± 0.18 / +0.02 ± 0.25, 0.00 ± 0.07
  - QCD, gluon channel: −0.22 ± 0.16 / +0.13 ± 0.24 / +0.31 ± 0.38, 0.22 ± 0.10
  - quark, T_R = 0: +0.03 ± 0.11 / −0.06 ± 0.17 / −0.23 ± 0.22, −0.02 ± 0.06
  - quark, C_A = 0: −0.18 ± 0.08 / 0.00 ± 0.12 / −0.10 ± 0.17, 0.15 ± 0.05
  - quark, abelian: −0.24 ± 0.09 / −0.30 ± 0.14 / −0.30 ± 0.18, 0.08 ± 0.05
  At 3e-5 and 1e-5 all six are consistent with zero (χ² 6.4/6 and 5.1/6);
  largest deviation the abelian part (−2.1σ, −1.7σ). So NLO 2+1 with the
  invariant measure reproduces DISENT within the statistical precision, about
  1.4–3.8% of the NLO coefficient per configuration. Still open: the
  geometric measure's rising residual (unexplained; it must not hide
  something the invariant measure shares) and more statistics for the
  abelian part.

## 2026-10-02 (night) — the geometric residual was the driver drop, not the measure

**Result.** With the fixed DISENT (branch `2026-10-disent-xcut`, zero
Jacobian instead of the abort), the geometric (Breit-frame) measure converges
like the invariant one. The "rising geometric residual" came only from runs
that used the unfixed DISENT with the driver-side drop of aborted events.
Offset per Born (QCD, all channels, x = 0.01, Q² = 400; scratchpad `edge.py`),
τ_cut = 1e-4 / 3e-5 / 1e-5, slope per ln(1/τ_cut) over 5e-4 … 1e-5 (points
treated as independent, indicative):

| geometric measure | aborted | 1e-4 | 3e-5 | 1e-5 | slope |
|---|---|---|---|---|---|
| unfixed + drop, CUTOFF 1e-6 (`cut-old-1e-6`, 100 runs) | 1.6% | +1.94 ± 0.11 | +2.60 ± 0.16 | +3.06 ± 0.22 | 0.58 ± 0.04 |
| unfixed + drop, CUTOFF 1e-8 (`fix-qcd`, 90) | 0.5% | +0.80 ± 0.11 | +1.09 ± 0.17 | +1.14 ± 0.23 | 0.29 ± 0.04 |
| unfixed + drop, CUTOFF 1e-10 (`cut-old-1e-10`, 100) | 0.16% | +0.17 ± 0.11 | +0.20 ± 0.16 | +0.58 ± 0.22 | 0.17 ± 0.04 |
| fixed, CUTOFF 1e-6 (`cut-new-1e-6`, 100) | 0 | +0.02 ± 0.11 | +0.12 ± 0.16 | −0.32 ± 0.22 | 0.06 ± 0.04 |
| fixed, CUTOFF 1e-8 (`cut-new-1e-8`, 100) | 0 | −0.00 ± 0.11 | +0.03 ± 0.16 | −0.10 ± 0.22 | 0.06 ± 0.04 |
| fixed, CUTOFF 1e-10 (`cut-new-1e-10`, 100) | 0 | −0.02 ± 0.11 | −0.07 ± 0.16 | +0.17 ± 0.22 | 0.08 ± 0.04 |
| fixed, Born with min s_ij(geo) < 0.01 removed (`edge-geo-s0.01`, 90) | 0 | −0.05 ± 0.12 | −0.22 ± 0.17 | −0.10 ± 0.24 | 0.07 ± 0.05 |
| invariant, fixed (`inv-qcd`, 90; for comparison) | 0 | −0.13 ± 0.11 | −0.06 ± 0.16 | +0.11 ± 0.22 | 0.12 ± 0.04 |

The bias scales with the abort rate and vanishes with the fix at every
cutoff. All `inv-*` sets used the fixed DISENT (0 aborted events) and all
geometric `fix-*` sets the unfixed one with the drop (e.g. 310591 aborted
events in `fix-qcd`), so yesterday's measure comparison was confounded with the
DISENT version. (The `cut-new` sets were made for the DISENT cutoff validation,
and I had only read their O(αs²) coefficients.)

**Mechanism** (`main:src/libdisent.f`). GENCOL stores the X it draws in
VIRTHR; COLFOR and the initial-state-spectator branch of GENFOR
(`OMIT.EQ.1`) take it back with GETCOL and weight with its full-range
Jacobian. These X-dependent pieces (K+P, COLFOR, IS-spectator reals with their
dipoles) are therefore unbiased over X < 1 − CUTOFF and are not thinned by the
abort, while the X-independent ones (virtual, final-state-spectator reals) are
lost with probability p. Neither treatment of the slicing's Born-level weights
matches this: dropping them (the driver drop) leaves −p × [IS-branch real
below τ_cut + its dipoles + K+P + COLFOR], whose real-below-cut part grows like
ln²(1/τ_cut); keeping them (`-keepall`) leaves +p × [virtual + FS-branch
reals below cut]. Only the fixed DISENT is consistent.

**Corrections of the record.**
- 1 Oct ("Later on 2026-10-01"): "tau2_nlo now drops their Born and
  below-cut weights" and the "paired reruns with the fix" / "decomposition
  with the fix" refer to the driver drop, which does not remove the bias. The
  slopes quoted there (QCD 0.32, quark 0.34, gluon 0.45, abelian 0.17, C_A = 0
  0.23, …) are driver-drop artefacts, as is the reading "a single-log mismatch
  or slow power corrections".
- 1 Oct evening: "geometric … still rising (slope 0.31 ± 0.07): unexplained"
  is explained above; the geometric measure does not have a residual. The
  invariant measure remains the working choice for its frame independence,
  not because the geometric one fails.

**Edge test** (fixed DISENT, 90 × 2M each, IR-safe cut on the Born: all pairs
beam/jet/jet at s_ij > s_min in the Breit frame, `-smin`, applied to DISENT and
slicing alike). Born fraction with geometric min s_ij < 1e-2 / 1e-3 / 1e-4
before the cut: printed per run (`Born fraction …` line); with s_min = 1e-3,
3.8% of the remaining Born weight has s_ij < 1e-2. Offsets at 1e-4 / 3e-5 /
1e-5: s_min = 0.01 +0.11 / −0.04 / +0.19 (± 0.14 / 0.21 / 0.29), s_min = 0.001
−0.00 / −0.17 / +0.07; invariant with s_min = 0.01 +0.03 / −0.10 / −0.51
(± 0.13 / 0.20 / 0.28). All consistent with zero, as is the uncut fixed run,
so the near-collinear Born configurations do not produce a visible NLO
effect at these τ_cut.

**Soft function near the edge** (AK, after talking to Rudi Rahn: a single log
in the soft function when two legs become collinear). Bell et al. 2312.11626
§4: when a third direction a approaches leg i of the dipole (i,j)
(n_ai = 2δ → 0), gluons at angle ~√δ around n_i (the "correction region")
shift the renormalised NLO dipole coefficient by −π²/3 (finite, approached
like √δ); the dipole spanned by i and a tends to the 0-jettiness value. At NNLO
the correction diverges: c⁽²⁾ = (4π²/9) T_F n_f ln(2δ/n_ij) − (11π²/9) C_A ln(2δ/n_ij)
+ const, i.e. −(π²/3) β₀ ln(2δ/n_ij): the running-coupling dressing of the NLO
offset. That is the single log. Check of our NLO soft function (GSTW form,
`soft_from_s`, scratchpad `softlim/softlim.f90`): for m → i the extra term
4[I₀ ln(s_jm/s_ij) + I₁] of g(i,j) tends to −6.57975 = −2π²/3, i.e. −π²/3 for
each ordering of the pair as in Bell et al. (constants differ from theirs only
by a universal Laplace-space shift), approached like √δ (0.049, 0.016, 0.005,
0.0016 at δ = 1e-4 … 1e-7), stable to δ = 1e-9, for back-to-back and generic
(i,j); g(j,i) changes by O(δ); g(i,m) = ln² s_im − ζ₂ + O(δ) (0-jettiness).
So the NLO soft function is right and numerically stable in this limit. For
NNLO the log matters for Born configurations with nearly collinear
directions; the invariant measure avoids small s_ij for resolved jets.

**arXiv:2604.13167** (Buonocore, Delto, Melnikov, Monni, Pikelner, Vita), read
in full. Dipole part of the N-jettiness soft function = inclusive soft function
(T_N evaluated on the total soft momentum; analytic from the fully differential
soft function, known to three loops) + Δ, with Δ = 0 at NLO, a finite
five-dimensional tree-level integral at NNLO (eq. 4.25: the double-soft
maximally non-abelian matrix element times ln[T_N(k̂₁+ξk̂₂)/T_N(k̂₁,ξk̂₂)]),
and NLO-like at N3LO. The geometry enters only through c_Ai = n_A·n_i/n_i·n_j,
c_Aj and the azimuth of n_A (eqs. 4.3–4.5); the inclusive part is
(n_ij/2)^{2ε} × the 0-jettiness one + a finite three-dimensional integral of
t^{4ε} − t₀^{4ε} (eq. 4.18). Tripoles: compact two-dimensional formula
(4.43–4.45), not needed for DIS 2+1 (three coloured legs, tripoles vanish by
colour conservation through N3LO, as for pp → V+j). Timing quoted: 0.1–1 s per
dipole at 1%, 1–10 s at 0.1%. No discussion of collinear hard directions; in
their variables Bell's correction region is z ~ c_Ai, so a numerical
implementation must sample z logarithmically when c_Ai is small. For DIS 2+1
with the invariant measure the soft function depends on the three ŝ_ij; with
the overall homogeneity that is a two-dimensional grid of the three dipoles,
which their method makes cheap. This is the route for the NNLO soft function.

**VBF (proVBFH notes for details).** Fixed-scale coefficient of the old
code's missing initial-state region: B(≥3 jets) = (2.63 ± 0.01)e-4 pb per
e-fold of k_T² (proVBFH-cs `cs_estimate 2`, `runs/estimate-coll-fixmh`).

**Per channel, geometric, fixed DISENT** (`geofix-qQCD`, `geofix-gQCD`, 45 × 2M
each), offset at 1e-4 / 3e-5 / 1e-5: quark +0.21 ± 0.17 / +0.24 ± 0.25 /
+0.12 ± 0.36 (the driver-drop set `fix-qQCD` had +1.23 ± 0.44 at 1e-5), gluon
+0.02 ± 0.22 / +0.33 ± 0.32 / +0.21 ± 0.45 (`fix-gQCD`: +1.92 ± 0.56).
Born weight with geometric min s_ij < 1e-2 / 1e-3 / 1e-4: quark 2.5% / 0.42% /
0.029%, gluon 7.2% / 0.76% / 0.019%.

**Slopes, properly.** The slopes in the table above treat the τ_cut points as
independent although they come from the same events. Fitting a slope per run
(over τ_cut ≤ 5e-4) and taking the mean and its error over runs
(scratchpad `slopes.py`): fixed DISENT, QCD total: cut-new 1e-6 / 1e-8 / 1e-10
−0.01 ± 0.05 / +0.01 ± 0.05 / +0.07 ± 0.06, edge s_min 0.01 / 0.001
+0.03 ± 0.05 / +0.05 ± 0.06, invariant +0.11 ± 0.06, invariant with edge cut
−0.05 ± 0.06; channels: geometric quark +0.09 ± 0.10, gluon +0.14 ± 0.12,
invariant quark +0.03 ± 0.06, gluon +0.22 ± 0.10, abelian +0.11 ± 0.05. Driver
drop: `fix-qcd` +0.25 ± 0.06, `cut-old-1e-6` +0.56 ± 0.06.
**Correction (same night):** I first wrote that all fixed-DISENT sets are
consistent with no slope, but had left out the C_A = 0 quark set (`inv-qA0`):
+0.17 ± 0.04 over 5e-4 … 1e-5 (4σ). That slope is the approach from below
(−0.78 ± 0.05 at 5e-4, −0.36 at 2e-4, −0.18 at 1e-4, 0.00 ± 0.13 at 3e-5,
−0.10 ± 0.18 at 1e-5), i.e. power corrections: over τ_cut ≤ 1e-4 its slope is
+0.04 ± 0.07 (abelian −0.03 ± 0.07, gluon +0.23 ± 0.14, QCD invariant
+0.10 ± 0.10). So a slope over 5e-4 … 1e-5 is not a clean test of a log
mismatch for the invariant sets, which converge from below; the offsets at
3e-5 and 1e-5 are all consistent with zero. The driver-drop sets rise
monotonically and cross zero (4σ and 9σ slopes, offsets +1.1 and +3.1 at 1e-5).

## 2026-10-02 (night) — NNLO DIS 1+1 three ways (AK: slicing+P2B, pure slicing, DISENT+P2B)

Purpose (AK): judge how realistic NNLO 2+1 by slicing and N3LO 1+1 by
slicing+P2B are, from the analogue one order lower.

**Set-up** (`slicing/nnlo11.f90`, commit 3c595f2; x = 0.01, Q² = 400,
√s = 318 GeV, photon exchange, μ_R = μ_F = Q, NNPDF30_nlo_as_0118). One DISENT
run (fixed DISENT plus its two-parton Born given to USER, scratch copy
`libdisent_born.f`) gives for lab-frame jet observables (anti-k_t R = 1,
p_t > 5 GeV, −1 < η < 2.5; Born jet p_t = 15.55 GeV, η = −0.335; bins: total,
≥ 1 jet, leading-jet p_t and η, ≥ 2 jets) the O(αs²) coefficient per Born:
1. DISENT + P2B: DISENT's O(αs²) 2+1 weights minus their Born projection, plus
   C₂^incl O(Born) from hoppet (C₂^incl = −34.2727, C₁^incl = −2.49506 per Born
   in (αs/2π)^k, stable to 4e-5 in the grid; `nnlo11_ref.py`);
2. τ₂ slicing + P2B: the same with the 2+1 part sliced (invariant measure),
   absolute cuts τ₂ > c or relative cuts τ₂ > ρ τ₁;
3. pure τ₁ slicing: τ₁ = min over partitions [2x P·p(beam) + m²(jet)]/Q²
   (recoil-free, thrust-like at leading power), DISENT's O(αs²) 2+1 weights with
   τ₁ > τ_cut plus the NNLO leading-power cumulant times the Born.

**NNLO leading-power cumulant** (`slicing/mcfm_lp/`): MCFM 10.3's 0-jettiness
pieces with the DY structure: beam a = quark beam function (xbeam1bis/xbeam2bis,
z-integrated with adaptive Gauss–Kronrod, e_q² weighted), "beam b" = quark jet
function (MCFM's `jetq` is in αs/4π: its own 1-jettiness assembly uses J1/2,
J2/4 — my first version missed this and had half the double log), soft = qq̄
0-jettiness soft function (equal to the DIS hemisphere soft function at
O(αs²), Kang–Labun–Lee 1501.04110), hard = `hardqq` at Qsq = −Q² (real logs:
the spacelike form factor). Q_B = Q_J = Q, so L = ln τ_cut at μ = Q. Check: the
NLO part equals the τ₁ᵃ cumulant of `tau1b_nlo.py` (validated deterministically
on 1 Oct) to 2e-6 relative at τ_cut = 1e-1 … 1e-6.

**Production**: 240 × 2M events (thservs, nice 10, 24 CPU-min per job, of which
the 2+1 leading-power weights take 90%).

**Findings** (240 runs unless stated; per Born, (αs/2π)²):
- Slicing + P2B against DISENT + P2B (difference with paired errors), largest
  pull over the 15 non-trivial bins: absolute τ₂ > 1e-4: 4.5σ (p_t 16–17),
  1e-5: 2.9σ (1.1σ with CUTOFF 1e-10); relative ρ = 3e-3: 2.2σ (3.3σ in one
  bin of the 1e-10 sample), 1e-3: 1.1σ; with CUTOFF 1e-10 also ρ = 3e-4, 1e-4
  (2.1σ, 1.7σ). (About 165 comparisons, so a few pulls above 2σ are expected.)
  The errors of the difference are 0.1–6 per Born where DISENT + P2B has
  0.01–0.6. Larger cuts show power corrections, largest next
  to the Born p_t (abs 1e-3: −38, −56 in the p_t 14–15, 16–17 bins; ρ = 3e-2:
  −17, −27). Smaller relative cuts (ρ ≤ 3e-4, i.e. absolute cuts down to 1e-7
  and below for small τ₁) drift (Born p_t bin +28 (3.4σ) and +83 (7.7σ)).
- Pure τ₁ slicing, total (exact −34.27): power corrections +25 to +31 at
  τ_cut = 1e-1 … 1e-2 (80–90% of C₂), consistent at 1e-3 (+3 ± 10), then
  +41 ± 24, +43 ± 50, +364 ± 119, +1357 ± 225 at 3e-4 … 1e-5. Slice by slice
  (leading power vs DISENT above the cut) the mismatch switches on below 1e-4
  and grows roughly like 1/τ₁ (≈ 0.014/τ₁ per Born), which no wrong log
  coefficient can produce (it would show at 1e-3, where the result agrees to
  ±10). The variance errors agree with the run-to-run scatter and the medians
  with the means (no heavy tails). NLO, same events: −2.540 ± 0.030 at 1e-3,
  −2.76 ± 0.10 at 1e-4 against −2.4951 (2.4σ); the slice 1e-5…3e-5 agrees
  at NLO (−0.14 ± 0.3). Hypothesis: a DISENT technical cutoff (CUTOFF = 1e-8)
  effect relative to τ₁.
- **Cutoff test (confirmed).** Same set-up with CUTOFF 1e-6 (30 runs) and
  1e-10 (90 runs); pure τ₁ slicing, total, minus exact at τ_cut = 3e-4 / 1e-4 /
  3e-5 / 1e-5: 1e-6: +253 ± 50 / +802 ± 110 / +2562 ± 214 / +5279 ± 401;
  1e-8 (240 runs, final): +28 ± 21 / +29 ± 44 / +401 ± 103 / +1330 ± 196;
  1e-10: +48 ± 37 / −31 ± 81 / +103 ± 247 / −400 ± 429 (and +14.7 ± 13.9 at
  1e-3). Slicing + P2B, Born p_t bin, relative cut ρ = 1e-4: +2063 ± 20,
  +85 ± 10, +11.5 ± 16.7. The deviation grows roughly like CUTOFF^(1/4)
  (factor ≈ 4 at 1e-5 from 1e-8 to 1e-6), the scaling of a sampling
  probability: most likely DISENT's z and x cuts in GENFOR/COLFOR (events
  whose real emission falls in the cut region lose it), which matter once the
  slicing variable gets close to the region they remove. With 1e-10 pure τ₁
  slicing is consistent with the exact NNLO total at every τ_cut ≤ 1e-3, and
  slicing + P2B agrees down to ρ = 1e-4. The NLO trend at 1e-8 was a
  fluctuation (1e-10: −0.011 ± 0.043, −0.005 ± 0.084, +0.29 ± 0.15 at 1e-3,
  3e-4, 1e-4).
- **DISENT + P2B itself is cutoff independent** within the statistics: 1e-10
  vs 1e-8 largest pull 1.9σ in 16 bins (final: χ² = 18.9/15 at O(αs²), largest
  NLO pull 1.6σ); 1e-6 vs 1e-8 3.0σ in the Born p_t bin
  but not monotonic (1e-10 is on the same side as 1e-6). disorder's default
  1e-8 is fine for P2B; the subtraction integrand (O − O_proj) suppresses the
  near-1+1 region where the cutoff acts.
- The τ₁ < 1e-7 skip of the 2+1 leading-power weights has no effect (same
  seeds with 1e-7 and 1e-10: identical sums, the same 156273 evaluations).
- Pure τ₁ slicing in fiducial bins: bins whose τ₁ range lies above the cut are
  DISENT exactly; bins next to the Born p_t have large fiducial power
  corrections (+48, +72 at 1e-3, +13, +25 at 3e-4: the leading power puts all
  events below the cut at the Born kinematics).
- **Cost** (errors per Born from 240 × 2M events): DISENT + P2B 0.006 (≥ 1
  jet), 0.59 (Born p_t bin), 0.31, 0.23, 0.13, 0.10 (p_t 16–17, Born η bin,
  η −0.2…0.2, ≥ 2 jets); slicing + P2B with τ₂ > 1e-4: 0.040, 2.2, 1.2, 1.4,
  0.94, 0.83; with ρ = 1e-3: 0.060, 5.3, 3.8, 2.1, 1.04, 0.74; pure τ₁ at 1e-3:
  8.7 in every bin that contains the Born, DISENT's errors elsewhere. So slicing
  + P2B needs 15–150× the events of DISENT + P2B for the same error, plus the
  leading-power CPU (now about 10× DISENT per event, mostly the soft-function
  angular integrals; tabulable on a two-dimensional grid).

- **Second kinematic point** x = 0.05, Q² = 1000 (`nnlo11 -relbins`, commit
  11d2796: p_t bins scaled and η bins shifted to the Born jet, p_t = 28.3 GeV,
  η = 0.96; jet cuts unchanged; DISENT CUTOFF 1e-10; 180 × 2M). hoppet:
  C₂ = −13.7365, C₁ = −2.23893 (dy = 0.05 was off by 0.2% here; converged at
  dy ≤ 0.0125; x = 0.01 unchanged). Slicing + P2B, largest pull over 15 bins:
  τ₂ > 1e-5: 1.8σ; ρ = 3e-4, 1e-4: 1.6σ, 1.5σ; ρ = 1e-3: 2.8σ (≥ 2 jets, four
  bins above 2σ); ρ = 3e-3: 7.2σ; τ₂ > 1e-4: 6.8σ. So the relative cut is not
  uniformly better: for 2+1 configurations with large τ₁ it is looser in
  absolute terms than τ₂ = 1e-5. Pure τ₁ slicing, total − exact: +12.3, +14.9,
  +10.6 ± 1.3 at 1e-1, 3e-2, 1e-2 (about 100% of C₂), +6.4 ± 4.0 at 3e-3,
  −23.5 ± 10.4 at 1e-3 (2.3σ), consistent below (−1.5 ± 22.4, −6.6 ± 49.3,
  −55 ± 126, −290 ± 266). NLO converges (+0.007 ± 0.007 at 1e-2).
- Page: https://claude.ai/artifact/GD64FDm4dLzUzoYRrHf7Ao (version 4).

**What it says** (AK's question):
- N3LO 1+1 by slicing + P2B: the one-order-lower analogue works. The sliced
  NLO 2+1 inside P2B reproduces the subtraction result in all bins at both
  points with τ₂ > 1e-5 (or ρ ≲ 3e-4 with DISENT's cutoff at 1e-10), at 4–12×
  the statistical error. The N3LO version needs
  the NNLO 2+1 part by τ₂ slicing (two-loop γ*qg hard function, NNLO beam, jet
  and three-direction soft functions, and an NLO 3+1 calculation above the
  cut), but inside P2B it is an N3LO correction and needs much less relative
  precision.
- NNLO 2+1 by pure slicing: hard. Already for 1+1 the NNLO power corrections
  are 75–90% of C₂ at τ_cut ~ 1e-2, the fiducial ones next to the Born are
  larger (+48, +72 per Born at 1e-3), and at 1e-3 the cancellation against the
  leading power costs 15× (Born p_t bin) to 1400× (≥ 1 jet) the subtraction
  error. For 2+1 the Born region is present in every bin, the soft function
  depends on the Born, and the above-cut NLO 3+1 must be numerically clean far
  below τ_cut (DISENT at 1e-8 is not below τ₁ ≈ 3e-5).
- **Results for plotting** (AK asked, 2 Oct morning): every bin, method, cut,
  order and DISENT cutoff of both points in
  `slicing/results/nnlo11/nnlo11_results.csv` (long format: value, error,
  difference to DISENT + P2B with its paired error), figures and
  `plot_nnlo11.py` in the same directory: `nnlo11_<point>_jets` (C₂ per
  leading-jet p_t and η bin, four methods, differences below),
  `nnlo11_<point>_pulls[_cut1e-10]` (slicing + P2B minus DISENT + P2B for every
  bin and cut), `nnlo11_pure_tau1_total` (pure τ₁ slicing, total, against
  τ_cut for the three DISENT cutoffs). Bins: lab frame, anti-k_t R = 1, jets
  with p_t > 5 GeV, −1 < η < 2.5; total, ≥ 1 jet, leading-jet p_t (5, 10, 14,
  15, 16, 17, 20, 40 GeV at x = 0.01), leading-jet η (−1, −0.6, −0.4, −0.2, 0.2,
  1, 2.5), ≥ 2 jets; at x = 0.05 the p_t edges scaled by 28.33/15.553 and the η
  edges shifted by +1.293.

## 2026-10-02 (morning) — tabulated soft function and beam coefficients

AK: "do the tabulation first and get that to work and tell me what speed-up it
leads to" (scale variations and W/Z not needed for now).

**What** (`slicing/mod_slicing_scet.f90`; options in `tau2_nlo` and `nnlo11`,
off by default):
- `-softtable FILE`: `soft_G` interpolates G(a,b) = I0(a,b) ln a + I1(a,b), the
  only expensive piece of the soft function (six calls per three-parton Born),
  from a table on (w, v): v = ln b, w = u + asinh(u/ε(v)), u = ln a,
  ε(v) = e^{v/2}/(1 + e^{v/2}). G has a near-logarithmic singularity at a = 1
  on the scale √b; in w it is smooth, so 0.025 in w and 0.05 in v suffice
  (|u|, |v| ≤ 21, outside direct; 2825 × 840 nodes, 19 MB). Built once at GK
  tolerance 1e-10 (80 s) into FILE (temporary file + rename, safe for
  concurrent jobs) and read by later runs (0.05 s). Catmull–Rom bicubic.
  Grid studies (scratch `gridtest4.f90`): maximum error of G 8e-4 at
  (0.05, 0.1), 2.2e-4 at (0.025, 0.05), at the worst points next to a = 1;
  typical errors much smaller.
- `-beamtable`: the beam coefficients c0, c1, c2 (`beam_coeffs`, 64-point
  convolutions, 11.7 µs per call) tabulated at the run's fixed Q, cubic in
  ln(η/(1−η)), spacing 0.001 (0.25 s); used when the event's Q agrees with
  that Q to 1e-6. Error ≤ 1.7e-6 of the largest coefficient (x = 0.01, Q = 20
  and x = 0.05, Q = 31.6; 1.5e-6 at 0.001, 5e-6 at 0.002, 4e-5 at 0.005).
- `test_scet` part 4 checks both against direct evaluation.

**Checks**: without the options, `nnlo11.dat` and `tau2_nlo.dat` are bitwise
identical to the production binaries (nnlo11 at 11d2796, tau2_nlo_fixed3, both
measures). With both tables, same events (200k, x = 0.01, Q² = 400): nnlo11's
sliced O(αs²) pieces (E2s, D2s) differ from direct evaluation (GK 1e-6) by at
most 4.7e-6 per Born, 4e-7 of the MC error of that run (about 2e-5 of the
error of the 240 × 2M production); everything not involving the
leading power is identical. tau2_nlo outputs: at most 8.6e-6 (geometric),
1.7e-6 (invariant) relative to direct at 1e-9; in units of the MC error of
that run, at most 8.7e-6 (geometric) and 1.3e-6 (invariant) for the below-cut
sums and 1.1e-6, 2.8e-7 for the differences to DISENT; above the cut
identical. Direct at 1e-6 is closer to
1e-9 (1e-8, 4e-8): the table is less accurate than GK at 1e-6, but both are
far below the statistical errors. (A table with 0.0125 in w would cut the
error by about 4 at 76 MB.)

**Speed** (thA371a, one job alone, 200k events, CPU seconds, two passes
agreeing to 1%):

| | direct, GK 1e-6 (production) | direct, 1e-9 (default) | soft table | soft + beam tables | speed-up |
|---|---|---|---|---|---|
| nnlo11 | 71.5 (below cut 68.6) | – | 5.3 (2.4) | 3.7 (0.54) | 19× |
| tau2_nlo, geometric | 33.6 | 60.4 | – | 2.24 | 15× (27× vs 1e-9) |
| tau2_nlo, invariant | 29.6 | 55.5 | – | 2.25 | 13× (25× vs 1e-9) |

The below-cut part of nnlo11 is 127× faster (440 → 3.4 µs per evaluation);
the runs are now dominated by DISENT itself (10–15 µs per event). Memory
37 MB instead of 12 MB.

- The soft table does not depend on Q or the PDFs and works as it is for runs
  integrated over x and Q². The beam table is per Q; for varying Q it needs a
  second dimension, or hoppet convolutions on its grid (the coefficients are
  convolutions of the PDFs with fixed kernels).
- The NLO 2+1 and NNLO 1+1 results above were obtained with direct
  evaluation; nothing to redo.

## 2026-10-02 (afternoon) — start of DIS 3+1 at NLO / 4+1 at LO

AK: "push both and start on LO 4+1 / NLO 3+1" (after the VBF scale
variations were handed over to a cluster session). Branch `2026-10-dis31`
(off `2026-10-tau-slicing`); plan in `docs/dis31-plan.md`.

- Reference for the matrix elements: NNLOJET v1.0.2's public core library
  exports the Z/γ* + partons functions needed, with their colour-ordered
  subleading pieces: B2g0Z, Bt2g0Z, C0g0Z, D0g0Z (3+1 trees), B3g0Z,
  Bt3g0Z, Btt3g0Z, C1g0Z, D1g0Z (4+1 trees), B2g1Z, Bt2g1Z, Btt2g1Z, C0g1Z,
  D0g1Z (3+1 one loop), and a DIS crossing `Bt2g1ZDIS`
  (`~/work/disorder-comparisons/nnlojet-v1.0.2/build/libnnlojet_core.so`).
- Source for our own implementation: MCFM 10.3's Z+2 jet routines
  (`src/Z2jet`), crossed to DIS.

**3+1 trees, photon exchange (2 Oct, evening).** MCFM's Z+2 jet tree
routines, crossed to DIS (`dis31/mcfm`, rules in its README), agree with
DISENT's MATFOR, summed over the labellings of the outgoing partons, up to
one constant: γ* g → q q̄ g to 3e-14 (20 random points), γ* q → q g g +
q Q Q̄ (five flavours, identical-quark interference included) to 1e-13 for
d and u. Piece by piece (q g g, boson on the incoming line, on the pair,
identical-quark interference) each agrees separately with the ERT
functions. Pitfalls: MCFM's `dot` clashes with DISENT's `DOT` (inlined);
the four-quark routine has the leptons fixed in slots 3, 4; the
identical-quark interference needs MCFM's own (phase-consistent)
construction. Harnesses: `dis31/tests`.

**me31: 3+1 trees with absolute normalisation (2 Oct, night).**
`dis31/me31.f90`: |M|^2 for a given flavour assignment of the four partons
(photon exchange), averaged over the incoming lepton and parton, in DISENT's
MATFOR normalisation (α = 1/137, divided by (αs/2π)²), no final-state
symmetry factors. Against MATFOR, summed over the labellings with 1/2 for
identical pairs: ratio 1 to 7e-14 for every incoming flavour −5..5 (50
random points; `dis31/tests/harness_me31.f90`, stops with an error above
1e-10). The normalisation follows from MCFM's colour factors and the
averages: |M|²/MATFOR = (αs/2π)² with α = 1/137.
Timing (thA371a, one point): MATFOR 2.6 µs (all flavours, fixed labels);
me31 41.6 µs for all 71 flavour assignments of one labelling, dominated by
recomputing the spinor products 71 times and by amplitudes that differ
only by charges. To do in the Born interface: spinors once per point,
q g g, g → q q̄ g and the four-quark A, B once per labelling, charges
afterwards.

**Correction: me31's four-quark channels were wrong pointwise (2 Oct,
late evening).** The two entries above claim agreement with MATFOR for the
four-quark channels. That is true for the sum over the labellings of the
outgoing partons, but that sum is blind to the charge-odd e_q e_Q term
(photon on the incoming line times photon on the pair). The term is odd
under Q ↔ Q̄ and cancels in the sum. me31 had this term with the wrong
sign: `ampqqb_qqb(2,1,5,6)` reads the incoming line opposite to MCFM's
orientation (MCFM's `qqb_z2jet` uses (1,2) for that channel), and
reversing a quark line flips the relative sign of the two photon couplings.
For identical quarks, me31 had also dropped the A·B cross terms (of the
same odd type). The "piece by piece" agreement claimed above was for the
symmetrised sums and therefore did not test these terms.
Found because the 4+1 trees (below), built with MCFM's own conventions,
did not factorise onto me31 in the q ∥ g limits of the four-quark
channels. The squares with only A or only B did factorise, which isolated
the problem to the interference.
Settled independently: `dis31/tests/fd31.py` evaluates γ* q → q Q Q̄ from
Feynman diagrams (explicit Dirac matrices and helicity spinors, physical
momenta, explicit SU(3) colour sums; two diagrams with the photon on each
line, and for identical quarks minus the same with the two quarks
exchanged). Fixed me31:
- non-identical: e_q A − e_Q B;
- identical: D = A − B and E = A_e − B_e from `ampqqb_qqb(5,1,2,6)`. The
  interference is (2/N) Re D(j,swap(j),j₃) E*(j,swap(j),j₃) over equal
  quark helicities, with the helicity labels found by a least-squares fit
  of the 32 possible D·E* terms to the Feynman diagrams. Exactly four terms
  have coefficient 1 and all others 0.
Now me31/FD = 1 to 1.3e-14 pointwise (d → d u ū, u → u u ū, d̄ → d̄ ū u;
10 points). The MATFOR sum test still gives 1 (2.6e-14). Consequence for
physics: none for flavour-blind observables (the term integrates to zero
under Q ↔ Q̄), but the subtraction needs the pointwise matrix element.

**4+1 trees, photon exchange (2 Oct, late evening).** `dis31/me41.f90`:
|M|² for a flavour assignment of the five partons (DIS layout P(4,8): 1
incoming parton, 2–5 outgoing, 6 q, 7 and 8 the leptons), normalisation
as me31 with one more power of αs/2π. From MCFM 10.3:
- `xzqqggg` (q q̄ g g g) for q → q g g g and g → q q̄ g g;
- `msq_ZqqQQg` (q q̄ Q Q̄ g) in a photon version `msq_gqqQQg` (with
  `makemb_photon`, `nagyqqQQg`), with the line charges as arguments, for
  q → q Q Q̄ g, identical quarks, and g → q q̄ Q Q̄.
The crossing follows MCFM's `qqb_z2jet_g` calls with MCFM's own sign
conventions; rules in `dis31/mcfm/README.md`.
Check (`dis31/tests/harness_lim41.f90`): me41 against me31 in
single-collinear limits, for every channel (20 limits, 3 random points
each):
- final state: q ∥ g, g ∥ g, g → q q̄ (CS FF map, y → 0, averaged over
  φ and φ + π/2 for the gluon splittings);
- initial state: q → q g, g → q q̄ (IF map, u → 0; averaged matrix
  elements with the averaged AP kernels).
All ratios tend to 1 linearly in y; at y = 1e-10 the largest deviation is
1e-4.
Timing (thA371a): 24 µs per call for q g g g and q q̄ g g, 33 µs for the
four-quark channels (MATFOR: 2.6 µs). To optimise in the integrand:
- spinors once per point;
- the e_q², e_Q², e_q e_Q pieces once per momentum assignment, not per
  flavour;
- the non-identical case without the exchange amplitudes (8 → 2 calls to
  `nagyqqQQg`).
If MCFM's trees remain the bottleneck, consider our own helicity-summed
trees with FORM (AK).
Next: NNLOJET pointwise (B3g0Z/C1g0Z/D1g0Z) for the 4+1 trees, then
colour- and spin-correlated 3+1 Borns and the CS dipoles.

**3+1 and 4+1 trees against NNLOJET, pointwise (2 Oct, night).**
Harnesses in `~/cernbox/disorder-comparisons/dis31_nnlojet` (outside the
repository; README there). NNLOJET v1.0.2's public Z+partons functions,
photon only (`igamma_proc = 1`), at the same random DIS momenta. Charge
structure resolved into the e_q², e_Q², e_q e_Q coefficients (three charge
assignments at the same momenta).
- me31 four-quark (d → d u ū) against C0g0Z, as NNLOJET's DIS real calls it:
  all three coefficients 1.000000000000. This is a third independent
  confirmation of the corrected sign of the charge-odd term.
- me41 against NNLOJET, all to 7e-14 including the absolute normalisation:
  - q → q g g g and g → q q̄ g g (B3g0Z, Bt3g0Z, Btt3g0Z);
  - d → d u ū g and g → d d̄ u ū (C1g0Z, Ct1g0Z, Ctt1g0Z; all three
    coefficients);
  - identical quarks d → d d d̄ g and g → d d̄ d d̄ (+ D1g0Z, Dt1g0Z).
- One wrong turn on the way: I first took the colour weights from
  NNLOJET's `FullC1g0Z` (+(Ct − Ctt)/N²) and found disagreements of tens of
  percent. NNLOJET's DIS process uses −(Ct − Ctt)/N² (`qcdnormDIS.f`), and
  with that everything agrees. Other pitfalls are listed in the README:
  the `astore` amplitude cache, quark types in /CZFlav/, and NNLOJET's
  pair-slot order.
- NNLOJET's DIS symmetrises its incoming-quark four-quark channels over
  Q ↔ Q̄ (sC1g0Z), which drops the charge-odd term. That is fine for
  flavour-blind observables.
Status: the photon-exchange trees for DIS 3+1 and 4+1 (me31, me41) are
validated pointwise in every channel by three independent references:
- Feynman diagrams (3+1 four-quark);
- DISENT's MATFOR (3+1, symmetrised);
- NNLOJET (3+1 four-quark, all of 4+1);
plus all single-collinear limits 4+1 → 3+1.

**Colour- and spin-correlated 3+1 Borns (2 Oct, night).** `dis31/born31.f90`:
- `born31_cc` returns |M|² (equal to me31) and ⟨T_i·T_k⟩ for all parton
  pairs (all partons outgoing, CS conventions).
- `born31_sc` contracts a Born gluon with a vector n (colour-summed and
  colour-correlated), for the gluon-splitting dipoles.

Construction:
- The colour matrices ⟨c_m|T_i·T_k|c_n⟩ for the q q̄ g g basis (T^A T^B,
  T^B T^A) and the four-quark basis (direct, exchange) are computed exactly
  with explicit SU(3) (`dis31/tests/colour.py`); colour conservation is
  checked there to 1e-15.
- They are combined with MCFM's colour-ordered amplitudes: `subqcd`, and
  `subqcdn` with one gluon contracted with n (ported from MCFM `src/W2jet`
  with `spinork` and `checkndotp`).

Assignments fixed by the tests:
- `subqcd`'s ordering (A, B) belongs to (T^A T^B)_{q q̄}.
- `subqcdn` has the opposite assignment, and its contraction is normalised
  to half the polarisation sum. The polarisation-sum identity showed both:
  |M|² alone cannot see the ordering, because the metric is symmetric.
- An incoming antiquark keeps the quark's roles (C swaps q ↔ q̄ and the
  colour order).
- The gluon channel needs MCFM slots (2,5,1,6). My first version
  hard-coded the quark channel's (2,1,5,6); the soft test caught it.

Tests (`dis31/tests/harness_born31.f90`):
- msq = me31 to 4e-16;
- colour conservation to 1e-15;
- polarisation sums to 3e-14;
- soft-gluon limit of me41 against −8π² Σ p_a·p_b/(p_a·q p_b·q) ⟨T_a·T_b⟩,
  nine channels (quark, antiquark, gluon-initiated, four quarks, identical),
  to 2e-4 at λ = 1e-5, converging linearly;
- collinear limits at fixed azimuth (three φ, no averaging), which need the
  spin correlations, to 5e-4 at y = 1e-9: final g → g g and g → q q̄;
  initial q → g + q and g → g + g.

The initial-state q → g kernel is CF[−g x − 4(1−x)/x k⊥k⊥/k⊥²]: that sign
of the k⊥ term is the one that matches the azimuthal dependence and, when
averaged, gives CF(1 + (1−x)²)/x. To keep in mind for the IF dipoles.

**CS dipoles for DIS 4+1 → 3+1 (2 Oct, night).** `dis31/dip41.f90`:
`dip41_list` returns every dipole with its mapped Born momenta, flavours
and value, ready for integration with a jet function on the mapped Born.
- Types: FF D_{ij,k}, FI D_ij^a (the incoming parton as spectator), IF
  D^{ai}_k (incoming q → q g, q → g q, g → q q̄, g → g g). No II dipoles
  (one coloured incoming parton).
- Formulas from CS section 5 at ε = 0. Colour and spin correlations come
  from born31, with the spin vectors (z̃_i p_i − z̃_j p_j) and
  (p_i/u − p_k/(1−u)). They are orthogonal to the Born gluon (checked
  algebraically), as born31_sc requires.
- Normalisation as me41: 8π α_s → 16π².

Test (`dis31/tests/harness_dip41.f90`): me41 / Σ dipoles → 1, 3 random
points, worst 1.9e-4 at the smallest parameter:
- 11 final-state collinear limits at fixed azimuth (q ∥ g, g ∥ g,
  Q ∥ Q̄, identical quarks, gluon-initiated);
- 7 initial-state limits (q → q, q → g, q̄ → g, g → q, g → g);
- 7 soft-gluon limits.
Dipoles whose mapped Born is itself unresolved (an invariant below
1e-3 W²) are left out, as in DISENT's test_subtraction. Their Born still
holds the limit pair, and the jet function removes them in a calculation.
My first run without that filter gave ratios around 0.5; that was the test
set-up, not the dipoles.

**3+1 one loop (3 Oct, early morning).** `dis31/virt31.f90` with the BDK
one-loop amplitudes ported from MCFM 10.3 (`dis31/mcfm/loop`).
- The 51 source files are the dependency closure of the routines used,
  computed from MCFM's object files.
- Left out: the exact top loops, and the boson-on-loop (vector and axial)
  pieces. The latter vanish for a photon by Furry's theorem: γ*gg through a
  quark loop has the symmetric colour factor Tr(T^a T^b).
- MCFM's pole convention is epinv = epinv2 = 1/ε, with the double pole
  given by epinv·epinv2. The Laurent coefficients come from evaluating at
  1/ε = 0, 1, −1. My first extraction set epinv2 = 1/ε² and gave a zero
  double pole.

`virt31_ren` returns 2Re⟨M0|M1⟩ renormalised in MS-bar, in the HV scheme,
with (4π)^ε/Γ(1−ε) factored out (CS's normalisation), for n_f = 5.

Checks:
- tree = me31 in every channel;
- the poles are those of −⟨I⟩: the double pole −Σ_i C_i|M0|², and the
  single pole Σ⟨T_i·T_k⟩ ln(μ²/2p_i·p_k) − Σ_i γ_i|M0|² with born31's
  colour correlations; nine channels, two scales, ≤ 7e-13
  (`tests/harness_virt31`).
- Finite part against NNLOJET's one-loop functions, using its DIS
  real-virtual colour weights (qcdnormDIS): all channels (q → q g g,
  g → q q̄ g, four quarks with each charge structure, identical quarks) at
  μ² = 130 and 500. Ours − NNLOJET is one constant × tree per channel type:
  −(π²/12)ΣC_i + C_F (q q̄ g g) and −(π²/12)ΣC_i + 2C_F (four quarks).
  - The π² term is NNLOJET's normalisation (e^{−γε}-type instead of
    1/Γ(1−ε)); its code comment calls this "a correction for C(ε)".
  - The rest is the scheme conversion of MCFM's raw output. FDH → HV is
    −Σγ̃_i with γ̃_q = C_F/2, γ̃_g = C_A/6. For q q̄ g g the DRED → MS-bar
    coupling conversion (+N/3) is still to be applied (MCFM does it in
    qqb_z2jet_v's subuv); a6routine has already applied it for four quarks.
  - Net conversions in virt31_ren: q q̄ g g −C_F × tree (and the UV pole
    −2β₀/ε), four quarks −2C_F × tree. Two independent channel types fix
    γ̃_q and γ̃_g consistently.
- Harnesses for the NNLOJET comparison:
  `~/cernbox/disorder-comparisons/dis31_nnlojet/harness_v31{,q}.f90`.

This completes the matrix-element ingredients of NLO 3+1: trees (me31,
me41), correlated Borns (born31), dipoles (dip41) and one loop (virt31).
Still to do: I, K, P (the integrated dipoles, from CS with born31's
correlations), the 3+1 phase space and the integration.

**Integrated dipoles I, K, P (3 Oct, morning).** `dis31/iop31.f90`, from
CS section 8 (one incoming hadron; eqs. 8.25, 8.38, 8.39 with 7.27, 7.28,
8.32–8.35; MS-bar, K_F.S. = 0, n_f = 5). The paper is in the scratchpad
(hep-ph/9605323).
- ⟨I(ε)⟩ as Laurent coefficients, in the same normalisation as virt31_ren
  ((4π)^ε/Γ(1−ε) factored out).
- K + P per Born point reduces to three numbers:
  B = |M_b|², G = Σ_i γ_i/T_i² ⟨T_i·T_b⟩ and
  L_P = Σ_i ⟨T_i·T_b⟩/T_b² ln(μ_F²/2p_b·p_i).
  Then K^{a,b}(x) + P^{a,b}(x) = K̄^{ab}(x) B
  + δ^{ab}[(1/(1−x))₊ + δ(1−x)] G + P^{ab}(x) L_P, with K̄ and P split into
  regular, plus-distribution and δ parts (`iop31_kernel`).
- DISENT's KPFUNS has the same structure for 2+1, but there the
  correlations are proportional to the Born and its δ term contains part
  of DISENT's own I convention, so it is not reused.

Tests (`tests/harness_iop31`, nine channels, two scales):
- the poles of V + I cancel to 1e-13 per tree;
- ∂L_P/∂ln μ_F² = −B to 1e-14 (CS 8.42);
- for four quarks G/B = −2 exactly.
Not yet tested: K + P against an integrated real minus dipoles, which needs
the phase-space integration (next).

**NLO 3+1 integrated: first results and LO validation (3 Oct, night).**
AK: "Keep working over night and see if you can finish the 3+1 NLO … if
you hit [the cluster] point you can stop and write instructions."

`dis31/nlo31.f90`: an integrator for σ(e p → e + ≥ 3 jets) at LO and NLO
(photon exchange).
- Breit-frame phase space: Q² (log), y, η (log), the hadronic system by
  sequential two-body decays.
- Inclusive kt jets (R = 1, E-scheme) in the Breit frame; our own VEGAS.
- Parts: lo, vi = V + I, kp = K + P (x-convolution with the
  plus-distribution integrals below ξ done analytically), r = R − Σ
  dipoles, with each dipole's jet function on its mapped Born.
- The flavour sums use that all pieces are bilinear in the quark charges
  and C-symmetric. Per point there are 1–3 evaluations per channel
  topology instead of about 65 (Born) and 80 (real) assignments; checked
  against the explicit sums to 1e-12 (`nlo31 chk`).
- `virt31_finite_only` evaluates only the finite part in production:
  identical result, 2.8 times faster.
- Standalone build: `dis31/build_nlo31.sh`.

Set-up (as NNLOJET's epLJJ runcard `dis31/validation/nnlojet_epLJJ_3j.run`):
27.5 × 920 GeV, photon, α = 1/137, 150 < Q² < 15000, 0.1 < y < 0.9,
p_T,jet > 5 GeV, ≥ 3 jets, μ_R = μ_F = Q, NNPDF30_nlo_as_0118.

| | nlo31 [pb] | NNLOJET [pb] |
|---|---|---|
| LO | 84.852 ± 0.092 | 84.739 ± 0.099 (R channel, 10 M points) |
| V + I | 28.881 ± 0.072 | |
| K + P | 31.092 ± 0.035 | |
| R − dipoles | −30.025 ± 0.455 (4 seeds × 6 M) | |
| NLO correction | 29.95 ± 0.46 (K = 1.353) | RV + RR: needs the cluster |

LO agrees: 0.8σ in total, and every Q² bin within about 0.5%. This tests
the normalisation, phase space, jets and cuts. Pointwise checks added in
`tests/harness_iop31`: the renormalisation-group structure of V + I
(2β₀ ln μ² |M0|² at fixed α_s, 7e-14).

NNLOJET pitfalls found (also in the runcard header):
- `dis_frame = BREIT` gives zero (R) or NaN (RV) in v1.0.2. It also calls
  setfixT(.false.), and the DIS process clusters in the Breit frame by
  default anyway.
- `beam1 = EM beam2 = P …` gives zero for epLJJ; `collider = ep` works.

NNLOJET's RV converges (≈ 56 pb; its antenna terms distribute differently
from ours, so only RV + RR compares). Its RR is very slow (≈ 0.1 s per point
on a loaded core), so the reference needs the cluster:
`dis31/validation/CLUSTER-INSTRUCTIONS.md`.
Technical cut on R − dipoles (smallest 2p_i·p_j / W², default 1e-9;
4 seeds × 6 M points each):
- 1e-7: −29.35 ± 0.50 pb;
- 1e-9: −30.02 ± 0.46 pb;
- 1e-11: −29.06 ± 0.48 pb.
No trend; the largest difference is 1.4σ. All 12 seeds together give
R − D ≈ −29.5 ± 0.3 pb and an NLO correction ≈ 30.5 ± 0.3 pb.

**NLO 3+1 validated against NNLOJET (3 Oct, day).** Runs on the MPP machines
(AK's rules: thservs nice 10, desktops nice 19 with at most half the cores;
dispatcher `dis31/validation/dispatch.sh`).

NNLOJET v1.0.2, runcard `dis31/validation/nnlojet_epLJJ_3j.run`:
- RV: 15 jobs × 300k points from the pilot grid, 56.03 ± 0.14 pb;
- RR: 220 independent jobs (warm-up 3 × 10k, production 40k each),
  −26.04 ± 0.66 pb. These are equal weights with error = scatter/√N.
  - The per-job values have heavy tails on both sides (−65 to +67 pb;
    median per-job error 5.0 pb, scatter 9.8 pb).
  - Inverse variance gives −25.12 ± 0.32 pb, which is biased.
  - Symmetric trimming of 1–5% gives −26.3 ± 0.5 pb.

nlo31: 36 r seeds (12 from the night run plus 24 new, 1 M × 6 each), 6 vi,
2 kp; the vi and kp from the new runs.

| | nlo31 [pb] | NNLOJET [pb] |
|---|---|---|
| NLO correction (≥ 3 jets) | 30.34 ± 0.29 | RV + RR 29.99 ± 0.68 (0.5σ) |

| Q² [GeV²] | nlo31 | NNLOJET | pull |
|---|---|---|---|
| 150–200 | 7.417 ± 0.131 | 7.178 ± 0.288 | +0.8 |
| 200–300 | 7.786 ± 0.122 | 7.606 ± 0.259 | +0.6 |
| 300–500 | 6.624 ± 0.091 | 7.143 ± 0.477 | −1.1 |
| 500–1000 | 4.944 ± 0.082 | 4.500 ± 0.170 | +2.3 |
| 1000–3000 | 2.858 ± 0.055 | 2.766 ± 0.133 | +0.6 |
| 3000–15000 | 0.715 ± 0.023 | 0.794 ± 0.069 | −1.1 |

χ² = 9.2 for 6 bins (p ≈ 0.16). The NNLOJET bin errors are themselves
uncertain because of the tails. With the LO agreement (0.8σ), this
validates the NLO 3+1 (photon exchange) at the 2% level of the correction,
i.e. about 0.6% of the NLO cross section.

**Correction (combination of the r seeds).** The night entry's "R − D ≈
−29.5 ± 0.3 pb, NLO correction ≈ 30.5 ± 0.3 pb" (and the −30.025 ± 0.455
of 4 seeds) used inverse-variance weights of the VEGAS errors.
- For R − dipoles these errors are correlated with the values: low
  fluctuations come with small errors. One new seed is −37.5 ± 1.9 pb, a 4σ
  pull.
- Inverse-variance weighting is therefore biased: for the 24 new seeds it
  gives −29.37 against −29.75 ± 0.40 pb with equal weights.
- `dis31/combine_nlo31.py` now uses equal weights with scatter errors and
  prints inverse variance as a check.
- Corrected values for the night runs (equal weights):
  - technical cut 1e-7: −29.31 ± 0.65 pb;
  - 1e-9: −29.81 ± 0.39 pb;
  - 1e-11: −28.76 ± 0.69 pb;
  - all 12 seeds: −29.29 ± 0.34 pb.
  The technical-cut independence still holds.

Operational note: the first nlo31 batch on the thservs was lost. The
launch line ended in `> /dev/null 2>&1` after the job's own redirect, so the
last redirect won and the output was discarded. The batch was rerun.

## 2026-10-03 (evening) — NNLO DIS 2+1: measure, hard function (branch 2026-10-nnlo21)

AK: "If it all looks right I suggest you move on to the 2+1 NNLO (after
pushing and publishing)." Plan: `docs/nnlo21-plan.md`. Branch
`2026-10-nnlo21` off `2026-10-dis31`, with the DISENT fix 5262826
cherry-picked (slicing must not use the unfixed DISENT). The notebook
conflict was resolved by keeping both entries.

**Measure: geometric in the jets' rest frame** (`-measure cm`, measure = 2
in `slicing/mod_tau2_run.f90`). For each partition of the outgoing partons
into beam / jet 1 / jet 2, let u be the four-velocity of P_J1 + P_J2. Then
T_π = Σ_beam P·p_k/P·u + Σ_jets (u·P_J − |P_J|_u) (covariant), and T₂ = min_π T_π.
- In every singular limit the frame is the Born's partonic CM frame, where
  the two jets are back to back. The NNLO soft function then depends on one
  angle and is known numerically: Bell, Dehnadi, Mohrmann, Rahn
  2312.11626 (grids), and the CEMW fit in MCFM.
- **Change of route against the 2 Oct entry** (invariant measure plus our
  own 2604.13167 soft function). With the invariant measure the soft
  function depends on two variables and is not tabulated anywhere; the
  new measure needs only published results.
- Unit test: boost invariance to 1e-15; the soft limit equals
  min_i n_i·k in the Born CM frame. Jet masses are computed as
  m²/(u·P + |P|), avoiding a cancellation at small mass.
- NLO 2+1 against DISENT (x = 0.01, Q² = 400 GeV², 60 × 2M events,
  `slicing-runs/cm-x0.01-Q400`), all τ bins, (below + above − DISENT)/DISENT:

  | τ_cut | (sum − DISENT)/DISENT |
  |---|---|
  | 5e-4 | −0.080 ± 0.006 |
  | 2e-4 | −0.029 ± 0.008 |
  | 1e-4 | −0.029 ± 0.011 |
  | 3e-5 | −0.024 ± 0.016 |
  | 1e-5 | −0.012 ± 0.022 |

  Each τ bin is within 1σ of zero at 1e-5, and power corrections approach
  from below, as for the Breit-frame measure.

**Two-loop hard function** (`nnlo21/hard21.f90`, README there).
- Two-loop helicity coefficients for (2+1)-jet DIS: Gehrmann, Glover,
  0904.2665. Their arXiv Fortran covers eight DIS regions, but the helicity
  sum of the quark channel also needs the q ↔ q̄ partner of each region
  (as MCFM's `iperm` loop), which is in none of the eight.
- NNLOJET v1.0.2 (GPL) has all 16 regions of the analytic continuation
  with their partners (`B1gNZ.f`, `helcoeff`).
  - Its region coefficients equal 0904.2665's in all eight shared regions:
    all two-loop coefficients exactly, and the one-loop ones after
    a = a_NJ + (11/24)(L13+L23), c = c_NJ − (1/3)(L13+L23).
  - The shift was fitted at 12 points and then checked at 4–16 points per
    region, 30 coefficients, ≤ 1e-6.
- Catani → SCET: C1 = Ω1 + I1 Ω0 and C2 = Ω2 + I1 Ω1 + (I1²/2 + R⁰) Ω0 at
  μ² = Q², complex logs ln(−s_ij/μ² − i0).
  - Derived in `nnlo21/scheme_conversion.py`.
  - Symbolically equal to MCFM's timelike `schemeconvC0` and
    `schemeconv2lM0`. MCFM's factors 2 and 4 belong to its α_s/4π
    coefficients.
- Tests (`nnlo21/tests/harness_hard21.f90`):
  - tree ∝ MATTHR (ratio × Q⁴ constant);
  - **one loop equal to the DISENT-based hard function of the NLO slicing
    to 3e-11 (quark) and 1e-9 (gluon) at 40 random points**;
  - two loop smooth along angular scans, except a 0.2% step at
    2p₁·p₂ = Q² (NNLOJET's displacement of v by 1e-3 at v → 1); negligible
    after integration.

Surveyed for the rest (plan):
- beam functions (MCFM `xbeam*`, `I2qq`, `I2gg`) and jet functions
  (`SCET1j/jet.f90`) at NNLO;
- soft function: MCFM `soft1.f90`, analytic in general y_ij except the
  non-abelian constant (CEMW fit, valid for back-to-back 1, 2). The quark
  channel maps onto MCFM's "qgq" (jets q, g back to back, beam q) and the
  gluon channel onto "qag". Cross-check against 2312.11626's grids.
- assembly: MCFM's `assemblejet`, with our second jet function in place of
  its second beam.

**Note (3 Oct): the Gehrmann–Glover DIS two-loop files miss half of the
quark-channel helicity sum.** Recorded so that it is not forgotten.

Gehrmann, Glover (arXiv:0904.2665) give the hadronic current for one
helicity configuration (q+, g+, q̄−) as coefficients α, β, γ(x, y, z).
- Parity flips all helicities with the same coefficients.
- The other gluon helicity follows from charge conjugation plus parity: the
  same current with the quark and antiquark momenta swapped (x fixed,
  y ↔ z).
- So |M|² summed over helicities needs the coefficients at P and at
  P(p₁ ↔ p₂), as MCFM's `Zampqqbgsq` does with its `iperm` loop.

In e⁺e⁻ → 3 jets (and Z+jet) both points lie in the same analytic region.
In DIS they do not:
- **Gluon channel:** the gluon is incoming on leg 3; the swap stays in the
  region, and 0904.2665's regions 5–8 cover it.
- **Quark channel:** the incoming parton moves from leg 1 to leg 2, i.e. from
  (y > 0, z < 0) to (y < 0, z > 0). The paper's eight Fortran regions
  (its "1d, 2c, 3b, 4d" for lepton–quark) contain only one of each pair. The
  partners 1c, 2b, 4b, 3d are missing.
- The text's labelling of the lepton–quark process (incoming quark = −p₂)
  contradicts the signs of the invariants it states (s₂₃ > 0, s₁₂, s₁₃ < 0
  imply the incoming parton on leg 1). This hides the issue.
- With only those files, half of the quark-channel helicity sum would be
  evaluated outside the range of the 2dHPL representation. `tdhpl` as
  distributed in MCFM 10.3 only prints a warning there (its `stop` is
  commented out).

Fix used in `nnlo21/hard21.f90`: NNLOJET v1.0.2's `helcoeff`/`makecoef`.
- It has all 16 regions and evaluates each together with its partner
  (`kinregion(1:2)`).
- Its coefficients equal 0904.2665's in all eight shared regions: all
  two-loop coefficients exactly, the one-loop ones after
  a = a_NJ + (11/24)(L13+L23), c = c_NJ − (1/3)(L13+L23).
- The one-loop hard function then equals the DISENT-based one to 1e-9 at
  random points. That check uses both members of every pair, including the
  four regions 0904.2665 lacks.

## 2026-10-04 (morning) — NNLO 2+1 pilot: the jets'-frame measure was not IR safe from four partons on

**Correction of the 3 Oct evening entry (measure 2).** "T_π in the rest frame
of the jets of π, T₂ = min_π T_π" is **not IR safe for four or more partons**:
- A partition that takes two collinear partons as the two jets has a nearly
  massless jet pair. Its frame is infinitely boosted, and its T_π → 0.
- With three partons that happens only in a genuine singular region, so the
  NLO tests (DISENT, 3 Oct) were not affected. With four partons, real
  events near a 3+1 collinear limit get τ₂ → 0 while the dipoles' mapped
  Borns keep their finite τ₂.
- Found in the first NNLO τ_cut pilot (54 jobs, x = 0.01, Q² = 400 GeV²,
  `nnlo21/runs/tcut1`):
  - in the all-bins cell at τ_cut = 2e-2, r = −1030 ± 420 pb/GeV² against
    lo = 21;
  - the NNLO sum drifted with τ_cut.
- Diagnostic (`diag_r.f90`): near singular limits the real's acceptance
  differed from that of its largest dipole in 33 of 570 events (τ₂(4) ≈ 2e-3
  against τ₂(3) ≈ 3e-2).

New definition (`tau2_jetframe`, `nlo31` `tau2cm`):
1. the partition that minimises the geometric T₂ in the Breit frame
   (measure 0, IR safe) fixes the frame u = (P_J1 + P_J2)/m;
2. then the exact minimum over partitions, all evaluated in that one frame.

It is the ordinary geometric 2-jettiness in a frame that tends to the Born's
partonic CM frame in every singular limit, so the soft function is unchanged.
- Check: no acceptance mismatches in 141 events with s_min < 1e-8 W². Eleven
  remain at 1e-5 W², at near ties of two Breit partitions (a measure-zero
  boundary, as for any jet algorithm).

Reruns:
- the NLO validation against DISENT with the new measure
  (`slicing-runs/cm2-x0.01-Q400`, 60 × 2M);
- the above-cut NNLO parts (`nnlo21/runs/tcut2`).

The below-cut runs (b1, b2) are unchanged, since the Born frame does not
change.

Also: the pilot's waiter counted 56 jobs, but the directory listing it came
from included two binaries; there were 54 jobs. So the results sat uncombined
overnight. The dispatcher's queue pattern matched the `sliced21` binary in
the run directory and looped with phantom launches; it was stopped.
Starting two dispatchers at once overfilled thserv05/06 (35 and 57 jobs);
the excess was killed by PID and re-queued.

**NNLO τ_cut test with the fixed measure (4 Oct, day).** Set-up: x = 0.01,
Q² = 400 GeV², τ_zQ bins; `nnlo21/runs/tcut2` (above the cut) with tcut1's b1,
b2 (below the cut). Combination: `nnlo21/combine_tcut.py`. Page:
https://claude.ai/artifact/ChKctwxDdEau7H9DGMxwyW

- NLO with our own codes (b1 + lo) against DISENT, all bins: 16.4 ± 0.5,
  17.4 ± 1.3, 18.3 ± 1.2 at τ_cut = 5e-4, 2e-4, 1e-4, against 16.87. Fine.
- NNLO (b2 + vi + kp + r), all bins: near zero and flat from 2e-2 to ≈ 2e-4
  (−5 ± 1, −3 ± 3, 4 ± 5, 5 ± 8, −11 ± 11, −19 ± 15, 9 ± 15). Then it rises:
  81 ± 35, 342 ± 59, 1100 ± 91 at 1e-4, 3e-5, 1e-5. Bin [0.4, 0.5): 0.8,
  0.5, 0.8, 0.7, −3.6, 5, 9, 20, 46, 107 (LO 3.4).
- Technical cut of the real part (drops events with s_ij < c W²). NNLO at
  τ_cut = 1e-5, all bins:
  - c = 1e-7: 1073 ± 140 (15 seeds);
  - c = 1e-9: 1100 ± 91;
  - c = 1e-12: 840 ± 69.

  r itself moves by −260 ± 100 from 1e-9 to 1e-12 at τ_cut = 1e-5, so the
  technical cut has some effect, but it is not the cause (1e-7 = 1e-9).
- Not a wrong log coefficient either: a wrong L² term that gives +1100 at
  1e-5 would give several hundred between 2e-2 and 1e-3, where the sum is
  flat. The mismatch switches on below 1e-4 and grows faster than any log,
  roughly like τ^(−0.7).
- The per-seed distributions of r have no heavy tails (mean = median).
- Hypothesis: sampling. The above-cut 4+1 events with τ₂ just above τ_cut
  need two small invariants at once. `nlo31`'s sequential decays sample the
  mass ratios and cosines uniformly, so these configurations are hardly ever
  generated. Every seed then misses the same region, which gives a bias
  without tails.
- Test: `logmap`, a logistic map of all mass ratios and decay cosines with
  both ends logarithmic down to 1e-12 (mode 1 only; mode 0 unchanged,
  checked bit-identical).
  - On lo it agrees with uniform sampling to 0.5% at large τ_cut, which
    checks its Jacobian. At small τ_cut it has more variance.
  - r with logmap: 30 seeds running (`runs/tcut4`).

**Sampling and dipole-stability checks (4 Oct, afternoon).**
- **Log-mapped `r`** (30 seeds, `runs/tcut4`): at τ_cut = 1e-5, all bins,
  −2321 ± 118 (median −2463), against uniform −2734 ± 57 (technical cut 1e-12).
  - Undersampling would need ≈ −3800 to remove the drift. The log map moves
    the result the other way, by 3σ, and is skewed.
- **Real ME stability under a random Lorentz transformation:** median 1e-9,
  99% quantile ≤ 1e-2 even at s_min/W² ~ 1e-11. Stable.
- **Dipole sum: correction of a statement made to AK the same afternoon**
  ("the problem is our dipoles").
  - In the corners the log map populates (several small invariants at once,
    double unresolved), individual dipoles are 1e4–1e7 times the real and
    cancel among themselves. Their relative precision (1e-5 … 1e-9) is fine,
    but the absolute rounding exceeds R. In 1e-9 < s_min/W² < 1e-6,
    |Δdipole| > 0.1 R occurs in ≈ 86k of 445k log-mapped points, against 1
    case with uniform sampling (total number of uniform points not
    recorded).
  - No dipole bug. `harness_dip41` (fixed azimuth, resolved Borns) stands.
    The pointwise ratio tests that seemed to fail summed dipoles whose mapped
    Born was itself unresolved.
  - Consequence: the log map is unsuitable (it samples rounding-dominated
    corners). The uniform runs, which show the drift, hardly visit them.
- So the small-τ_cut drift is still unexplained. Next: quark/gluon channel
  split (`runs/chan1`, PDF masks `pdfmask` in nlo31 and sliced21/lp21), and
  an RG check of the below-cut cumulant.

**Cause of the small-τ_cut drift: VEGAS adaptation bias above the cut (4 Oct, night).**
- **Quark/gluon split** (`runs/chan1`, `chan2`; PDF masks):
  - The drift is present in both channels: NNLO at 1e-5 is 406 ± 36 (q) and
    582 ± 39 (g), after a roughly flat plateau.
  - Additivity q + g = all holds within 0.5σ for b1, b2, kp at every τ_cut,
    and within 1.5σ for vi. It fails for lo: −3.9 (2.0σ), −9.9 (2.3σ),
    −19.4 (3.3σ) at 1e-4, 3e-5, 1e-5. r: −50, −88, −157 (1.2–1.7σ, same
    sign).
  - Fixed in the process: `sliced21` with a mask had dropped off-diagonal
    beam terms (Born × c/f0 with f0 = 0). It now multiplies the unit
    matrix element by c directly. Unmasked results are unchanged.
- **Direct test:** VEGAS adapting to the slice 1e-5 < τ₂ < 3e-5 of lo (option
  `vslice`, 11th argument).
  - Iterations 1–6 (four seeds): 38 … 197; the 6-iteration result is
    158.6 ± 9.1.
  - With 12 iterations, iterations 7–12 give 169.3 ± 1.6, against DISENT's
    171.9.
  - The early iterations are strongly biased low. Even late ones fluctuate
    107–272.
- **Conclusion.** The near-2+1 region of the 3+1 phase space (relative
  measure ~1e-5 under nlo31's sequential decays) is poorly mapped. VEGAS's
  inverse-variance combination of iterations is biased low there, worst for
  r (4+1, large cancellations), whose magnitude is underestimated at small
  τ_cut, so the NNLO sum rises.
  - The below-cut pieces pass every check (degree-4 polynomial structure,
    additivity); the dipoles and the real ME are fine.
  - The plateau 2e-2 … 2e-4 stands.
- **Fix (= the efficiency work):** a slicing-adapted phase space for mode 1,
  with 3+1 and 4+1 events as 2+1 Born × emissions sampled logarithmically
  (DISENT style, multichannel over the dipole-like maps). Then redo the
  τ_cut test.

**Slicing-adapted phase space `psmc` for the above-cut parts (5 Oct).**
- `dis31/psmc.f90`, option `psmc` of nlo31 mode 1 (design in
  `docs/nnlo21-plan.md`): flat sequential-decay channel plus Catani–Seymour
  emission channels (FF, FI, IF for all labels; 12 for 3+1, 30 for 4+1 on top
  of the full 3+1 mixture). y and 1−x are sampled logarithmically down to 1e-10,
  z̃/u with a logistic map. The weight is the inverse of the mixture density,
  with every channel density from the exact inverse CS map.
- Checks:
  - phase-space volume against the flat generator, n = 3 and 4, within 1–2σ
    above the technical cut (`disorder-comparisons/nnlo21/pstest`);
  - lo slice 1e-5 < τ₂ < 3e-5 (vslice): 172.08 ± 0.12 against DISENT's 171.9
    (uniform: 169.3 ± 1.6 from late iterations only); all lo cells agree with
    DISENT × 1/x within 1–1.5σ;
  - mode 0 unchanged (bit-identical).
- r, all bins, four seeds of 200k × 4 (`runs/rvar`), against the uniform tcut2
  (30 × 1M × 6):

  | τ_cut | uniform (tcut2) | psmc |
  |---|---|---|
  | 2e-2 | 18.9 ± 1.3 | 18.3 ± 4.6 |
  | 1e-3 | −156 ± 11 | −136 ± 19 |
  | 2e-4 | −619 ± 14 | −636 ± 34 |
  | 1e-4 | −929 ± 34 | −957 ± 46 |
  | 3e-5 | −1704 ± 57 | −1935 ± 49 |
  | 1e-5 | −2474 ± 83 | −3270 ± 64 |

  Agreement on the plateau; below 1e-4 the uniform r was too small in
  magnitude by −28, −231, −796, most of the NNLO drift (81, 342, 1100). This
  confirms the VEGAS-bias diagnosis. psmc is also about 1.6 times faster per
  point (fewer ME calls in the cut-away region), and iterations are stable
  (χ²/it 0.3).
- Full rerun of lo, vi, kp, r with psmc: `runs/tcut5` (60 r, 12 vi, 6 lo, 4 kp).

**NNLOJET dijet NNLO reference: set-up (5 Oct).** AK chose a ZEUS-like set-up;
the NNLOJET warmups start once the τ_cut test (tcut5) passes.
- Runcard `nnlo21/validation/nnlojet_epLJJ_zeus2j.run`: cuts of NNLOJET's
  ZEUS dijet example (1703.05977 / ZEUS 1010.6167: 125 < Q² < 20000 GeV²,
  0.2 < y < 0.6, Breit-frame jets_et > 8 GeV, lab −1 < η < 2.5, ≥ 2 jets,
  m₁₂ > 20 GeV), with photon exchange only, α = 1/137, E-scheme
  recombination (V4, not ZEUS's E_T scheme, so that our jet code matches),
  μ_R = μ_F = Q.
- NNLOJET's E_T (`v1_et`, ObsHelper.f90) is sqrt(p_T² + m²) of the jet, not
  E p_T/|p| as the first version of the runcard comment said (corrected);
  jets are ordered by Breit-frame p_T.
- NNLOJET's order of the cuts (`driver/core/ecuts.f`, `ecuts_dis`): kt
  clustering and jets_et in the Breit frame, then jets outside the lab η
  window are dropped, then njets and m₁₂ of the two leading jets.
- Parse test: LO 102.2 ± 2.0 pb (tiny run). Warmups prepared (not started) in
  `disorder-comparisons/nnlo21/nnlojet_zeus/warmup`: LO, V, VV 1M[5];
  R, RV 2M[5]; RR 4M[5].
- Our side still needs: sliced21 integrated over (x, Q²) with the jet
  selection on the projected 2+1 Born jets; nlo31 mode 0 with ≥ 2 jets,
  τ₂ > τ_cut, psmc, lab-frame boost for the η cut; a τ_cut scan in this set-up.

**psmc r: catastrophic weights (5 Oct, late morning).** In `runs/tcut5`, 38 of
the 60 r jobs (1M × 6 each) have an iteration with a point of weight ~1e150
or more (e.g. seed 815: −2.3e150 ± Inf in iteration 3). VEGAS's grid is then
destroyed and every later iteration is zero. lo, vi and kp (3+1) are clean.
- The psmc 4+1 phase-space weight alone is bounded under uniform random
  numbers (2e6 points with s_min > 1e-9 W²: all below 2.6e4, none non-finite;
  `pstest/wscan.f90`). So the huge values come from the integrand (real minus
  dipoles) at extreme points that the adapted VEGAS grid reaches, or from the
  combination with VEGAS's own jacobian.
- Consequence: the psmc r numbers above (`runs/rvar`, 4 seeds without
  blow-ups) are **not yet trustworthy**. A rare wrong or unbounded weight
  biases every seed, not only those where it shows. The comparison with the
  uniform r (most of the drift removed) stands only as a hint until the
  cause is found.
- Debugging: seed 815 rerun with a dump of the first point with |weight| >
  1e12.
- **Cause found:** the first spike of seed 815 (|weight| 1.4e12, s_min/W² =
  1.25e-9, partons 2, 4, 5 mutually collinear, τ₂(real) = 1.9e-9, cut) comes from
  the FF dipole (i, j; k) = (2, 5; 4), whose spectator is collinear with the
  emitter pair. p̃₂₅ = p₂ + p₅ − y/(1−y) p₄ nearly cancels (E ≈ 2.637 − 2.636
  GeV), so the mapped Born has an almost zero-energy parton. In floating point
  its invariants turn negative (2p̃·p̃/W² = −1.2e-8), a light-cone energy in
  `tau2cm` vanishes, every partition fails and T stayed at `huge`. `accept`
  then let this dipole pass every τ_cut while the real was cut: an unmatched
  dipole of 1e24–1e26 (the real is 1e17–1e19). VEGAS then adapts to the spike
  and the following iterations reach 1e150.
- This is a robustness bug of the τ₂ measure (`tau2cm`), not of psmc or the
  dipoles. psmc only reaches such configurations; the uniform map practically
  never did. Fix: a non-finite or `huge` T means a degenerate configuration
  and gives T = 0 (exactly, T is tiny there). The same guard in
  `slicing/mod_tau2_run.f90` (`tau2_jetframe`, used with DISENT's mapped
  kinematics). Results change only at such failed points.
- Testing: seed 815 rerun with the fix (dump of any |weight| > 1e12).

**ZEUS-like dijets, first pieces (5 Oct).**
- `nlo31` mode 2 (above the cut, ZEUS selection, 15 observable bins × 10
  τ_cut) and `sliced21` mode 2 (below the cut, integrated over Q² and y; beam
  tables on a grid in Q, `lp21_grid_build`/`lp21_grid_load`, Δln Q = 0.1,
  h = 0.02 in ln(ξ/(1−ξ)), 29 nodes, ~1.4 CPU-min each). Modes 0 and 1 checked
  bit-identical after the changes. Combination: `nnlo21/combine_zeus.py`.
- Grid against direct evaluation: at most 3e-4 (up), 7e-4 (down), 3e-3 (gluon)
  of the largest coefficient. A grid twice as fine in both directions changes
  b1 by ≤ 4e-5 relative and b2 by ≤ 0.0024 pb at τ_cut = 1e-5 (b2 = 1.04e4 pb)
  with identical random numbers: negligible.
- LO (sliced21 b0, 8 seeds) against NNLOJET LO: total 103.29 ± 0.03 against
  103.17 ± 0.04 pb (+0.1%, 2.5σ); bins within 0–2.5σ, ours systematically
  ~0.1% higher. NNLOJET's single-run errors are not reliable (its two LO runs
  differ by 3σ in Q² 1000–2000), so several NNLOJET seeds are needed before
  reading anything into 0.1%. The selection (lab η direction, E_T definition,
  cut order) is right: a flipped η window would change bins by tens of %.

**NNLO τ_cut test with psmc and the τ₂ fix (5 Oct, afternoon; preliminary).**
`runs/tcut6`: r (60 seeds, same seeds as tcut5) with the fixed `tau2cm`; no
blow-ups, iterations stable. The seed-815 rerun with the fix: iterations
−2941, −3188, −3276, −3300 (before: −7.7e4 ± 7.4e4, then −2.3e150), no point
with |weight| > 1e12. lo and kp rerun with the fixed binary: bit-identical to
tcut5 (the bug never hit the 3+1 parts); vi rerunning.
Combination (r tcut6, lo/vi/kp tcut5, b1/b2 tcut1), all τ_zQ bins:

| τ_cut | NLO b1+lo (DISENT 16.885 ± 0.04) | NNLO b2+vi+kp+r |
|---|---|---|
| 2e-2 | 6.56 | −3.3 ± 0.2 |
| 1e-2 | 7.44 | 2.7 ± 0.3 |
| 5e-3 | 9.18 | 8.8 ± 0.5 |
| 2e-3 | 13.12 | 11.9 ± 0.7 |
| 1e-3 | 15.00 | 9.9 ± 0.9 |
| 5e-4 | 15.87 | 8.2 ± 1.2 |
| 2e-4 | 16.60 | 7.9 ± 1.5 |
| 1e-4 | 16.88 ± 0.15 | 8.2 ± 1.7 |
| 3e-5 | 17.14 ± 0.11 | 9.0 ± 2.3 |
| 1e-5 | 17.12 ± 0.13 | 0.4 ± 2.9 |

- The drift (81, 342, 1100 at 1e-4, 3e-5, 1e-5 with the uniform r) is gone.
  NNLO is flat at ≈ 8.5 ± 1 from 1e-3 to 3e-5, with power corrections above
  2e-3. The 1e-5 point is 2.8σ low; the technical cut moved r at 1e-5 in the
  uniform runs, so `runs/tcut7` repeats r with s_min > 1e-11 W² (30 seeds).
- **Correction:** the earlier statement "NNLO near zero and flat from 2e-2 to
  ≈ 2e-4" (4 Oct, with the uniform r) is superseded. That plateau came from
  the biased uniform r with large errors (±8–15). With psmc the NNLO
  coefficient in all bins is ≈ 8.5 pb/GeV² (NLO coefficient 16.9), and the
  values at 2e-2…5e-3 are power-correction dominated.
- vi rerun with the fixed binary: bit-identical to tcut5 (like lo, kp). The
  preliminary combination above is therefore final for lo, vi, kp, r (60).
- Technical cut 1e-11 instead of 1e-9 (`runs/tcut7`, 18 of 30 seeds, the rest
  stopped): r errors grow ~100-fold (±244 against ±2.5 at τ_cut = 1e-5); paired
  shifts +465 ± 245 (1e-5), +274 ± 144 (3e-5), and −4 ± 2 even at 2e-2 where
  the cut cannot matter. The region 1e-11…1e-9 is rounding-dominated (as with
  logmap, 4 Oct); 1e-9 stays the working cut. Not decisive for the 1e-5 point.
- r at 1e-5 over the 60 seeds is Gaussian (mean = median, halves agree,
  bootstrap error = quoted). b2 has only 4 seeds. More statistics:
  `runs/tcut8` (60 more r seeds, 12 b1, 12 b2).
- Dispatcher: sizing from the 1-minute load average overfilled thserv06 twice
  (load 38); excess jobs requeued by PID. To change (after this dispatcher
  exits): use the instantaneous number of running processes.
- **Final (tcut8 added: r 120 seeds, b1/b2 16 seeds), all bins:** NNLO 10.1 ±
  0.7, 8.6 ± 0.9, 8.1 ± 1.1, 7.95 ± 1.4, 8.1 ± 1.9, 0.9 ± 2.4 at 1e-3, 5e-4,
  2e-4, 1e-4, 3e-5, 1e-5 (LO 61.1, so +13% of LO); NLO b1+lo 16.85 ± 0.17 …
  17.07 ± 0.19 against DISENT 16.885. Plateau 5e-4 … 3e-5. The 1e-5 point stays
  ~3σ low with doubled statistics; r is Gaussian there. Next: technical cut
  1e-10 and 3e-9. Page v3: https://claude.ai/artifact/ChKctwxDdEau7H9DGMxwyW
  (plot `nnlo21/plot_tcut.py compare`).

**psmc edge and the technical-cut tests (5 Oct, evening).**
- **Correction:** the 10⁻¹¹ technical-cut run (`runs/tcut7`) was attributed to
  rounding ("errors ×100, rounding-dominated"). But psmc's log maps stop at
  10⁻¹⁰ (y, 1−x, both ends of z, u), so the region 10⁻¹¹ … 10⁻⁹ W² was mostly
  populated only by the flat channel, with large weights. Undersampling explains
  the errors as well as rounding does; the run cannot tell them apart. (AK's
  question whether 10⁻¹⁰ is low enough.) The 10⁻¹⁰ half of `runs/tcut9` sits at
  the edge as well; the 3·10⁻⁹ half is unaffected.
- With the working cut 10⁻⁹ the edge is below what passes the cut (all pair
  invariants ~ y, yz, y(1−z) × a hard scale must exceed 10⁻⁹ W²), and the flat
  channel keeps the density positive everywhere, so no bias.
- `psmc_set_edge(e)`, nlo31 12th argument (default 10⁻¹⁰, bit-identical; explicit
  10⁻¹⁰ also identical). Volume test with edge 10⁻¹²: mix against flat within
  1.5σ (n = 3), 1σ (n = 4).
- `runs/tcut10`: edge 10⁻¹² with technical cut 10⁻⁹ (control), 10⁻¹⁰, 10⁻¹¹; 30
  seeds each (801–830, as tcut6).

**ZEUS-like dijets at NLO: DISENT references disagree (5 Oct, evening).**
- `slicing/dis2j_nlo.f90`: DISENTFULL driven directly (as tau2_nlo), ZEUS
  selection of nlo31/sliced21 mode 2. DISENT's weights include 1/NEV and pb
  (NRM): the cross section is the sum over events (first version divided by
  NEV and multiplied by GeV→pb again; corrected). LO 103.2–103.5 ± 0.13
  (4 seeds × 10M), = NNLOJET 103.17 and sliced21 b0 103.29. NLO coefficient:
  −59, −54, −66, −61 (±3–4) total; m12 45–65: −49.5 ± 0.6 (LO 13.1).
- disorder itself (`analysis/zeus_dijet_analysis.f`, FastJet-free, same
  selection; `disorder -p2b -nlo|-nnlo -Q2min 125 -Q2max 20000 -ymin 0.2
  -ymax 0.6 -Ehad 920`, 4 seeds × 10M): LO (−nlo) 103.35 ± 0.03, but LO + NLO
  (−nnlo) −559 ± 5 total; −4.9 in ptavg 30–60 (LO 2.66), −104 in m12 45–65.
- Our slicing (b1 + lo, mode 2): +9 … +10 at τ_cut ≤ 1e-4.
- So the two DISENT drivers differ from each other by ~10 at O(α_s²), and both
  are far from slicing; neither looks physical (−380% NLO in a bin). Not the
  lepton azimuth (DISENT fixes the lepton plane and generates the hadronic
  azimuth; reals are built from their 3-parton parents). Cause not found.
- Tie-breaker: NNLOJET NLO (R + V) production, `nnlojet_zeus/nlo`.

**ZEUS NLO resolved; DISENT analysis pitfall (5 Oct, evening).**
- AK: DISENT passes counter-events (and collinear terms) with partons that are
  exactly soft or exactly collinear, which can upset an analysis. That was it:
  both my kt routines checked a parton's beam distance and then its pairs, so a
  parton exactly collinear to the beam (p_T = 0: d_iB = 0, but also d_ij =
  min(p_T²) ΔR² = 0) was merged into a jet when it came second, unlike the
  nearby real events (finite p_T, d_ij = p_T² ΔR² ≫ d_iB). The subtraction
  stopped cancelling. Fix: all beam distances first (ties to the beam), pairs
  with a zero-p_T parton skipped. The disorder analysis also had a stale
  parton→jet map (a fourth momentum added to three-parton jets). anti-k_t
  (nnlo11's obs11) is safe: it drops p_T = 0 partons explicitly.
- ZEUS-like dijets, NLO coefficient [pb] (LO 103.29):
  | | total | m12 30–45 | ptavg 15–22 | Q² 125–250 |
  |---|---|---|---|---|
  | NNLOJET (R 20 + V 10 seeds) | 10.49 ± 0.07 | 6.76 ± 0.06 | 5.07 ± 0.02 | 7.46 ± 0.06 |
  | DISENT driver (fixed) | 10.23 ± 0.22 | 6.94 ± 0.15 | 5.00 ± 0.04 | 7.31 ± 0.21 |
  | disorder −nnlo minus −nlo | 11.03 ± 0.18 | 7.17 ± 0.07 | 5.10 ± 0.04 | 7.78 ± 0.12 |
  | slicing b1 + lo, τ_cut 1e-4 | 9.8 ± 0.5 | 6.84 | 5.11 | 7.85 |
  NNLOJET LO (4 seeds) 103.29 ± 0.05 = sliced21 b0 103.29 ± 0.03: the earlier
  0.1% tension came from NNLOJET's single-run error.
- **Correction:** I first attributed disorder's +0.8 pb to its α_s (3-loop at
  −nnlo, 2-loop at −nlo). Numerically that is only ≈0.1% of the LO (≈0.1 pb):
  the difference stays unexplained (2.8σ, small; open).
- Higher-statistics slicing NLO in this set-up: `runs/znlo2` (32 b1 seeds with
  the soft table, 12× faster, agrees with direct to 5e-8; 32 lo seeds).

**NLO/NNLO 1+1 for lab-frame observables against disorder (AK, 5 Oct).**
- `nnlo11 -integrated`: DISENT over 125 < Q² < 20000, 0.2 < y < 0.6; τ₁ with
  each event's x; pure τ₁ slicing with LP_k(τ_cut; x, Q)/Born from a table
  (dis_tau1_lp on 237 × 57 nodes in ln(x/(1−x)), ln Q, step 0.05; interpolation
  error ≤ 3e-4 (c1), ≤ 1e-3 (c2, at τ_cut 0.1), typically 1e-5…1e-4); absolute
  bins (leading-jet p_T 5…100 GeV, rapidity −1…2.5). Linked with a copy of the
  current DISENT that gives the Born to USER (sed of the commented lines).
  Fixed-point mode bit-identical to the 2 Oct production binary.
- Reference: `analysis/lab11_analysis.f` (obs11's definitions: lab frame,
  proton +z, partons with E ≤ 0 or p_T = 0 dropped, anti-k_t R = 1 with
  rapidity, p_T > 5, −1 < y < 2.5), `disorder -p2b -nlocoef / -nnlocoef`.
- α_s: disorder's own coupling (nf = 5, 2-/3-loop from α_s(M_Z)) differs from
  LHAPDF's by ≤ 0.08% (O(α_s)) and ≤ 0.44% (O(α_s²), at Q = 11 GeV).
- Runs: `runs/lab11` (48 slicing seeds × 20M, DISENT cutoff 1e-10; 16 + 16
  disorder seeds × 10M). Combination: `slicing/lab11_combine.py`.
- **ZEUS NLO with statistics (`runs/znlo2`, 32 + 32 seeds):** our NLO
  coefficient converges to NNLOJET: total 12.52, 11.23, 10.95, 10.78 (± ≈ 0.27)
  at τ_cut 1e-3, 1e-4, 3e-5, 1e-5 against 10.49 ± 0.07; at 1e-5 all 15 bins
  within 1.5σ. The excess falls ≈ 2.6× per decade (≈ √τ): fiducial power
  corrections of the recoil-free projection, concentrated next to the cuts
  (p̄_T 8–15 next to E_T > 8, m12 20–30 next to m12 > 20, low Q²); bins away
  from the cuts agree from τ_cut ≲ 1e-3. For the NNLO in this set-up: τ_cut ≈
  1e-5 or a recoil-aware projection. Page v5/v6.
- **1+1 lab-frame validation (`runs/lab11`):** τ₁ slicing against disorder,
  pulls over 16 bins: NLO ≤ 1.2σ at 1e-5 (a few 2.7–2.9σ at 3e-4); NNLO mostly
  ≤ 2σ from 1e-3 down; the NNLO total (inclusive, exact) +0.2σ at 1e-5. Open:
  ≥ 2 jets at NNLO −2.9σ at every τ_cut (both sides DISENT NLO 2+1 there:
  DISENT cutoff 1e-10 vs disorder's 1e-8, and α_s ≈ −1.4σ); disorder with
  cutoff 1e-10 running (`runs/lab11b`). Plot: `plot_tcut.py lab11`.
- **≥ 2 jets resolved:** disorder with DISENT cutoff 1e-10 (`runs/lab11b`):
  67.25 ± 0.23 (1e-8: 67.69 ± 0.21) against slicing 67.04 ± 0.10: −0.9σ. The
  NNLO total moves by −0.007 (3σ of its tiny error) with the cutoff.
- **Strength of the 1+1 NNLO test (qualifies the statement above):** the
  slicing error on the NNLO total is ±5.7, 12, 32, 71, 190, 310 pb at τ_cut
  3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5 (coefficient −42 pb; disorder ±0.002). So
  at NNLO the agreement holds at the 10–30% level of the coefficient in most
  bins (tighter in the forward bins, e.g. y 1–1.5 ±0.1 at 1e-3); the pulls
  below 1e-4 carry little information. NLO is a sharp test and passes. A sharp
  NNLO 1+1 test needs cluster statistics (or variance reduction): the same
  cost problem of pure slicing as at the fixed point (2 Oct).

**The 1e-5 point: psmc's sampling edge (5 Oct, evening).** `runs/tcut9`,
`runs/tcut10` (30 r seeds per variant, seeds of tcut6; NNLO with that r
swapped in, all bins):

| technical cut, psmc edge | 2e-4 | 1e-4 | 3e-5 | 1e-5 |
|---|---|---|---|---|
| 1e-9, 1e-10 (tcut6/8) | 8.1 ± 1.1 | 8.0 ± 1.4 | 8.1 ± 1.9 | 0.9 ± 2.4 |
| 1e-9, 1e-12 | 7.6 ± 3.1 | 5.9 ± 3.4 | 9.9 ± 4.7 | 7.1 ± 5.2 |
| 1e-10, 1e-12 | 10.5 ± 4.1 | 10.5 ± 4.1 | 11.1 ± 6.0 | 17.0 ± 8.2 |
| 3e-9, 1e-10 | 8.0 ± 1.6 | 8.7 ± 1.9 | 6.7 ± 2.8 | −8.3 ± 3.1 |
| 1e-10, 1e-10 | 22 ± 12 | 31 ± 19 | 49 ± 34 | 73 ± 62 |
| 1e-11, 1e-12 | 79 ± 35 | 131 ± 57 | 257 ± 115 | 593 ± 286 |

- With the log-map edge at 1e-12 the 1e-5 point is on the plateau (7.1 ± 5.2):
  the deficit was undersampling at the edge 1e-10 (AK's question whether 1e-10
  is low enough: it is not, for τ_cut = 1e-5). A higher technical cut (3e-9)
  makes it worse (−8.3 ± 3.1); 1e-10 is consistent but noisier; with an edge
  above the cut (1e-10, 1e-10) the flat channel alone covers the gap, with
  large weights.
- **Correction of today's correction:** at technical cut 1e-11 the results are
  noise-dominated even with the edge at 1e-12, so for that run rounding in R −
  dipoles was the right explanation after all (tcut7).
- Working set-up for small τ_cut: edge 1e-12, technical cut 1e-9. The NNLO
  τ_cut test passes from 5e-4 to 1e-5.

**Correction: what the 1+1 validation was meant to be (AK, 5 Oct evening).** I
ran pure τ₁ slicing for the (N)NLO 1+1 (`runs/lab11`). AK wanted the
one-order-lower analogue of the N3LO plan: P2B on our own τ₂-sliced NLO 2+1.
The τ₁ results stand as a separate check, but they are not the requested test.

**NNLO 1+1 by P2B + τ₂ slicing (requested test).**
- NNLO 1+1 (lab-frame jets, `lab11` bins) = inclusive NNLO structure function ×
  O(Born) (disorder in inclusive mode, `-nnlocoef`, same analysis) + Σ_events w
  [O(event) − O(1+1 Born at the event's x, Q², y)] over our NLO 2+1 = `sliced21
  b1` (below τ₂ cut) + `nlo31 lo` (above), both in the new mode 3 (P2B, lab
  frame from a tetrad of Breit-frame vectors; projected Born q + xP).
  Reference: `disorder -p2b -nnlocoef` (DISENT's NLO 2+1 in place of ours).
- First check without slicing, NLO 1+1 = disorder inclusive `-nlocoef` +
  `sliced21 b0` mode 3 (6 seeds) against `disorder -p2b -nlocoef`: total exact;
  13 of 16 bins within 2σ; ≥ 1 jet +2.6σ, y −1…−0.5 +2.4σ, y 1.5…2.5 −3.4σ
  (0.4%, P2B part only). Being checked: more b0 seeds and disorder with DISENT
  cutoff 1e-10 (`runs/lab11c`).
- Production: `runs/p2b2` (48 b1 seeds × 3M×6 with soft table; 48 lo seeds ×
  20M×6 with psmc edge 1e-12; 16 disorder inclusive seeds).
- Aside (AK, after talking to K. Melnikov): concern about the efficiency of the
  one-loop 2→4+H (our NLO 3+1, vi). Per point (thserv, one core): vi ≈ 2 ms,
  r ≈ 1.4 ms (flat; psmc ≈ 1.6× faster); NNLOJET RV > 6 ms, RR > 3 ms (first
  warmup iterations). A proper CPU × error² comparison is part of the
  efficiency study.

**NNLO 1+1 by P2B + τ₂ slicing: first results and efficiency (5 Oct, ~21:45).**
- NLO 1+1 with high statistics (64 b0 × 45M; disorder 64 × 75M, cutoff 1e-10):
  agreement below 0.01% of the jet rate in every bin (total −92.4060 ± 0.0014
  vs −92.4058 ± 0.0007 pb). (AK: per-mille precision on the NLO coefficient;
  enough here, more on the cluster.)
- NNLO 1+1 (11 b1 + 48 lo seeds, disorder 16 × 10M): consistent within errors
  at τ_cut 1e-4 … 1e-5; errors up to ±2% of the jet rate next to the Born p_T.
  More statistics running (b1 144, lo 288 seeds; disorder 128 × 40M).
- Efficiency, CPU × error² at τ_cut 1e-4 (thserv jobs, optimal b1/lo split):
  ours / disorder = 2–7 (p_T bins around the Born jet), 9–12 (p_T 20–50),
  30–100 (≥ 1 jet, central/backward y), 100–170 (forward y, ≥ 2 jets, p_T
  5–8). Cause: b1 and lo each grow like ln² τ_cut and are sampled
  independently. AK: after the validation, try correlated sampling (the
  above-cut 3+1 and its below-cut 2+1 Born from the same point).
- **Full statistics (21:50; 144 b1 + 288 lo seeds; disorder 122 × 40M, cutoff
  1e-10):** NNLO 1+1 by P2B + τ₂ slicing against disorder + DISENT, χ² over the
  15 jet bins: 1822, 83, 28, 9.3, 15.1, 7.4, 10.2 at τ_cut 5e-3, 1e-3, 5e-4,
  2e-4, 1e-4, 3e-5, 1e-5 (largest pulls 39, 6.8, 3.1, 1.7, 2.1, 1.6, 2.0). So
  agreement for τ_cut ≤ 2e-4; power corrections above (mostly ≥ 2 jets and p_T
  5–8). Errors at 1e-4 ≤ 0.22% of the jet rate per bin (disorder ≤ 0.13%).
- Efficiency, CPU × error², ours/disorder: τ_cut 1e-4: 2–5 (p_T 8–20), 11–28
  (≥ 1 jet, p_T 20–100, central/backward y), 52–163 (forward y, ≥ 2 jets, p_T
  5–8); τ_cut 1e-5: 17–396. Page v11.

**Correlated sampling (5 Oct, night; AK: after the validation).**
- `sliced21` part `c1` (modes 2, 3): the NLO 2+1 coefficient b1 + lo sampled
  together. psmc builds every CS-channel 3+1 point from a 2+1 Born (unit numbers
  r(4:6)); `psmc_born` returns that Born and its density, b1 is evaluated there
  with weight 1/g_Born, lo (nlo31 `born_part`) at the 3+1 point from the same r.
  Each term unbiased; the Born-to-Born fluctuation of the large logs should
  cancel. Options: `C1_NEMIT = M` (M emissions per Born, fresh channel/emission
  numbers), `C1_NEMIT = -1` (stratified: one emission per psmc channel, 12 CS +
  flat, weighted by the channel probabilities).
- Regressions: modes 0–2 bit-identical after splitting b21_part into kinematics
  and `b21_eval`. Mode 3 VEGAS target changed from the total row (identically
  zero in P2B, so VEGAS never adapted there) to ≥ 1 jet.
- Correctness: c1 reproduces the independent b1 + lo in all bins (pulls ≤ 1.6).
- Efficiency, CPU × error² relative to the independent b1 + lo (optimal split;
  timing favours c1 by ~1.5–2× since thA371a is faster than a thserv):
  - VEGAS-adapted c1 (target ≥ 1 jet): 20× (1e-4) / 90× (1e-5) better in ≥ 1
    jet, much worse elsewhere: the adaptation starves the other bins (the
    independent production effectively did not adapt). Not a fair comparison.
  - unadapted c1, M = 1: 2–20× worse: b1 costs ~30× lo per point, so 1:1
    points starve lo (optimal independent split ~7 lo per b1).
  - unadapted c1, M = 8: about break-even (0.1–2 at 1e-4; up to 7 at 1e-5 in
    ≥ 1 jet). Not the hoped one to two orders of magnitude.
  - stratified (M = −1, one emission per channel; correct: pulls ≤ 1.3 except
    one 2.5 among 15): 0.1–5 at 1e-4 (most ≈ 1), 0.5–6.6 at 1e-5 (most 1–3).
    After the ~1.5–2× CPU bias: break-even at 1e-4, up to ~3× better at 1e-5.
- **Conclusion:** correlated sampling works and helps where the logs are
  largest (small τ_cut), but only by a factor of a few, not orders of
  magnitude. The remaining variance is not the Born-level fluctuation of the
  logs; plausibly the genuine spread of hard emissions and of the P2B difference
  O(event) − O(Born). Next candidates: VEGAS adaptation on a balanced target
  (sum of |bins|) instead of one bin; larger τ_cut with NLP corrections.

**ZEUS-like dijets at NNLO, our side (5–6 Oct night; preliminary, 405/428 jobs).**
`runs/znnlo` (r 283, vi 59, kp 16, b2 48 seeds; psmc edge 1e-12, technical
cut 1e-9). NNLO coefficient, total [pb]: 11.0 ± 3.1, 23.1 ± 4.9, 32.2 ± 6.4,
48.7 ± 10.1, 67.6 ± 15.1, 84.7 ± 27.8, 86.1 ± 43.8 at τ_cut 2e-3, 1e-3, 5e-4,
2e-4, 1e-4, 3e-5, 1e-5. **No plateau.** The rise sits in the bins next to the
cuts (p̄_T 8–15: 5.6 → 72.9; m12 20–30: 3.6 → 54.4; low Q²); bins away from the
cuts (p̄_T 22–60, m12 65–120) are small and flat within errors.
- Same pattern as the NLO fiducial power corrections (√τ), but much larger; at
  NNLO the recoil-free projection's fiducial power corrections can carry up to
  ln³ τ, so they need not be small at 1e-5.
- The apparent limit (≈ 90 pb, close to the LO) would be an implausibly large
  NNLO correction for HERA dijets (NNLOJET: a few % with μ² = (Q² + p_T²)/2).
  So neither the plateau nor the number is trusted. Possible causes: (a) huge
  fiducial power corrections, (b) a problem in mode 2 at NNLO (vi, kp, r in mode
  2 have not been validated against anything; mode 1 has, at the fixed point).
- What decides it: NNLOJET's NNLO for this set-up (cluster). Meanwhile this
  points to a recoil-aware projection (or NLP corrections) for jet observables
  with cuts. To discuss with AK.
- **Balanced VEGAS target (`C1_VBAL=1`, with stratified emissions):** VEGAS
  adapts to √(Σ_b f_b²) of the point's combined (b1 + lo) contribution over all
  jet bins at τ_cut = 1e-4. 12 seeds on thA371a: CPU × error² 10–40× better
  than the independent sampling in most bins at 1e-4 (25–65× at 1e-5); p_T
  30–100, forward y, ≥ 2 jets 1–5×. After the CPU bias still ~5–30×. Pulls
  against the independent b1 + lo mostly fine, but +3.1 (p_T 20–30 at 1e-4) and
  −2.3/−2.6 (y 1–1.5): possible bias from VEGAS's inverse-variance combination
  of iterations with an adapted target, to be checked. Running: 48 seeds on
  thserv09/15 (`runs/corrv`) for statistics and a timing on the same CPU as the
  disorder reference.
- **Balanced correlated sampling, 48 seeds on thserv09/15 (`runs/corrv`, all
  timings on thservs):**
  - τ_cut 1e-4: pulls against the independent b1 + lo ≤ 1.8σ (no sign of
    bias). Gain over the independent sampling 5–22× (≥ 1 jet, p_T 8–20, y −1…0),
    4–11× (y 0…1, p_T 20–30), 1.4–2× (p_T 30–100), 1.8–6× (p_T 5–8, forward y,
    ≥ 2 jets).
  - **Against disorder (CPU × error², ours/disorder):** 0.2–0.97 in the core
    bins (≥ 1 jet 0.63, p_T 8–11 0.29, 11–15 0.21, 15–20 0.36, y −1…−0.5 0.73,
    −0.5…0 0.97): P2B + τ₂ slicing with correlated sampling is as fast as or
    faster than disorder there. 1.8–2.8 (y 0…1, p_T 20–30), 5–6 (p_T 30–100),
    25–90 (p_T 5–8, y 1…2.5, ≥ 2 jets: bins filled by hard emissions far from
    the Born, where the correlation does not help).
  - τ_cut 1e-5: gains up to 80×, but 4 of 15 bins at 2.1–2.9σ against the
    independent result: test with equal-weight iterations (`VEGAS_EQUAL=3`,
    `runs/corrve`).
- **Correction (1e-5 tension):** the "4 of 15 bins at 2.1–2.9σ" came from
  comparing with the thserv subset of the independent runs only (52 of 288 lo
  seeds, used for the timing). Against all independent seeds the balanced
  correlated result has χ² = 12.9/15 (1e-4) and 17.0/15 (1e-5), largest pull
  2.0σ: no sign of bias. Equal-weight iterations (`VEGAS_EQUAL=3`, `runs/corrve`):
  also unbiased (χ² 17.2, 20.9/15) but the poorly adapted early iterations
  inflate some errors enormously (y −1…0: ±60 instead of ±2); keep 1/s².
- **Correlated NNLO 2+1 (`sliced21 c2`, mode 1, 6 Oct ~01:15):** one psmc Born
  feeds b2 (below the cut), vi and kp (3+1 from psmc_gen(3) with that Born) and
  r (4+1 from psmc_gen(4), whose 3+1 step takes the same Born numbers); 16 VEGAS
  dimensions, balanced target over the τ_zQ bins. ~2 ms per point (vi and r).
  Test: 48 seeds × 1.8M on thserv09/15 (`runs/corr2`) against the independent
  τ_cut runs (tcut5/6/8; compare at τ_cut ≥ 3e-5, the old r used edge 1e-10).
- **Correlated NNLO 2+1 (`c2`): does not work as implemented.** 48 seeds ×
  1.8M (thservs): CPU × error² ~10× *worse* than the independent parts, and a
  systematic excess: total 25.0 ± 6.9, 41 ± 10, 53 ± 12 at τ_cut 1e-3, 2e-4,
  1e-4 against the independent 10.1 ± 0.7, 8.1 ± 1.1, 8.0 ± 1.4 (pulls +2.2 …
  +3.8; per bin up to +4.4). Heavy positive tails (largest seeds 261, 375) but
  the median is shifted too (17 vs 10 at 1e-3, 29 vs 8 at 1e-4): a bias, not
  only tails. Candidates: VEGAS adaptation/iteration weighting on a heavy-tailed
  16-dimensional integrand with the balanced target, or a bug in how c2
  combines the parts. Next: c2 with one part at a time against the independent
  part, and unadapted. Not to be used until understood. (The NLO-level `c1`
  is unaffected: χ² 12.9/15 against the independent result.)
- **Mode-2 diagnostic for the ZEUS NNLO (`ZFIX=1`, `runs/zfix`):** Q² and y
  squeezed to a window of relative width 1e-3 around x = 0.01, Q² = 400 with
  the mode-1 observable; σ_window / (x w² Q²) must equal mode 1. b0: 61.1
  (= mode 1). Parts against mode 1 (tcut runs), pulls at τ_cut 2e-2, 5e-3,
  1e-3, 1e-4: b2 +1.2, +0.9, −0.4, −1.0; kp −1.1, −0.9, +0.6, +1.0; r +0.6,
  −0.1, +1.2, 0.0; vi −1.0, −2.3, −1.9, −2.4 (vi ~1–2% low, 3 seeds; even if
  real it moves the NNLO total by ~2 pb). So mode 2's own code paths (Q²/y
  integration, per-event psmc, beam grid in Q) reproduce mode 1: the ZEUS
  no-plateau is most likely the fiducial power corrections (a), not a mode-2
  bug (b). NNLOJET's NNLO remains the decisive check; a recoil-aware projection
  is the likely remedy.

**Fiducial power corrections: literature (6 Oct morning).** Notes in
`docs/literature-fiducial-power-corrections.md` (Ebert–Tackmann 2020;
Campbell–Neumann–Vita 2408.05265, P2B-improved slicing in MCFM; Alioli et al.
2504.11357, GENEVA incl. Z+jet with 1-jettiness). The remedy for our ZEUS
drift is P2B-improved slicing: below the cut add ∫ dσ_{3+1} (O − Õ) with a
projection to the 2+1 Born. New for us: the projection for jet final states in
DIS. Next: NLO test in the ZEUS set-up.

**P2B-improved τ₂ slicing at NLO, ZEUS-like dijets (6 Oct morning).**
- What: CNV eq. 2.14 in `nlo31` mode 2 (env `P2BSLICE=1`): below τ_cut the
  3+1 events contribute O(event) − O(projected 2+1 Born); b1 unchanged.
- Projection (`project21`), new:
  - from the τ₂-minimising partition (`tau2cm` now records it in `t2code`):
    K = J₁ + J₂, K̃ = q + c·p_in with K̃² = K² (x, Q², y unchanged);
  - the jets go through the CS Lorentz map K → K̃ and are made massless back
    to back in the K̃ frame along ΛJ₁;
  - q has no transverse component in the Breit frame, so beam recoil changes
    K² only at O(k_T²).
- The code is generic in the number of partons, so vi, kp and r (each real
  and dipole configuration with its own partition) use the same treatment.
- Regression: without `P2BSLICE` the output (lo and r) is bit-identical to
  the 5 Oct binaries.
- Runs: `runs/zp2b`, 32 lo seeds (3M × 6, desktops), combined with the 32 b1
  seeds of `runs/znlo2`.
- Result: the NLO coefficient against NNLOJET (10.49 ± 0.07 pb), as % of LO:
  - total, plain → P2B: +2.77 → +0.58 (2e-3), +1.96 → +0.35 (1e-3),
    +1.40 → +0.20 (5e-4), +0.91 → −0.01 (2e-4), +0.28 → +0.10 (1e-5);
  - every bin agrees from τ_cut ≤ 5e-4 (plain needs about 3e-5).
- Fits with the full seed covariance (τ_cut ≤ 2e-3, 7 points):
  - plain is pure √τ (χ² 0.6/6; τ ln τ χ² 30): the 5 Oct "fiducial" reading
    is confirmed;
  - with P2B the √τ coefficient is 5.5× smaller (64 → 12 pb) and τ ln τ fits
    equally well (χ² 6.6/6); the remainder is largest at low p̄_T and m₁₂.
- **Correction:** I first reported (to AK, same morning) that P2B removed only
  ≈ 40% and left a √τ residual. That combination had averaged the plain and
  the P2B lo seeds ("lo 64 seeds" in the combine output). With the P2B seeds
  only, the numbers are as above.
- Plot `plot_tcut.py zeusp2b`; fit script `runs/zp2b/fit_pc.py`; page v17,
  section 6a.
- NNLO with P2B running: `runs/znnlop2b` (r 300, vi 64, kp 16 seeds,
  statistics of `runs/znnlo`; b2 reused).

**P2B slicing at NNLO: r blew up, fixed (6 Oct, 08:00–08:30).**
- First P2B NNLO run (`runs/znnlop2b`): vi and kp fine, but r reached
  1e14…1e18 in 290 of 300 jobs. The desktop jobs died silently, most likely
  from NaN in VEGAS after the overflow; the thserv jobs were killed by PID.
- Diagnosis (debug env `P2BDEBUG=1`, 200k-point runs on thA371a):
  1. `project21` failed for many mapped dipole Borns: either no partition
     (degenerate, T set to 0) or K² ≤ 0, c ≤ 0, or a degenerate jet. Õ := O
     there (weight 0) left the real's O − Õ unmatched.
     Fix: fallback partition, with the two most energetic partons as jets
     and the rest in the beam (the IR-consistent limit).
  2. Still blowing up: dipoles of reals deep in double-unresolved regions
     (τ₂(real) ≈ 1e-12) have numerically garbage mapped Borns (partons of
     ±1e5 GeV that cancel, dipole values 1e17–1e28, against a real of 1e-4).
     Plain slicing never saw them (real cut, those dipoles T = 0); with P2B
     they enter with O − Õ ≠ 0.
     Fix: with P2B, drop the whole event (real and its dipoles) if
     τ₂(real) < techcut (1e-9). The dropped piece is O − Õ over a region of
     size 1e-9, which is power suppressed.
- After the fixes, 6 seeds × 200k: stable (−8000 ± 500). Against the plain r
  (299 seeds): +160 at 1e-3, +210 at 5e-4, +140 at 1e-4 (P2B corrections of
  plausible size).
- Regression: without `P2BSLICE`, r is bit-identical to the 5 Oct binary.
  With P2B, lo, vi and kp are bit-identical to the first P2B binary (the
  fallback does not trigger there), so the NLO result (`runs/zp2b`) and the
  vi and kp runs stand.
- r rerun with the fixed binary: `runs/znnlop2b_r` (300 seeds, thservs).

**P2B r at full statistics: second fix (6 Oct, 10:20–11:20).**
- The 08:30 fix was incomplete. In `runs/znnlop2b_r` (300 × 1M × 6), about
  60% of the jobs still reached 1e10 or more in some iteration (a few up to
  1e96); the remaining jobs were killed by PID.
- Debug runs at 1M points: the same class of numerically garbage mapped
  dipole Borns (partons of ±5e6 GeV that cancel, dipoles up to 1e29). They
  also occur with τ₂(real) up to ≈ 2e-8, so the real-τ cut alone does not
  catch them. The 200k-point tests at 08:30 were too small to see them.
- Fix: with P2B, an event is dropped (real and all its dipoles) if any mapped
  dipole Born has no τ₂ partition (the tau2cm guard) or a negative-energy
  parton; the drops are counted.
- Test, 8 seeds × 1M × 3 (including the worst seeds): stable. Dropped
  13–35 of ≈ 1.45M real events per job (≈ 1e-5).
- r per τ_cut against plain (299 seeds), P2B / plain:
  - 2e-2: 43 ± 5 / 20 ± 2
  - 1e-2: 11 ± 12 / −25
  - 1e-3: −666 ± 15 / −664
  - 1e-4: −2995 ± 35 / −2953
  - 1e-5: −8560 ± 44 / −8494
  So the P2B correction in r is small, visible only at large τ_cut.
- **Correction:** the "+160 at 1e-3, +210 at 5e-4, +140 at 1e-4" of the
  08:30 entry came from 6 short runs (200k × 2) and are not confirmed: at
  full statistics r agrees with plain there within ≈ 20 pb.
- Production r: `runs/znnlop2b_r2` (300 seeds). `runs/znnlop2b_r` is
  invalid (kept only for the record). The vi seed lost on thserv19 was
  rerun as `znnlop2b_r/s301`.

**ZEUS NNLO with P2B: full result and the τ_cut = 1e-5 drop (6 Oct, afternoon).**
- Combination: b2 48 + vi 64 + kp 16 + r 300 (`runs/znnlop2b_r2`, drop rule
  P2BDROP = 0); `runs/znnlop2b_r2/nnlo_p2b.json`.
- NNLO coefficient, total (pb), P2B / plain:
  - 2e-2: 11.1 ± 0.9 / −0.4 ± 1.8
  - 2e-3: 17.8 ± 1.6 / 10.7 ± 3.0
  - 2e-4: 39.0 ± 2.7 / 48.2 ± 9.5
  - 1e-4: 39.9 ± 3.4 / 66.5 ± 14
  - 3e-5: 33.7 ± 4.4 / 81.9 ± 26
  - 1e-5: 4.9 ± 5.4 / 81.7 ± 41
  P2B errors are 3–8× smaller.
- There is still no plateau, and the drop at 3e-5 and 1e-5 also appears in
  bins away from the cuts (p̄_T 30–60: 0.5 → −1.0 ± 0.4; m₁₂ 45–65:
  1.8 → −3.5 ± 1.4), where P2B and plain agree. The drop is ≈ 0.3% of b2 at
  1e-5 in every bin: a per-mille effect at small τ_cut common to both
  methods, made visible by the smaller P2B errors.
- Excluded:
  - The b2 beam grid: coarse and twice-finer grids agree to 2e-7 (5 Oct,
    `beamgrid/gchk_b2_*`), although `tabchkg` shows a pointwise coefficient
    deviation up to 2e-3.
  - The P2B terms of vi and kp: they level off towards 0 at small τ_cut.
- **Correction:** my first two diagnostics (dropped-event contribution
  P2BDROP = 2, and the P2B term P2BEXTRA = 1) were not reliable. VEGAS
  adapted to the 1e-5 column, where these integrands are ≈ 0, and weighted
  iterations by 1/variance of that column. On lo the P2BEXTRA diagnostic gave
  nonsense. With `VTARGET = 1` and `VEGAS_EQUAL = 2` it reproduces the true
  lo P2B term (−7.6 vs −7.2 at 2e-2 … −0.26 vs −0.19 at 1e-5). The conclusion
  "the dropped events contribute ≈ 0" is therefore not established; it is
  re-tested.
- Note: the lo P2B term (P2B minus plain lo) is −7.19 ± 0.10 pb at 2e-2,
  −1.67 at 1e-3 and −0.19 at 1e-5. The −3.6 etc. that I looked at on the
  morning of 6 Oct (exploratory output, not in this notebook) came from the
  mixed combination. The NLO results in the 07:00 entry are unaffected.
- Running:
  - P2B term of r (P2BEXTRA, 8 seeds);
  - psmc edge 1e-13 scan (`runs/zedge`: r 100, vi 16, kp 8);
  - next, the techcut slice [1e-10, 1e-9]·W² (TECHSLICE, new diagnostic).

**The τ_cut = 1e-5 drop is the technical cut (6 Oct, 16:00–18:00).**
- Slice diagnostic (`TECHSLICE`; VEGAS_EQUAL = 2; validated: with the slice
  opened to everything it reproduces production r within 1–2σ),
  `runs/zedge`. Events with 1e-10 ≤ s_min/W² < 1e-9, removed by the
  production cut, contribute to the NNLO total (ZEUS, P2B):
  - < 1 pb for τ_cut ≥ 1e-4, +7.6 ± 2.7 at 3e-5, +30.3 ± 5.6 at 1e-5;
  - positive in 18 of 20 seeds.
  - Mode 1 (x = 0.01, Q² = 400): +4.1 ± 2.2 (3e-5), +8.4 ± 1.7 pb/GeV² (1e-5).
- Next decade, 1e-11 … 1e-10:
  - ZEUS: +3.6 ± 1.7 (3e-5), +16.6 ± 4.3 (1e-5), heavy-tailed seeds;
  - mode 1: −0.5 ± 0.9, converged.
- Mechanism: near a single-unresolved limit the real can be just above
  τ_cut and a dipole just below. These are legitimate slicing contributions,
  in a window Δτ ~ s_min/Q², so the cut must satisfy techcut·W² ≪
  τ_cut·Q². With W²/Q² ≈ 1/x up to 500 in ZEUS, a cut relative to W² is too
  coarse at low x, and lowering it runs into rounding.
- **Correction of 5 Oct ("the NNLO τ_cut test passes from 5e-4 to 1e-5"):**
  - at the mode-1 point the 3e-5 and 1e-5 values with techcut 1e-9 are
    biased low (by ≈ 4 and ≈ 8 pb/GeV²);
  - the 5 Oct run with techcut 1e-10 (11.1 ± 6.0, 17.0 ± 8.2) already hinted
    at it within its errors;
  - the plateau 5e-4 … 1e-4 stands.
- Fix (`TECHMIS = 1`): in the technical region, keep per cell only the
  mismatch part Σ_i (F_real − F_i) D_i, which involves dipoles only, so
  nothing cancels and there is no rounding problem. Drop the cancelling part
  F_real (R − ΣD). Default bit-identical (plain, P2B, mode 1).
- Validation running (`runs/ztech`):
  - A: mismatch part in [1e-10, 1e-9] against the full slice (30.3 ± 5.6);
  - B: [1e-11, 1e-9] against the sum (46.9 ± 7);
  - C: mode 1 [1e-10, 1e-9] against 8.4 ± 1.7;
  - D: everything below 1e-9 down to the sampling edge, to add to the
    production r.
- **TECHMIS failed validation and was removed.**
  - A gave ≈ 3e6 pb (target 30 ± 6), and ≈ 600 pb even at τ_cut 2e-2,
    where the slice is 0; B and C were similar.
  - The split F_R (R − ΣD) + Σ (F_R − F_i) D_i is exact, but in the
    technical region both pieces are individually huge and opposite: several
    dipoles are singular in the same invariant and map to different Borns.
  - Jobs stopped by PID.
- **Plan B (`TECHQ = c`): technical cut s_min < c·Q².** The mismatch
  window scales with s_min/Q².
  - Mode 1 converged once the cut reached 1e-10 W² = 1e-8 Q².
  - Rounding sets in below ≈ 1e-11 W²; c = 1e-8 stays ≥ 2e-11 W² in ZEUS.
  - Diagnostic `TECHDIFF = 1`: weight θ(s_min ≥ c Q²) − θ(s_min ≥ reference).
    Default bit-identical.
- Runs in `runs/ztechq`:
  - V: mode 1, 9.9e-9 Q² (= 1e-10 W² there) against 1e-9 W²; must equal
    the measured slice 8.4 ± 1.7;
  - P: ZEUS P2B, 1e-8 Q² against 1e-9 W², 60 seeds; to add to the
    300-seed r;
  - C: ZEUS, 1e-9 Q² against 1e-8 Q²; convergence, must be ≈ 0.
- **Plan B (cut c·Q²) failed its check.**
  - V: −34.6 ± 3.1 at 1e-5 (expected ≈ +7.9).
  - C (1e-9 Q² against 1e-8 Q²): +66 ± 18.
  - For events with small partonic W² (= Q²(η/x − 1), ≲ 10 Q²) a Q² cut is
    coarser than the production cut and removes large legitimate
    contributions.
  - My own setup error, caught before the results: V's comparison region
    is not exactly the slice region, because the "W²" of the cut is the
    partonic W², which varies per event. The right target was the sum of
    the slices (7.9 ± 1.9), not 8.4.
- Where the removed contributions sit (slice per Q² bin, at 1e-5):
  - Q² 125–250: 17 ± 4 (first decade), 12 ± 4 (next), not converged;
  - Q² ≥ 500: next decade ≈ 0, converged.
  So low x, i.e. large W²/Q², needs a finer cut, and small W² must not get
  a coarser one.
- **Plan C:** cut s_min < min(cw W², cq Q²) = min(1e-10 W², 1e-8 Q²)
  (`TECHMIN`, reference `TECHREF`, `TECHDIFF`).
  - It is never coarser than 1e-10 W², and at most 2e-11 W² at the lowest x.
  - Exact machinery check: TECHDIFF over [1e-10, 1e-9] W² is bit-identical
    to the TECHSLICE run with the same seed. Default bit-identical.
- Runs `runs/ztechm`, expectations fixed beforehand:
  - V′: mode 1, ≈ 7.9 ± 1.9 (1e-5) and ≈ 4.0 ± 2.2 (3e-5);
  - P′: ZEUS correction; in the Q² ≥ 500 bins ≈ the first slice;
  - C′: min(1e-11 W², 1e-9 Q²) against the new cut, ≈ 0;
  - D2: drop-rule bias, i.e. the plain weights of events with
    τ_real < 1e-9 (my earlier measurement of it was flawed).
- **Plan C passes all checks** (`runs/ztechm`; expectations fixed beforehand):
  - V′ (mode 1): 8.7 ± 3.5 (1e-5), 3.9 ± 1.8 (3e-5), 2.4 ± 2.1 (1e-4),
    against 7.9 ± 1.9, 4.0 ± 2.2, ≈ 0.
  - P′ in the Q² ≥ 500 bins at 1e-5: 3.9 ± 0.7, 0.8 ± 0.2, 0.6 ± 0.1, 0.05,
    against the first slice (+ next decade) ≈ 2.4, 1.4, 0.9, 0.1.
  - C′ (min(1e-11 W², 1e-9 Q²) against the new cut): −0.1 ± 0.3 (1e-4),
    1.0 ± 3.8 (3e-5), 6.1 ± 3.6 (1e-5). Converged at 3e-5; at 1e-5 a
    residual up to ≈ 6 ± 4 pb is allowed.
  - D2 (drop-rule bias, 13 seeds): −2.3 ± 2.1 at 1e-5, 0 elsewhere.
- **Corrected ZEUS NNLO coefficient (P2B)**, = 300-seed P2B result + P′
  (60 seeds), `runs/ztechm/nnlo_p2b_corr.json`. Total (pb): 17.9 ± 1.6
  (2e-3), 23.9 ± 1.9, 30.9 ± 2.2, 40.1 ± 2.7 (2e-4), 41.7 ± 3.4, 42.4 ± 5.0,
  43.2 ± 6.1 (1e-5). Flat from 2e-4; all bins flat within errors from 2e-4.
  No NNLOJET NNLO reference yet.
- Plot `plot_tcut.py zeusnnlop2b`; page v18, section 6b, with corrections
  (section 2 note, invalid diagnostics, the failed technical-region ideas).
- Small open item: psmc edge 1e-13 against 1e-12 shifts vi by +2 ± 0.9 pb
  at small τ_cut (`runs/zedge`); r and kp unchanged.
- Next: production with the new cut built in (cluster instructions to
  update: TECHMIN = "1d-10 1d-8", P2BSLICE = 1), the mode-1 test redone with
  the new cut, and plain slicing with the new cut for comparison.

- Pushed `2026-10-nnlo21` with AK's OK (6 Oct, evening, up to 2697d78), for the cluster production.

**Fixed-point τ_cut test with the new technical cut (6 Oct, night).**
- `runs/tcut11`: r 120 seeds with `TECHMIN="1d-10 1d-8"` (= 1e-10 W² at
  this point, since W² ≤ 99 Q²), vi/kp/lo with psmc edge 1e-12.
- **The direct run is contaminated:** rare numerically garbage events,
  single seeds with r = −0 or +2108 instead of ≈ −3260 at 1e-5 (as in the
  1e-11 run of 5 Oct). Mean −3148 ± 60 against median −3261. The combined
  NNLO from it (14.8, 18.2, 29 ± 9, 44 ± 15 …) is not usable.
- The difference runs (`runs/ztechm` P′, C′, V′) and the default-cut
  productions are clean (mean ≈ median ≈ trimmed).
  → Procedure: production with the default cut plus a separate `TECHDIFF`
  correction (as for ZEUS). The cluster instructions said "TECHMIN
  directly"; **corrected**, not yet pushed.
- Fixed point, robust version (b1/b2 16 seeds, lo/vi/kp tcut11, default-cut
  r, plus the slice correction), all bins, pb/GeV²:
  - r tcut6/8 (edge 1e-10): 10.2 ± 0.7 (1e-3), 8.7 ± 1.0, 8.5 ± 1.2,
    8.1 ± 1.4 (1e-4), 13.1 ± 2.9 (3e-5);
  - r tcut10 control (edge 1e-12, 30 seeds): 6.0 ± 3.4 (1e-4),
    14.9 ± 5.2 (3e-5), 15.5 ± 5.5 (1e-5).
  The plateau ≈ 8.5 from 5e-4 to 1e-4 stands. At 3e-5 and 1e-5 the values
  are 1.3–1.5σ above it: consistent, with a slight upward tendency that
  needs the cluster statistics.

**Cluster results (7 Oct, `2026-10-nnlo21-cluster`, merged): our ZEUS NNLO
disagrees with NNLOJET.**
- NNLOJET (11,931 jobs): LO 103.2810 ± 0.0005, NLO coefficient 10.5797 ±
  0.0048, **NNLO coefficient 11.12 ± 1.19 pb**.
- Ours (projected slicing, default cut + TECHDIFF, ≈ 10× the MPP
  statistics): 11.3, 11.8, 12.6, 18.1, 24.0, 33.2 ± 1.7, 40.7 ± 1.8,
  48.1 ± 2.0, 49.3 ± 5.8, 50.1 ± 6.0 from τ_cut 2e-2 to 1e-5.
  - The excess at 1e-4 is 18–45% of LO in every bin, also far from the
    cuts. The NLO agrees with NNLOJET.
- **Corrections of my statements:**
  - "Flat from 2e-4 at ≈ 41 pb" (6 Oct, page v18) is wrong: the coefficient
    keeps rising and is ≈ 4× NNLOJET.
  - "psmc edge 1e-13: r unchanged" is wrong: with 300 seeds r moves by
    +10.8 ± 3.4 (1e-3), +22 ± 7 (1e-4), +44 ± 11 (1e-5). My 100 seeds had
    ±8 and I read "within errors" as "unchanged".
  - The fixed-point NNLO plateau has no independent reference; it only
    shows τ_cut independence. The cluster gets 10.3–11.5 there (5e-4…1e-4),
    against my 8.5 ± 1.2.
- Not affected: the NLO (P2B and plain) and the NNLO 1+1 by P2B (D:
  χ²/15 of 22–33 at τ_cut ≤ 2e-4, inclusive total −42.1305 ± 0.0067
  against −42.1368).
- Bugs reported by the cluster, fixed here:
  - `kt_jets`: NaN momenta gave ii = 0 and a write to act(0), "double
    free". Now `accept` rejects configurations with non-finite momenta
    (and flags them as degenerate for P2B), and `kt_jets`/`njets` exit on
    ii = 0.
  - `hcacc = hcacc + sg…`: shape mismatch (160 vs 150 in mode 2). Now
    `hcacc(1:nv)`.
  - Bit-identical on normal events (plain, P2B, mode 1, lo).
- Reasoning towards the cause:
  - At the fixed point the NNLO is plausible (14–18% of LO), and mode 2 in
    a window there with the τ_zQ observable reproduced mode 1 (ZFIX, 6 Oct).
  - The excess is everywhere in ZEUS, including Q² 125–250 next to that
    point.
  - Suspect: the ZEUS jet selection as observable, i.e. a non-cancellation
    that appears only at NNLO.
- Test running (`runs/zwin`): the same window (x = 0.01, Q² = 400) with
  ZFIX = 1 (τ_zQ) against ZFIX = 2 (ZEUS selection, new option), all
  parts. A window is also small enough for an NNLOJET reference.
- Replays of the cluster's crash seed (4711) and blow-up seeds (2087, 2367)
  running, with the old and the fixed binary.

**Cause found (7 Oct evening): inverse-variance iteration weighting biases
r.**
- Window test (`runs/zwin`; x = 0.01, Q² = 400, relative width 1e-3):
  - z1 (τ_zQ) against the fixed-point runs, part by part (normalised to b0):
    b1, b2, lo, kp agree; vi −0.9 ± 0.7; **r +4.7 ± 1.4 (5e-3), +15.2 ± 4.2
    (2e-4), +13.4 ± 5.3 (1e-4), +27.7 ± 8.6 (1e-5)**;
  - z2 (ZEUS selection) grows like the full ZEUS result. Statistics too
    low to separate it from z1.
- The cells are combined over iterations with weights 1/σ²_it of the VEGAS
  target. For a heavy-tailed integrand like R − D that is biased. In the
  target cell (total, 1e-5), equal-weight iteration average against the
  weighted cell:
  - window r: −3300 ± 13 against −3249 (+51);
  - fixed point (1M points): −3281 ± 3 against −3277 (+4);
  - **ZEUS P2B production r (300 seeds): −8620.4 ± 8.5 against −8570.6,
    bias +49.8 ± 7.1 pb (7σ)**;
  - b2, vi, kp, lo, b1: |bias| ≤ 0.5 pb.
- The bias makes r less negative, i.e. the NNLO too large: the direction
  and size of our excess over NNLOJET.
  - Corrected at that one cell: P2B NNLO 43.2 → −6.6 ± 10.7 pb at 1e-5,
    against NNLOJET 11.1 ± 1.2 (1.7σ).
  - The other cells need reruns.
- My earlier dismissal of the 6 Oct ZFIX "vi 1–2% low" as ≈ 2 pb was wrong
  in size (vi is large, so 1–2% is a big fraction of the NNLO). The new
  window data show vi consistent (−0.9 ± 0.7); the real problem was r.
- Fix: `VEGAS_EQUAL=2` (equal-weight iterations 2…itmx) for r; the
  TECHDIFF runs already used it.
- Reruns (`runs/zeq`):
  - ZEUS P2B r 300 seeds (same seeds as r2);
  - window r for z1 and z2, 80 seeds each;
  - fixed-point r 120.
- **Equal-weight reruns (`runs/zeq`, VEGAS_EQUAL = 2):**
  - ZEUS P2B NNLO total (b2, vi, kp existing + r 291 seeds + TECHDIFF): 10.6
    ± 1.5, 12.6 ± 1.9, 11.2 ± 2.4, 11.2 ± 2.7, 12.3 ± 3.1, 11.8 ± 4.4,
    13.7 ± 4.2, 14.6 ± 4.7, −1.7 ± 7.5, −8.5 ± 9.2 (2e-2 … 1e-5), against
    NNLOJET 11.12 ± 1.19.
    - Flat and in agreement down to 1e-4.
    - 3e-5 and 1e-5 low (1.7–2σ); r mean is ≈ 2σ below the median there
      (outliers).
    - Some high-p̄_T/m₁₂ bins are 2–4σ high at 2e-4…1e-4.
  - Window z1 against the fixed point (both equal weights): agree (one
    3σ column at 5e-3).
- Cluster instructions updated (section "UPDATE 7 Oct evening"): rerun all
  r with VEGAS_EQUAL = 2, outlier report, debug replays.

**Cluster report 8 Oct (merged; `nnlo21/cluster/REPORT-2026-10-08.md`).**
- Iteration-weighting bias confirmed in the cluster outputs: B r ≈ +40 pb;
  F r +3.9 … +9.5 pb/GeV² (2e-4 … 1e-5).
- ZEUS with equal weights:
  - plain means agree with NNLOJET down to 2e-3;
  - with flagged seeds left out, all τ_cut within 2.1σ;
  - below 1e-3 the plain means are undefined: 59 of 3500 r seeds contain
    garbage events;
  - high-p̄_T/m₁₂ excesses at large τ_cut (+6σ at 2e-3 in p̄_T 30–60).
- **Correction:** with equal weights the fixed point (F) is not flat:
  11.6, 10.1, 8.1, 6.3, 6.1 at 2e-3 … 1e-4. The "plateau ≈ 8.5" (6 Oct) and
  "10.3–11.5" (cluster) were lifted by the weighting bias.
- The edge shift (+22 pb) was mostly a weighting artefact (equal weights:
  median +7.7 at 1e-4).
- Garbage events (replays 3197, 4637): mapped dipole Borns with an
  almost-zero-energy parton (E ≈ 6e-4 GeV), dipoles 1e61 … 1e179, against
  ≈ 4e13 for the legitimate large terms. Finite and positive, so the
  non-finite guard missed them.
  - New guard: an event with |D| > 1e30 and |D| > 1e12 |R| (or a
    non-finite D) is dropped with all its dipoles, in all modes, and counted
    (`DIPGARB`).
  - Bit-identical on normal runs here; the ≈ 3e-4 of events it counts were
    already dropped as non-finite.
- Cluster instructions: rerun r (B, F) and all TECHDIFF sets with the guard
  and VEGAS_EQUAL = 2; report the tail statistics and F's flatness.

**NC/CC for NNLO 2+1, stage 1a: me31 with photon + Z (8 Oct, branch
`2026-10-nnlo21-ncc`).**
- AK: implement the NC/CC ingredients now; keep the coding style; option 1
  for the axial pair-flavour terms (as disorder's MATFOR and HOPPET).
  HOPPET's InitC3NNLO sets the singlet parts of C3 to zero ("no impact on
  Z case").
- New `dis31/ew31.f90`: helicity couplings c(hq, hl) = Q_f q_l + z_f z_l
  Q²/(Q² + MZ²) as in MCFM's qqb_z2jet (zcouple), disorder's MZ, MW; e∓,
  ν, ν̄ (ν: × √2 for the lepton spin average). Photon default.
- `me31`: m2 weighted per helicity; four-quark with helicity couplings per
  line.
  - Fixed against disorder's NC MATFOR (`harness_me31`): the DIS crossing
    puts the outgoing quark of the first line into MCFM's incoming-quark
    slot (helicity label exchanged for quarks). The pair line (5,6) is in
    MCFM's own orientation (not exchanged).
  - With photon + Z the interference of the boson on different quark lines
    is dropped (different flavours, and its pair = incoming flavour member
    in the identical-quark |D|², |E|²); the direct × exchange interference
    is kept.
  - Result: me31/MATFOR − 1 ≤ 1.5e-13 for e⁻, e⁺, ν, ν̄ (50 points, all
    flavours); photon ≤ 1e-13. (`-Zonly`, `-intonly` are not in ew31.)
- Size of the dropped terms (`dis31/tests/axsum31`, full 3+1 phase space,
  labels symmetrised, e⁻), against disorder's inclusive NC:
  - (2.4 ± 0.5)e-4 at x = 0.01, Q² = 400 (LO 104.3, NNLO coefficient −2.18);
  - 7.4e-8 at x = 0.1, Q² = 5000 (LO 3.3e-2);
  - −2.4e-9 at x = 0.3, Q² = 20000 (LO 3.9e-4).
  That is ≈ 2–6e-6 of LO and ≤ 1e-3 of the NNLO coefficient. For e⁺ the sign
  flips at low Q² (F3-like). Negligible.
- Next: born31 (production 3+1 Born, colour correlations), me41, dip41,
  virt31/iop31, the nlo31 flavour sums (from charges to couplings); then CC.
- **Stage 1b: born31, me41, the dipoles with photon + Z (8 Oct).**
  - born31: q g g / g q q̄ g with per-helicity weights; the parity-derived
    right-handed quark line weighted c(1,l)² + c(2,l')²; in the
    spin-correlated case 2 c(pq,pl)²; four quarks as me31.
  - me41: xzqqggg weighted as me31. `msq_gqqQQg`/`makemb_photon` (our
    photon version of MCFM's msq_ZqqQQg) take per-helicity couplings and a
    flag `nll`: the line-line interference within a pairing class is
    dropped (terms A, D, F, G from the boson on each line alone); the
    direct-exchange terms B, C, E are kept, as in me31. The leptons enter as
    (4,3): lepton label exchanged.
  - Tests with EW31 = "1 0", "1 1", "1 2" (e⁻, e⁺, ν):
    - harness_born31: msq = me31 4e-16, colour conservation 1e-15,
      polarisation sum 4e-14, soft and spin-correlated collinear limits as
      the photon;
    - harness_lim41: every me41 channel in its collinear limits against
      me31, 1.7e-4 at the smallest y, as the photon;
    - harness_dip41: me41/Σ dipoles − 1 = 8e-5, as the photon.
    The ν case (pure Z, maximal parity violation) checks the helicity
    conventions of every line.
- **Stage 1c: nlo31 flavour sums with photon + Z (8 Oct).**
  - The photon sums use unit charges and `solve3` (charge structures); with
    Z the couplings depend on helicity, so for `EW31 = "mode lepton"`
    (mode 1) nlo31 evaluates every channel with ew31's basis couplings
    (`ew31_basis` 1: c = 1 on the diagonal helicity pairs, 2: off-diagonal,
    on the line(s) of flavour `ew31_bf`, zero elsewhere) and weights per
    flavour with (c11² + c22², c12² + c21²)/2 of `ew31_cpl(-f)`. Exact with
    option 1 (linear in the squared couplings per line) and tree-level
    parity. Antiquarks from the conjugated evaluations (C invariance, as
    me31/me41 implement them).
  - Pitfall found on the way: the basis couplings must carry ew31's
    quark/antiquark helicity exchange. Without it the dipoles whose mapped
    Born has the line crossed (initial-state g → q q̄) got the wrong relative
    helicity; the bare 4+1 was fine. Real minus dipoles was off by up to 5 %
    (gluon and antiquark channels), opposite for e⁻ and e⁺.
  - Evaluations per point: Born 10 (photon 6), real 16 (photon 10); r costs
    1.8× the photon. vi stops with Z (no one-loop amplitudes yet).
  - Checks, e⁻, e⁺, ν (point by point, against the explicit sum over
    flavour assignments with symmetry factors):
    - `nlo31 chk`: lo and r agree to all 11 printed digits;
    - kp (scratch test of born_part against an explicit kp_conv sum): 13
      digits;
    - every 4+1 channel alone (me41 and real minus dipoles) ≤ 1e-8;
    - photon (no EW31): lo, kp, r and the mode-1 cells bit-identical to
      d03c2a5.
  - With MZ → 1e10 the NC sums differ from the photon at 1e-3 point by
    point: the photon keeps the line-line interference that option 1 drops
    (e_q e_Q for different flavours, and the identical-flavour one). It is
    odd under pair exchange and integrates to zero, so this is not a
    pointwise test.
  - `build_nlo31.sh` now compiles ew31 (since d03c2a5 it did not build).
- **Stage 1d: below the cut (sliced21, lp21, hard21) with photon + Z (8 Oct).**
  - The 2+1 Born has one quark line, so nothing is dropped there (option 1
    concerns only 3+1/4+1 and the two-loop closed quark loop, see below).
    DISENT's MATTHR splits the quark Born into S = (p1·p6)² + (p2·p7)² and
    O = (p1·p7)² + (p2·p6)² (C2 ± C3). In ew31 terms the quark Born is
    QQ/(S+O) (S w1 + O w2) with w = `ew31_w(f)`, and the gluon Born is GQ
    Σ_q (w1 + w2)/2 (the parity-odd part, odd under q ↔ q̄, dropped as in
    MATTHR).
  - hard21 returns the hard function per helicity class as well (optional
    `hk`). The one-loop (and two-loop) ratios differ between S and O, so
    with Z they need separate weights. Classes: amplitudes {1,2,5,6} = S,
    {3,4,7,8} = O.
  - lp21: the beam table is now unweighted in 5 classes (up q, up q̄, down
    q, down q̄, gluon); `lp21_born` takes the couplings wk(k, c) per helicity
    class and beam class and combines the classes with hard21's Born
    fractions. sliced21 computes the couplings (`class_weights`; photon
    e_q², NC from ew31) and splits b0 with S/(S+O). **The beam grids of
    mode 2 must be rebuilt** (`sliced21 mktab`); `lp21_grid_load` refuses
    the old 3-class files.
  - `nnlo21/tests/harness_born21` (with disorder's flags; `build_harness.sh`)
    against DISENT, 50 points, photon, e⁻, e⁺, ν:
    - hard21 classes S:O = S/O to 6e-15;
    - quark Born = MATTHR to 3e-15 for f = −5..5, gluon Born 8e-16;
    - one loop per class: H1/H0 = hard_fact + (QQ ± QQ3)/2 over the class
      Born (VIRTHR's C2 QQ + C3 QQ3, QQ3 from VIRT3PV) to 6e-13.
  - On the way: hard21 has a step of ~1e-5 relative in H1 within 1.3e-3 of
    the region boundaries s_ij = s45 (NNLOJET's displacement v → v − 1e-3
    for 1 − v < 1e-3, inherited; DISENT is smooth there). It is negligible
    after integration and was there before; the harness uses a fixed seed.
  - Photon: sliced21 b0, b1, b2 at (0.01, 400) identical to HEAD to all 9
    printed digits in all ten τ_cut cells; test_lp21 unchanged (2e-5, the
    beam-integral tolerance).
  - NC with MZ → 1e10 (scratch copy of ew31): b0, b1, b2 for e⁻ and e⁺ equal
    the photon to all 9 digits (exact test of the weights, the q/q̄ classes
    and the gluon sum).
  - With the Z at (0.01, 400), τ_cut = 1e-5: b0 +0.7 % (e⁻), +0.3 % (e⁺)
    against the photon; ν 0.33 pb/GeV² (pure Z).
  - Beam table: the q̄ classes have 4e-3 deviation relative to their own
    maximum (small sea at large ξ). The interpolation is linear in the
    classes, so q + q̄ has exactly the old error (1e-3).
  - Not done yet (stage 3): the two-loop N_F,γ term (boson on a closed
    quark loop) still uses photon charges with Z. With Z its vector part is
    Σ v_q, and its axial part is the analogue of the dropped pair-flavour
    terms. To be decided together with the 3+1 one-loop closed-loop terms.
  - Next: NLO 2+1 with Z (b1 + lo) against disorder's DISENT with
    `-includeZ` (tau2_nlo, as the photon validation of 1–3 Oct), once
    stage 2 is in, in one campaign.
- **Stage 2: one-loop 3+1 (virt31) with photon + Z (8 Oct).**
  - qq̄gg: each helicity term (hq, lh) of the BDK amplitudes weighted with
    c(hq, lh)² (ew31_cpl of the line; MCFM's label convention without
    exchange). Four quarks (`fourq_cpl`): cq(polq, 3 − polz) and
    cQ2(polb, 3 − polz) (the leptons enter as (4,3), as in me41). The
    mapping was fixed by a scan of the label exchanges against me31 (only
    these two and their parity mirrors give tree = me31 in every channel).
    Option 1 as in me31: in 2 Re(tree* loop) the cross terms of the boson on
    the two lines are dropped, direct × exchange kept. The boson on a closed
    quark loop stays out (vector: Furry; axial: only via the top–bottom
    splitting, analogous to the pair-flavour terms). `virt31_keepint`
    (diagnostic) keeps the line–line interference.
  - `harness_virt31` (EW31 env added), photon, e⁻, e⁺, ν, ν̄: tree = me31 to
    6e-14, double pole 6e-14, single pole (−⟨I⟩ with born31's colour
    correlations) 1.3e-12. `harness_iop31`: V + I poles cancel to 2e-13,
    RG and μ_F checks as before.
  - Finite part against NNLOJET's Z functions with ewcalc couplings set to
    ew31's (no width), e⁻ and ν, with ours keeping the line–line interference
    (`dis31_nnlojet/harness_v31z`, outside the repository; README there):
    with each point plus its mirror image (y → −y), and the gluon channel
    averaged over q ↔ q̄ as NNLOJET's sB2g1Z does, K equals the photon's
    (ν: twice), the single pole ratio is 1, and (v0 − K X0)/tree = −5.79471429
    (qq̄gg) and −1.71982418 (four quarks) at every point, the photon's scheme
    constants.
  - **Open:** point by point the reflection-odd (ε-tensor) parts of the
    one-loop interference differ between virt31 (MCFM's crossing) and
    NNLOJET. They are neither opposite nor absent in one of the codes. They
    exist only with Z (c_LL ≠ c_RR) and integrate to zero for observables
    symmetric under the reflection (everything we compute). A third
    reference would settle it.
  - nlo31 vi with Z enabled. The basis sums assume m2(1,1) = m2(2,2), which
    the one-loop odd part breaks, so they keep exactly the reflection-even
    part. Checked against the explicit flavour sum (scratch test with
    arbitrary flavour weights): equal to 9 digits for P + mirror, e⁻, e⁺, ν.
  - Photon: nlo31 vi bit-identical to the previous commit.
  - Still to do for full NC: the two-loop closed-loop N_F,γ term in hard21
    with Z (stage 3), CC, and the NC validation campaign (NLO 2+1 against
    DISENT `-includeZ`; NNLO against NNLOJET with Z).

## 2026-10-08 (evening) — NC stage 3, closed-loop terms, dropped-term notes

- **Correction (virt31, also photon):** virt31's comment said the photon
  couples to no closed quark loop in q q̄ g g (Furry). Wrong: with three
  gluons on the loop (two external, one to the quark line) the vector term
  survives (BDK's A6^v, colour (N − 4/N), ∝ Σ e_q; MCFM's `mqqb_vec0`).
  virt31 (ported from MCFM's colourchoice-0 path) had left it out. Now
  included (`fvs.f`, `fvf.f`, `a64v` from MCFM, unchanged), photon and Z
  (loop coupling Σ_f v_f, `ew31_cv`). Checked against MCFM's own
  `xzqqgg_v` at DIS points: our closed-loop term / MCFM's c·c_V-weighted
  `mqqb_vec0` equals our main term / MCFM's `mqqb`, at every point, all
  three channels, photon, e⁺, ν. Finite (the pole checks of harness_virt31
  and harness_iop31 unchanged, ≤ 4e-13). Size for the photon: −0.03 % of
  vi in NLO 3+1, 3e-6 of vi at (0.01, 400): no visible effect on any
  earlier result. NNLOJET's DIS process lacks it too (the photon
  finite-part agreement of 2 Oct was therefore without this term on both
  sides).
- nlo31: the closed-loop term is linear in the line's coupling (c·c_V), so
  it is weighted apart (`virt31_vloop`): photon −3 e_q times the d value;
  NC basis values with unit loop coupling, weights `ew31_wl`. Checked
  against the explicit flavour sum: photon point by point (11 digits), NC
  point + mirror (11 digits). A stale `virt31_vloop` from a previous q q̄ g g
  call leaked into four-quark channels during development (found by this
  check, fixed before commit: reset in `virt31_ren`).
- **Stage 3, two-loop N_F,V term (hard21's G):** linear in the coupling
  ratio nfz per helicity amplitude; lp21 evaluates hard21 with the term off
  and with ratio 1 and combines per class with nfk(k, c) from sliced21
  (`class_weights`): Σ_class c(h,l) c_V(l) / Σ c(h,l)², c_V the vector
  coupling of the loop (photon: Σ e/e_q and 1/11, as before). The axial
  part is dropped. Photon b1, b2 at (0.01, 400) identical to before (all
  digits); NC with MZ → 1e10 equal to the photon (e⁻, e⁺).
- **Dropped terms:** `docs/nc-dropped-terms.md` (for the paper): every
  term, whether it exists, kept/dropped, size. Open: the axial
  closed-loop terms (rows 5, 10, 12) need a numerical estimate; HOPPET's
  singlet treatment for F2/FL with Z to be documented.
- NC validation started (x = 0.1, Q² = 5000, e⁻ and ν): NLO 2+1 (sliced21
  b1 + nlo31 lo) against DISENT with `-includeZ` (tau2_nlo, new flags),
  136 jobs via the dispatcher (`runs/ncz1`).
- **NC validation of NLO 2+1 against DISENT with Z (8 Oct, evening)**,
  x = 0.1, Q² = 5000 (y = 0.49), τ_zQ bins, `runs/ncz1`: sliced21 b1 (8
  seeds) + nlo31 lo (20 seeds, psmc) against tau2_nlo's DISENT reference
  (40 × 3M events), (b1 + lo)/DISENT − 1 in %:
  - e⁻, total bin [0.05, 0.5): −0.13 ± 0.31 (τ_cut 2e-4), +0.27 ± 0.35
    (1e-4); all six bins within ~1σ at 1–2e-4 except [0.2, 0.3) −1.7 ± 1.1
    and [0.4, 0.5) +7.6 ± 3.0;
  - ν (pure Z), total bin: −0.32 ± 0.29 (2e-4), +0.09 ± 0.33 (1e-4); the
    same two bins −2.4 ± 1.1 and +5.2 ± 2.3. The two leptons share the lo
    seeds, so these are one correlated fluctuation (or an effect of these
    bins at x = 0.1), not two.
  - Below 3e-5 the sums drift up by ~1 % (2σ): the lo runs used the old
    technical cut 1e-9 W², known to bias τ_cut ≤ 3e-5 (5–6 Oct).
  - **Correction:** the first ν reference (s097–s136) was wrong: my new
    `-neutrino` in tau2_nlo did not set `Zonly` (set_parameters does: no
    photon for neutrinos), so the neutrino coupled to the photon; ours/DISENT
    came out flat at −60 % in every bin. Fixed (tau2_nlo), rerun s137–s176.

## 2026-10-08 (late evening) — cluster DIPGARB reruns merged; axial estimate; HOPPET with Z

- Cluster (`origin/2026-10-nnlo21-cluster`, REPORT-2026-10-08b, merged into
  2026-10-nnlo21 by fast-forward and into -ncc): the garbage-dipole guard
  removes the 1e60–1e179 values but the 1e5–1e6 tail of r stays (52 of
  3500 B seeds flagged). Replays: real events deep in collinear/soft
  corners (s_min/W² ~ 1e-9, at the technical cut) whose dipole is orders
  of magnitude larger and has its mapped Born above τ_cut — not garbage.
  ZEUS NNLO plain means agree with NNLOJET to 2e-3 (11.8, 12.0, 12.4, 14.1
  ± 1.1–1.3 pb against 11.12 ± 1.19), undefined below 1e-3; DIAG within
  ≤ 2σ down to 1e-5. Fixed point still not flat (10.1 → 6.2 from 1e-3 to
  1e-4). Advice to AK: no more cluster jobs until the tail has a method fix.
- Axial closed loop of the one-loop 3+1 (row 5 of nc-dropped-terms):
  0.0041 ± 0.0002 pb = 2.7e-4 of vi (qq̄gg channels) in the HERA 3-jet
  set-up, e⁻, MCFM's BDK axial amplitudes (heavy-top approximation);
  m_t-independent (173 vs 500 GeV).
- HOPPET with Z (section 2 of nc-dropped-terms): option-1 flavour weights at
  all orders; C3 singlet zero at NNLO (consistent); at N3LO the fl11 terms
  of F2/FL carry MVV's photon normalisation 3⟨e⟩ but get the Z weights
  v² + a² — not the closed-loop coupling Σ v_f our 2+1 side uses. To be
  resolved before an N3LO 1+1 with Z by P2B.

## 2026-10-09 — the heavy tail of r: non-massless psmc momenta (fixed)

- **Correction of the 8 Oct conclusion** ("real events deep in collinear/soft
  corners … not garbage", cluster REPORT-2026-10-08b, and my reading of it):
  all 81 large events of the cluster replays (seeds 2042, 2049) have the real
  cut away (τ₂ ≈ 4e-9 … 3e-8) and an accepted dipole whose mapped Born has a
  **negative-energy parton** (e.g. E = 183, 173.7, −31.8 GeV; dipole 2e17).
  Real and dipole do not cancel at all there.
- Cause: psmc's emission (`emit`) built momenta with m²/E² up to 1e-6 (a
  few per mille of 4-parton events with m²/E² > 1e-12; scratch harness, x =
  0.01, Q² = 400, edge 1e-12). Two sources:
  - `tbasis` projected a unit vector off the two light-like momenta, losing
    ~1/θ² in the orthogonality when the emitter–spectator pair is itself
    nearly collinear (second emission);
  - even with an exact basis, one transverse direction then has components
    ~1/θ, so the new momenta are differences of large numbers.
  The tiny invariants then come out negative (the technical cut uses |p_i·p_j|,
  so such points pass), y > 1 in the FF mapping, negative energies in the
  mapped Born. With P2B such an event gets the plain-slicing weights (the
  6 Oct fallback), which accept the garbage dipole above τ_cut.
- Fix (`dis31/psmc.f90`): e1 from ε(a, b, axis) like e2, and `emit`/`tbasis`/
  `eps` in quad precision (`qp`), the result rounded. Largest m²/E² of
  4-parton events now 1e-9 (one in 1e6), negative-energy mapped Borns 485 → 143
  per 583k events; the rest are genuine near-zero-energy Borns (both partons
  collinear to the incoming one). No measurable cost in run time.
  The new basis is rotated with respect to the old one, so points differ; the
  phase-space volume agrees (300k points: 2.437 ± 0.031 vs 2.453 ± 0.030 for
  3 partons, 69.2 ± 1.3 vs 71.6 ± 1.3 for 4).
- Counting test (instrumented scratch copy of nlo31, 3 seeds × 3M points
  each, thserv05): dipoles with a negative-energy Born that are accepted while
  the real is cut, ZEUS P2B (mode 2): **453, 568, 809 (up to 8e30) → 0, 0, 0**;
  fixed point (mode 1): 0 → 0. All negative-energy Borns: 3.3–5.2k → 0.5–0.8k
  (mode 2), 1.9–2.6k → 0.35–0.48k (mode 1).
- Fixed point r, 48 seeds old vs new psmc (thservs/desktops): agree at every
  τ_cut (new − old = −0.9 ± 1.2, …, −9.7 ± 10.3 (1e-4), −15.6 ± 15.2 (3e-5),
  3.1 ± 18.6 (1e-5)), no outliers in either. So the fix does **not** explain
  the fixed point's non-flatness (it had no accepted garbage dipoles); that
  stays open.
- An exact local replay of the cluster blow-ups is not possible: the d638917
  build matches the cluster's iteration 1 bit for bit but drifts from
  iteration 2 (−8520.976 vs −8520.246), so VEGAS adapts differently.

## 2026-10-09 (afternoon) — HOPPET's N3LO fl11 terms with Z

- **Photon: exact.** MVV's fl11 decomposition (non-singlet e_q² · 3⟨e⟩,
  pure singlet ⟨e⟩²/⟨e²⟩ − 3⟨e⟩, gluon ⟨e⟩²/⟨e²⟩) reproduces the diagram
  coupling e_q Σe per flavour because 3e_q² − e_q = 2/3 for both charges
  (exact rational check nf = 3, 4, 5). Numerically (program
  `~/work/disorder-comparisons/hoppet-fl11/fl11z.f90`, outside the repo):
  HOPPET's photon fl11 = exact couplings to 1e-9. Closes the open photon
  question of 8 Oct.
- **Z, γZ: approximate.** HOPPET weights the fl11 sum with v_q² + a_q² and
  2 e_q v_q; the diagram couplings are g_V,q Σ g_V (Z), e_q Σ g_V + g_V,q Σ e
  (γZ), and (Σ g_V)², 2 Σe Σg_V for the gluon. For γZ even the sign of the
  fl11 term differs (x = 0.01, Q² = 400: −0.37 vs +1.66, rest of the N3LO
  coefficient −159). Size: ≲ 3·10⁻⁵ of F2^γZ, ≲ 10⁻⁵ of F2^Z (FL smaller),
  i.e. 1–2% of the N3LO coefficient of F2^γZ, ≲ 0.3% for F2^Z; below 1e-6
  in the e±p cross section at HERA. Table in `docs/nc-dropped-terms.md`
  section 2. For a consistent N3LO 1+1 with Z by P2B the fix belongs in
  HOPPET (weight the fl11 pieces with the diagram couplings).
- F3: HOPPET's C3 at N3LO has MVV's fl02 term through the valence
  function (`C%NS_V = cfN3LO_F3NS_val`). **Correction (same afternoon):**
  this entry first said C3 had no fl02 term at all (a grep missed the NS_V
  line) and asked whether W should use it; that was wrong, the 8 Oct text
  was right, the question is void.
- HOPPET bug found on the way: `InitCoefHolder` tests the optional
  argument `use_exact_cf` directly (`coefficient_functions_holder.f90`
  line 123) instead of its local copy `use_exact`, so a call without it
  segfaults (`structure_functions` always passes it, so HOPPET's own use is
  fine).
- Patch for HOPPET (for AK to port; outside the repo,
  `~/work/disorder-comparisons/hoppet-fl11/hoppet-fl11-Z-couplings.patch`,
  against d891044 = 2.3.0-272): the fl11 terms of the Z and γZ structure
  functions get the diagram couplings (new `fl11_split`, internal `fl11_Z`,
  `fl11_gZ`); photon, F3 and W unchanged. Tests, pristine vs patched build
  of the same commit (`sfcmp.f90`, `sfvfn.f90`): EM, F3^Z, F2^W+ at N3LO
  bit-identical; ΔF2^γZ, ΔF2^Z, ΔFL equal (exact − HOPPET) × (α_s/2π)³ of
  `fl11z.f90` to 4 digits (e.g. 2.952e-5 vs 2.0284 × 1.455e-5 at x = 0.01,
  Q² = 400); variable flavour number finite, n_f = 3 (Σe = 0) left at zero.
  Second patch `hoppet-initcoefholder-optional.patch` for the optional-
  argument bug.
- fl02 (F3, both bosons on the closed loop): HOPPET's couplings for Z and
  γZ are the loop sums, exactly (`fl02z.f90`, equal to 7 digits); no patch
  needed.
- HOPPET ctest: pristine 18/18; patched 17/18, `structure_functions_example`
  differs exactly in F1/F2 Z and γZ at N3LO (N3LO coefficient of F2^Z up to
  3%, F2^γZ up to 22% where it is small, F1 ≤ 1%), nothing else. With the
  regenerated reference (`hoppet-fl11-Z-reference-output.patch`) 18/18.
  `structure_functions_example_flavour` cannot see it (flavour-decomposed
  SFs stop at NLO).
- History/APFEL++: HOPPET's NC weighting of fl11 is from 1c110a8 (2016) and
  unchanged; the Feb 2023 fixes from the comparison with V. Bertone (1058001,
  8c4ac2a) changed the coefficient functions (δ(1−x) pieces, fl11 parts of
  C2S3C, NC/CC difference pieces), not the couplings. The comparison plots
  (now `benchmarking/apfelpp-checks`) are EM only. APFEL++ has the same
  per-flavour effective-charge weighting of fl11 for Z and γZ
  (`ElectroWeakCharges`, `DISNCBasis`), so it shares the issue. HOPPET's
  GitHub issues (7) contain nothing on this.
- yadism (0.13.11 / 8bbfdc8): fl11 coupling routine has the right structure
  (⟨Q_b⟩Q_b′) but a mode bug: γZ orderings give ⟨Z⟩Z_q + ⟨γ⟩γ_q instead of
  ⟨γ⟩Z_q + ⟨Z⟩γ_q (d quark: 0.0734 vs 0.0198), and ZZ includes ⟨a⟩a_q.
  Confirmed with yadism's own `partonic_coupling_fl11` in a scratch venv; fix
  tested (`~/work/disorder-comparisons/hoppet-fl11/yfix.py`). Size ≤ 6·10⁻⁵ of
  F2^γZ. fl02 in F3: correct in yadism and APFEL++.
- Bug-report pages for AK to pass on: APFEL++ (V. Bertone)
  https://claude.ai/artifact/76tsS1d8vZt4vWiwNpdmkS, yadism
  https://claude.ai/artifact/UaaDjNhibLSKLns1EXjHkC.

## 2026-10-09 (evening) — fixed-point non-flatness; 5-class beam grid; event-shape set-up

AK (Fri afternoon): free rein for the weekend up to (not including) the
clean-up: (1) fixed-point non-flatness, (2) 5-class beam tables, then CC,
NNLO-with-Z validation prep, the rest.

**Beam grid (5 classes, ZEUS/inclusive Q range).** `sliced21 mktab` for the
29 nodes with the NC-branch build (46b1e4b), thA371a:
`~/work/disorder-comparisons/nnlo21/beamgrid5/zeus5_*.tab`. `tabchkg`: max
deviation / max coefficient 1.1e-3, 5.3e-3, 4.0e-4, 6.7e-3, 3.3e-3 (up q,
q̄, down q, q̄, gluon).

**Fixed point (x = 0.01, Q² = 400): where the drift is.** From the cluster's
DIAG combination (F_g_diag): the NNLO falls 10.1 → 6.2 between τ_cut = 1e-3
and 1e-4 only in the two lowest τ_zQ bins (0.05–0.1: 3.31 → 0.43; 0.1–0.2:
3.24 → 2.22); the bins 0.2–0.5 are flat within errors. The same two bins
carry the NLO power corrections (0.05–0.1: 2.75, 4.03, 4.16 at 1e-3, 1e-4,
1e-5); NLO total 14.99 (1e-3) vs 16.68 (1e-4), 16.82 (3e-5), DISENT 16.885.
- P2B for mode 1 (new: nlo31 accept with `P2BSLICE=1` in mode 1, τ_zQ bins
  of the projected 2+1 Born; without P2BSLICE bit-identical, lo and r):
  - NLO, P2B − plain (`P2BEXTRA=1`, 12 seeds lo): 0.20 ± 0.10 at 1e-3,
    0.09 ± 0.06 at 1e-4. P2B does **not** remove the fixed-point drift
    (NLO P2B 15.19 at 1e-3): the drift is the ordinary leading-power
    slicing correction (O(τ_cut log^k)), not a recoil/fiducial effect.
  - NNLO, P2B − plain (vi 24, kp 16, r 96 seeds): −3.1 ± 1.6 (1e-3),
    −1.4 ± 0.6 (1e-4), −0.2 ± 0.6 (3e-5), −0.15 ± 0.58 (1e-5). NNLO with P2B:
    6.9, 4.7 ± 0.7 (1e-4), 8.4 ± 3.6, 3.7 ± 0.7 (1e-5); plain 6.2 ± 0.8
    (1e-4), 3.8 ± 1.2 (1e-5).
- Conclusion so far: the 1e-3 → 1e-4 drift is power corrections in the
  bins next to the 2-jet edge (consistent with the NLO); below 1e-4 the NNLO
  is ≈ 4–6 within errors, with no independent reference. The "plateau 8.5"
  of 6 Oct and the cluster's 10.3–11.5 are both superseded.
- A narrow-window NNLOJET reference is not feasible (NNLOJET samples the
  whole lepton phase space; Q²/y are only selectors). Instead: the same
  observable as an **event shape in the inclusive cuts** (AK's plan C):
  - our side: new `ZFIX=3` (mode 2, τ_zQ bins, 125 < Q² < 20000,
    0.2 < y < 0.6, no window); LO 969.46 ± 0.93 pb in 0.05 < τ_zQ < 0.5;
  - NNLOJET: `dis_thrust` is NOT our τ_zQ (normalised to Σ|p| in the current
    hemisphere, no energy cut: infrared-unsafe for 2+1 Borns with an empty
    current hemisphere). Added locally `dis_tauzq` = 1 − 2 Σ_{p_z>0} p_z/Q
    (`driver/core/EvalDIS.f90`, `Observables.f90`, originals kept as
    `.orig`; separate build `build-tauzq`); runcard
    `~/work/disorder-comparisons/nnlo21/nnlojet_tauzq/template.run`.

**W exchange (CC) for the NNLO 2+1 pieces (9 Oct, evening; local 80e28ab).**
- ew31: `ew31_mode = 2`: c(L,L) = Q²/(Q² + M_W²)/(2 sin²θ_W) (labels and
  lepton cases as for Z), unit CKM in (u,d), (c,s), no coupling for b (as
  disorder's CC, HOPPET and NNLOJET's default identity CKM);
  `ew31_out(f)`: the outgoing flavour of the boson line.
- me31/me41/born31/virt31: the boson line identified with ew31_out; four
  quarks with couplings per assignment of the outgoing quarks (W on the
  incoming line or on the pair), the two assignments interfering
  (four_quark_cc; me41 via msq_gqqQQg with the exchanged pairing's couplings,
  a backwards-compatible generalisation; gluon-initiated g → q q̄ Q Q̄ with
  the W on one pair).
- **Checks (all four leptons):**
  - me31 = DISENT's CC MATFOR to 1e-12 for every incoming parton (with
    `ew31_ccmatfor`, see below);
  - me41 against me31 in all single collinear limits (harness_lim41),
    born31 (= me31, colour conservation, polarisation sums, soft and spin
    limits), the sum of dipoles against me41 in all limits (harness_dip41):
    as for photon/Z (≤ 1e-4 at the smallest parameter);
  - virt31: tree = me31, poles = −⟨I⟩ (≤ 5e-13), incl. EXX, EXY and the
    W-on-the-pair identical-quark channel; the finite part has no
    independent W check yet (NNLOJET's DISW four-quark functions: not done);
  - nlo31's sums from representative channels (generation 1, multiplicities
    per incoming flavour) = brute-force sums over all final-state flavours
    (harness_ccsum, 3+1 and 4+1, 5e-16);
  - sliced21's 2+1 Born and one loop per helicity class = MATTHR/VIRTHR CC
    (harness_born21 -noNC -CC, 1e-15 / 6e-13);
  - photon and photon + Z bit-identical to the build before (lo, r).
- **DISENT's CC MATFOR omits a physical term:** for incoming flavours the W
  cannot hit directly (d, s, ū, c̄ for W⁻) the identical-quark
  interference of the two W-on-the-pair assignments (e.g. d → d d ū). Up
  to 10% of those (small) channels pointwise. NC MATFOR keeps the analogous
  term (D·E*). Our code keeps it (default); `ew31_ccmatfor` drops it for
  comparisons. Affects disorder's CC P2B at NNLO (differential only). To
  report to AK.
- sliced21/lp21: W classes (up q, down q̄ for W⁻; b out of the down
  classes via `lp21_nodn`, separate beam tables with a flag word; gluon class
  with the two W pairs; no closed loop). W beam grid
  `beamgrid5w/zeus5w_*.tab` (tabchkg ≤ 7.7e-3).
- LO in the ZEUS set-up against NNLOJET epNJJ (W⁻, identity CKM): total
  3.7554 ± 0.0011 vs 3.7508 ± 0.0028 pb (1.5σ); per Q² bin within 0.1%
  except 500–1000 GeV² +0.38% (2.5σ, 4 NNLOJET seeds). The photon shows the
  same pattern (one bin +0.42%, 2.5σ; total 103.286 ± 0.016 vs 103.281):
  NNLOJET's single-run bin errors are somewhat optimistic.
- LO with photon + Z (ZEUS) against NNLOJET V_NC = Z+GAMMA: Z/γ per Q² bin
  within NNLOJET's errors (−0.24 … +0.17%, up to 2σ; Z/γ = 1.36 in the top
  bin).
- Event shape τ_zQ (ZFIX = 3) against NNLOJET (local dis_tauzq): LO 969.46 ±
  0.93 vs 968.87 ± 0.55 pb; our NLO (b1 + lo, 16 seeds): 249–252 pb at
  τ_cut ≤ 1e-4, power corrections of −30 pb at 1e-3 (lowest bin), P2B flatter
  from 5e-4. NNLOJET NLO (V, R) and our W NLO against DISENT running.
- **W: NLO 2+1 against DISENT -CC (9 Oct, evening)**, x = 0.1, Q² = 5000,
  τ_zQ bins, `runs/ccw1` (as ncz1: sliced21 b1 8 seeds + nlo31 lo 20 seeds
  against tau2_nlo -noNC -CC, 40 × 3M events), (b1 + lo)/DISENT − 1 in %:
  - e⁻ (W⁻), total bin [0.05, 0.5): −0.25 ± 0.31 (τ_cut 2e-4), +0.26 ± 0.35
    (1e-4); ν (W⁺): −0.46 ± 0.27, −0.06 ± 0.31.
  - Per bin the same pattern as the Z check (ncz1): [0.2, 0.3) −2.1 ± 1.1
    (e⁻), [0.4, 0.5) +8.0 ± 2.4 (e⁻, rising with smaller τ_cut; ν −0.4 ±
    2.8). The Z runs show +7.6 ± 3.0 in the same bin with the same seeds,
    so this is the correlated fluctuation noted there (or a property of the
    bin at x = 0.1), not something W-specific. The DISENT reference does
    not include the B–B identical-quark interference that MATFOR CC drops,
    but it enters only at 4+1 (NNLO), not here.
- **W one loop against NNLOJET, pointwise (9 Oct, evening)**:
  `harness_v31w` (outside the repository, dis31_nnlojet): virt31 with W⁻
  (e⁻) against NNLOJET DISWm's four-quark RV functions (C-type per
  assignment, D0g1Wa/b-type interferences), ten channel types including
  the B–B identical-quark one (d → d d ū). Point + mirror: one K for all
  channels, pole ratio 1, finite constant −1.71982418 (the four-quark
  scheme constant of the photon check) everywhere. The W finite part is
  now independently checked.
- **Correction (reflection-odd one-loop parts, open since 8 Oct):** I wrote
  "neither opposite nor absent". That was wrong. NNLOJET's one-loop
  functions equal ours with all helicities flipped. For W the odd part is
  exactly −ours in all channels. For Z it is exactly −ours unswapped and
  exactly ours with the quark and lepton helicities both swapped (scratch
  ew31), for q q̄ g g and four quarks, e⁻ and ν. The 8 Oct q q̄ g g
  comparison is now contaminated by our closed-loop vector term, which
  NNLOJET's DIS lacks; it is removed for this check. A global flip only
  changes the sign of the ε-tensor (absorptive) term. NNLOJET's
  `D0g1Wab.f` comments: "Fill only the LLL amplitudes. These correspond
  to our RRR amplitudes because ? But it works...". Which sign is physical
  is not settled; it cancels for reflection-symmetric observables (all of
  ours).
- **Event shape τ_zQ at NLO against NNLOJET (9 Oct, night)**, photon,
  inclusive cuts (ZFIX = 3), patched NNLOJET (`dis_tauzq`), V + R from 8 +
  24 seeds (`runs/esnnlojet`) against our b1 + lo (16 seeds, `runs/esnlo`),
  pb, τ_cut 1e-4 (plain): [0.05, 0.1) 49.3 ± 1.4 vs 49.3 ± 1.0; [0.1, 0.2)
  85.5 ± 0.9 vs 84.9 ± 0.2; [0.2, 0.3) 52.3 ± 0.7 vs 52.7 ± 0.1; [0.3, 0.4)
  36.1 ± 0.6 vs 36.3 ± 0.1; [0.4, 0.5) 25.8 ± 0.6 vs 26.4 ± 0.1; total 248.9
  ± 2.0 vs 249.7 ± 1.1. Same at 3e-5. P2B agrees as well (from 5e-4).
- **W NLO, the [0.2, 0.3) and [0.4, 0.5) deviations:** 20 new lo seeds
  (521–540, `runs/ccw1x`) give the same deviations (+5.6 ± 2.4 % in [0.4,
  0.5) at 2e-4). The error is dominated by b1 (sliced21, 8 seeds, ≈2 % of
  the bin against 0.5 % from lo), and Z and W used the **same b1 seeds**
  (601–608). So the shared piece is b1, not lo as I wrote above. 16 new b1
  seeds (609–624) are running.
- **Axial closed loop of the 2+1 one loop (row 10 of nc-dropped-terms, 10
  Oct, night):** MCFM's Z + jet anomaly term (`A53`, large-m_t expansion)
  crossed to DIS e⁻ p and PDF-weighted over the 2+1 Born phase space,
  τ_zQ ∈ [0.05, 0.5): axial/Born = (α_s/2π)·c with c = −5·10⁻⁵ (0.01, 400),
  −3·10⁻⁴ (0.05, 1000), −4·10⁻³ (0.1, 5000), −7·10⁻³ (0.2, 10⁴); at most
  ~1.4·10⁻⁴ of the 2+1 Born, ~10⁻³ of the NLO correction. Scratch program
  outside the repository (`disorder-comparisons/axial21`). The tree/virtual
  normalisation relies on the leading-colour double pole of A51; it was not
  checked numerically (an order-of-magnitude estimate).
- **Resolved (10 Oct, 00:30):** the 16 new b1 seeds (609–624) alone give
  [0.4, 0.5) −1.1 ± 1.9 % and [0.2, 0.3) −1.0 ± 1.0 % at τ_cut 2e-4. With
  all 24 b1 and 40 lo seeds: every bin within 1.6σ of DISENT, total bin
  −0.40 ± 0.24 % (2e-4) and +0.04 ± 0.27 % (1e-4). The old deviations were
  an upward fluctuation of b1 seeds 601–608 in [0.4, 0.5). The Z check
  (ncz1) used the same seeds, so its +7.6 ± 3.0 % there is the same
  fluctuation. Lesson: use distinct seeds per validation run.

## 2026-10-10 — which sign of the reflection-odd one-loop term is right

AK: dig into the literature; report to NNLOJET only if their code is wrong.
**Result: NNLOJET is right; our virt31 has the wrong sign of the
reflection-odd (T-odd) part. This corrects 9–10 Oct, where I suspected
NNLOJET (its D0g1Wab.f comment) and put it into the bug report (now removed).**
- Literature: Körner, Melić, Merebashvili, hep-ph/0002302 (PRD 62 096011),
  eqs. (3.7), (3.12)–(3.14): the T-odd DIS structure functions H5, H8, H9
  (quark initiated, leading and subleading colour) and the sin φ, sin 2φ
  coefficients, ε_0123 = −1, ε(k, p, p′, q) ∝ +sin φ; they agree with the
  independent cut calculation of Hagiwara, Hikasa, Kai.
- 5-point (2+1): MCFM's BDK amplitudes (disorder's port in
  src/disent_virt3.f, labels LL) give T-odd = +0.0013 × KMM (with
  2/y²(D sin φ + E sin 2φ)) at every point, the same constant for both colour
  structures (`todd5`). NNLOJET's B1g1W/Bt1g1W agree with MCFM, sign
  included (`todd5nj`). So both 5-point implementations are physical.
- 6-point: collinear factorisation (final-state pair, CS FF map, y → 0) must
  give T6/tree6 → T5/tree5 (the one-loop splitting's absorptive part
  multiplies a T-even tree). virt31: u → d g g (d ∥ g) and u → d c c̄ (c ∥ c̄,
  averaged over the azimuth) both tend to −1.00000 (`todd6`). NNLOJET's 6-point
  = −ours exactly (harness_v31w/z) → +1: consistent.
- Cause: ew31's helicity labels ("z_f exchanged for quarks", 8 Oct) were
  fixed against MATFOR, a T-even comparison that fixes the quark helicity
  relative to the lepton but not the global flip; exchanging the lepton
  labels instead (equivalently both, SWAPHEL=ql) gives NNLOJET's sign. No
  effect on any result (reflection-symmetric observables; nlo31's flavour
  sums keep only the even part). Fixing ew31 is left for the clean-up (it
  would change the trees at the rounding level during the cluster
  validations).
- Harnesses: `disorder-comparisons/dis31_nnlojet/todd` (with the KMM pdf).

## 2026-10-10 — fixed-point non-flatness: per-bin power corrections

AK: look into the fixed-point (F) non-flatness. Analysis of the cluster's
F_g_diag per τ_zQ bin (`~/work/disorder-comparisons/nnlo21/fscale.py`;
NLO reference: DISENT per bin, `runs/tcut1/tcut1.json`).
- **NLO:** the relative power correction is set mainly by τ_cut/τ_lo (τ_lo =
  lower bin edge; the best of lower edge, centre, upper edge), but the lowest
  bin is 2–4× worse at equal τ_cut/τ_lo: near the 2-jet edge the corrections
  grow faster than 1/τ_zQ. All bins within ≈ 1.5% of their plateau for
  τ_cut ≲ 2e-3·τ_lo. Side finding: on the plateau the upper three bins sit
  −0.5 … −1.1% below the DISENT reference (1–3σ; total 16.68 ± 0.05 at 1e-4
  against 16.885 ± 0.04); to keep an eye on in the cluster's NLO.
- **NNLO:** every bin is flat for τ_cut ≲ 5e-3·τ_lo (and τ_cut ≥ 1e-4, above
  the technical-cut region). Value at the largest τ_cut in that window per
  bin, pulls against the smaller τ_cut with the nested-sample error
  √|σ_a² − σ_b²| (the τ_cut points share events):
  - r_max = 5e-3: 0.56 ± 0.48 (2e-4), 2.46 ± 0.23 (5e-4), 1.58 ± 0.12,
    1.21 ± 0.08 (1e-3), 0.80 ± 0.05 (2e-3); all pulls ≤ 2.2; total 6.61 ± 0.55;
  - r_max = 2.5e-3: total 6.01 ± 0.71, pulls ≤ 1.2;
  - common τ_cut: 6.3 ± 0.6 (2e-4), 6.2 ± 0.8 (1e-4).
  With r_max = 1e-2 and 2e-2 the lowest two bins fail (pulls 3–8), i.e.
  the window ends at about 5e-3·τ_lo.
- **Conclusion (replaces "fixed point not flat"):** each bin has a plateau
  once τ_cut ≲ 5e-3 times its distance to the 2-jet edge; the total drifts
  between 1e-3 and 2e-4 only because the lowest bins reach their plateaus
  last. NNLO at the fixed point ≈ 6.3 ± 0.6 pb/GeV² (all bins, 0.05–0.5).
  Below 1e-4 the upper bins dip by 1–2σ where the technical-cut correction
  rcorr grows (≈ 6 at 1e-5): the known technical-cut residual, not slicing.
  For a τ_zQ-binned observable a per-bin choice of τ_cut is legitimate (each
  bin is its own cross section) and gives smaller errors. A Born-dependent
  cut τ_cut·τ_zQ(Born)/τ_ref would be the same thing for binned τ_zQ; not
  implemented.
- **ZEUS (B set) does not follow this rule:** the bins converging last
  against NNLOJET are the hard ones (p̄_T 30–60: +6σ at 1e-3, m₁₂ 65–120:
  +4.8σ), not those nearest the dijet edge, plus the known fiducial √τ in the
  high-Q² NLO bins. From 1e-4 all bins are within ≈ 2σ (p̄_T 22–30 +3.2 at
  1e-4, +1.2 at 3e-5). Why the hard bins converge slowly is open.
- Still no independent reference for the fixed-point value: the E set on
  the cluster (τ_zQ bins in inclusive cuts against NNLOJET) is that test.

## 2026-10-11 (Sunday) — the three open questions

AK: look into (a) the hard ZEUS bins, (b) the fixed-point NLO offset, (c) the
technical-cut residual; not critical.

**(b) Fixed-point NLO offset: none (corrects 10 Oct).** Per bin with the
DISENT reference error: [0.2, 0.3) −0.036 ± 0.028 (1.3σ), [0.3, 0.4) −0.015 ±
0.024, [0.4, 0.5) −0.008 ± 0.017 at τ_cut 1e-4; the total at 1e-4 still holds
the lowest bins' power corrections (−0.07 each); at 3e-5 16.82 ± 0.08 against
16.885 ± 0.045 (0.7σ). My "1–3σ" compared that total. The ZEUS NLO is fine
too (−1.1 ± 0.8% at 1e-4 against NNLOJET; all bins within errors).

**(c) Technical-cut residual: understood.** Decade-by-decade TECHDIFF runs
(`runs/techd`, current build with quad psmc, psmc edge 1e-14, VEGAS_EQUAL=2):
- fixed point (mode 1, 40 seeds each), total at τ_cut 1e-5: d1 [1e-10, 1e-9] W²
  11.0 ± 5.5 (median 5.0; = the cluster's rcorr 5.96 ± 0.47), d2 [1e-11,
  1e-10] 1.4 ± 0.8 (median 0.5), d3 [1e-12, 1e-11] 0.19 ± 0.27 (median 0.002);
  all ≈ 0 for τ_cut ≥ 1e-4;
- ZEUS (P2B, 60 seeds), d3 [min(1e-12 W², 1e-10 Q²), min(1e-11 W², 1e-9 Q²)]:
  −0.24 ± 0.35 (1e-4), 0.1 ± 0.3 (3e-5), 8.3 ± 5.5 (1e-5, median 0.2: a few
  heavy-tailed seeds), against the cluster's d1 (rcorr) 2.2 / 16.6 / 46 and
  d2 (conv) −0.1 / 2.3 / 9.0.
- Each decade grows like 1/τ_cut (real and dipole on opposite sides of the
  τ_cut boundary in a zone of width ~ the mapping's τ₂ shift), and the series
  converges by ≈ 5 (ZEUS) to 10 (fixed point) per decade. Below 1e-12 W² the
  decades give τ_cut-independent offsets and the garbage-dipole count rises
  (double precision in the dipoles; quick scan, 50k points): do not go
  lower.
- Consequence: for τ_cut ≥ 1e-4 nothing to do. Below, d2 must be added as
  well (the cluster added only d1): ZEUS +2.3 ± 1.3 pb (3e-5), +9 ± 2 pb
  (1e-5); fixed point upper bins +0.5, +0.15, +0.06 at 1e-5 (part of the
  1e-5 dips of 10 Oct).

**(a) Hard ZEUS bins: not slow convergence but an upward drift in ln τ_cut,
present in most bins.** Bin-4 data, NNLO change per unit ln τ_cut in % of
the bin's LO (nested-sample errors), windows 5e-3→5e-4, 2e-3→2e-4,
1e-3→1e-4: total +1.3 ± 0.5, +0.6 ± 0.6, +1.0 ± 0.8; p̄_T 30–60 +3.3 ± 0.9,
+4.6 ± 1.4, +4.0 ± 1.2; p̄_T 22–30 ≈ +2.4; m₁₂ 65–120 ≈ +2.4; soft bins
+1 … +2 (≤ 2σ). With d2 added the total reads ≈ 11–12 (2e-2…1e-3), 14.8 (1e-4),
≈ 17.6 ± 3.7 (3e-5), ≈ 17 ± 6.5 (1e-5) against NNLOJET 11.12 ± 1.19. The fixed
point (τ_zQ, plain slicing, upper bins) shows no such slope (0 ± 2.5% of NLO).
The NLO in every ZEUS bin is flat (< 0.1% of LO per unit ln).
- Excluded: the beam-function grid (b2's τ_cut slope from the grid equals the
  direct evaluation to median 1e-6, max 6e-4, at hard and soft Born points;
  0.3% would be needed; `disorder-comparisons/nnlo21` scratch harness hb2).
- Not usable: the 7 Oct window runs (zwin z2; pre-psmc-fix, weighting bias).
- Scale of the problem: b2 and r each change by 4–7 × LO per unit ln τ_cut in
  the hard bins and cancel to the percent level, so a 0.3% effect in either
  suffices.
- Candidates: the P2B extra terms at NNLO (ZEUS uses P2BSLICE, the fixed
  point not), the jet observables, or the per-event μ = Q over the Q²/y range.
  The cluster's E set (mode 2, τ_zQ bins in the inclusive cuts, plain
  slicing, plus P2B-extra runs) separates these; the Z set (ZEUS dijets,
  P2B, fresh seeds) replicates the drift or not. Locally the precision needed
  (≈ 3500 r seeds) is out of reach.

## 2026-10-11/12 (night) — first runs with thsub: N3LO 1+1 by P2B, and the ZEUS P2B − plain term

AK: test the new submission tool (`~/bin/thsub`) with overnight runs; as one
option "some p2b n3lo jets in the lab frame to see what kind of
convergence we have". 2666 jobs (≈ 2.4k core-h) in one night; all done with
exit 0 except one r job on thA357 that vanished mid-run (process group gone,
users logged in there; thsub marked it `lost` as designed; resubmitted). The
scheduler filled all 19 thservs (nice 0, ≤ 32 of mine) and the free desktops
(nice 19) and exited at 03:01.

**ZEUS dijets: P2B − plain (`P2BEXTRA=1`; vi 100, kp 100, r 600 seeds,
`runs/zx`, `zx_combine.py`).** Its change per unit ln τ_cut in % of the bin's
LO is tiny in the drifting hard bins (p̄_T 30–60: −0.04 ± 0.03, −0.02 ±
0.02, 0.00 ± 0.01 in the windows 5e-3→5e-4, 2e-3→2e-4, 1e-3→1e-4; p̄_T 22–30
≤ 0.1; m₁₂ 65–120 ≤ 0.23 ± 0.14) against the NNLO drift of +2.4 … +4.6;
total −0.7 ± 0.4 against +1.3 ± 0.5. **The drift is not from P2B**: it sits
in the plain-slicing parts (b2, vi, kp, r) of the jet bins. Left: the jet
observables (slow power corrections in hard kinematics) or μ = Q per event
over the Q²/y range — the cluster's E set (plain, τ_zQ, integrated) separates
them.

**N3LO 1+1 lab-frame jets by P2B (`runs/n3p`, `n3p_combine.py`).** N3LO
coefficient = disorder inclusive `-n3locoef` (lab11 analysis, 16 seeds) +
sliced21 b2 (100) + nlo31 vi (200), kp (100), r (1200) in mode 3, plus the
technical-cut decades rc1 (150), rc2 (100) below 1e-4; first use of the
NNLO 2+1 parts in mode 3 (smoke tests finite).
- The total row is exact at every τ_cut (−0.40% of σ_LO(≥ 1 jet) = the
  inclusive N3LO coefficient; P2B leaves the total unchanged).
- Per bin, plain seed means: consistent with flat (pulls against 2e-4 within
  ≈ 2.5σ), but the errors are 0.5–4% (≥ 1 jet) and up to ±14% of LO in the
  central rapidity bins from 1e-3 on: ≥ 1 jet −1.4 ± 0.5 (2e-2) … −7.9 ± 4.1
  (2e-4).
- The errors come from r: single seeds 23–30 seed-σ out (seed 65035:
  −14552% and +14488% of LO in the two neighbouring central y bins = a
  large R − D cancellation binned into adjacent bins); b2, vi, kp are clean
  (0.02–0.3%). With r trimmed by 1% on each side (diagnostic only, biased):
  ≥ 1 jet −0.6% of LO, flat from 2e-2 to 1e-3 (±0.1–0.5), then −1.3, −1.5,
  −2.6 ± 1.1 at 1e-4; forward y and ≥ 2 jets drift up (to +5%).
- Next: replay of seed 65035 with P2BDEBUG (running) to see the tail events;
  then decide whether the tails are genuine mis-binning (more statistics, or
  smoother binning) or a defect like the psmc one of 9 Oct.
- **The r tails in mode 3 are a VEGAS-target effect; `VTARGET=0` (12 Oct,
  morning).** Replays of seed 65035 (bit-identical, also with a per-cell
  debug trigger) show no single event above 1e10: the ±1.95e5 pb in two
  neighbouring y bins come from the VEGAS weight. Mode 3 adapts only to the
  ≥ 1 jet row, to which an event whose leading jet migrates between y bins
  contributes nothing, so VEGAS never adapts there and the rare points carry
  large grid weights. New option `VTARGET=0`: VEGAS adapts to the sum of
  |cells| at the smallest τ_cut over all rows (sampling only; default
  bit-identical, checked). Same 167 seeds, default against VTARGET=0, r at
  2e-4 in % of σ_LO(≥ 1 jet): p_T 5–8 −159 ± 28 → −138 ± 2.0; p_T 8–11 +76 ±
  30 → +53 ± 5.7; y −0.5…0 +24 ± 25 → +69 ± 7.6; y 0…0.5 −119 ± 24 → −142 ±
  7.0; y 0.5…1 −95 ± 7.2 → −104 ± 1.6; ≥ 2 jets −345 ± 18 → −339 ± 2.8; ≥ 1
  jet −133 ± 9.2 → −131 ± 6.0; same CPU per job, all bins agree within 2σ.
  Some bins keep ~10σ seeds. Next: the N3LO r set with VTARGET=0; the same
  question for mode 2 (ZEUS, target = total row).

## 2026-10-12/13 (night) — fixed point in hard kinematics; E P2B extras and the VEGAS target

AK: after REPORT-2026-10-12 (Z/W jets NNLO rises in ln τ_cut, W 7–11σ; E plain falls),
swap tonight's runs; the cluster gets UPDATE 12 Oct (Ex6, G6). AK: add statistics if useful.

**Fixed point x = 0.1, Q² = 5000 (`runs/fph`; plain, mode 1, τ_zQ bins; photon and W⁻).**
- NLO (b1 50 + lo 50) against DISENT (photon: 40 × 3M tonight; W: `runs/ccw1`): converged within ≈ 1σ in
  every bin; totals at 1e-5 −0.24 ± 0.30% (photon), +0.37 ± 0.30% (W).
- NNLO (b2 50, vi 60 / 37, kp 50, r 500 / 700), % of the NLO, τ_cut 5e-3, 5e-4, 2e-4, 1e-4: photon total
  50.1, 51.0, 55.0, 54.4; W total 50.1, 54.1, 53.0, 51.0; slopes per unit ln τ_cut (windows 5e-3→5e-4,
  2e-3→2e-4, 1e-3→1e-4) photon +0.4 ± 0.5, +0.8 ± 0.6, +1.1 ± 0.8; W +1.8 ± 0.4, −0.2 ± 0.5, −1.5 ± 0.7.
  The upper τ_zQ bins (0.2–0.5) rise by ≈ 2–5% of the NLO per unit ln in both photon and W (2–4σ in
  single windows): not W-specific. Below 1e-4 the technical-cut decades are not included here.
- **The ZEUS W drift (≈ 25% of the NNLO per unit ln) is absent at this fixed point**, although the ZEUS
  W drift is 8σ in the Q² bin 2000–5000 that contains it, and the W drift is significant in every Q²
  bin ≥ 1000 (so not a √ŝ/Q effect at low Q²). It grows with p̄_T (−0.9 … −2.0% of LO per unit ln).
  → it comes with the jet observables. Test running: W with the ZEUS jet selection in the fixed-point
  window (ZFIX=2 at (0.01, 400), `runs/zw2`).
- Set-up errors caught before CPU was used (zw2): nlo31 in ZFIX=2 takes the window centre from its
  arguments (`2 0.01 400`, not `2 0 0`: empty window); W needs the W beam tables (`zeus5w`). 560 jobs
  cancelled while pending and resubmitted corrected.

**E P2B extras and the VEGAS target (`runs/exv` VTARGET=0, `runs/exw` VTARGET=-1, against the cluster's
bin5 default).** `VTARGET=0` (smallest τ_cut, all rows; right for the above-cut parts of mode 3) is the
wrong target for the P2B extras, θ(T ≤ τ_cut): the smallest τ_cut is the smallest region; per-seed spread
worse than the default in many cells (e.g. r-x bins 0.1–0.3 at 2e-4: 2533 vs 592, 2551 vs 56 pb). New
`VTARGET=-1` (all cells): per-seed spread against the default lo-x 10–40× smaller, kp-x 3–10× (1e-3),
r-x 2–4× (1e-3) and 1.8–20× (1e-4), vi-x mixed (lowest bin 5–10× better, two bins ≈ 2× worse).
Cluster instructions corrected (Ex6: default target, not VTARGET=0; e57dac8 pushed alone); the
`VTARGET=-1` code (03f6981) awaits AK's OK.
- **W with the ZEUS jet selection at the fixed-point window (ZFIX=2 at (0.01, 400), P2BSLICE=1;
  `runs/zw2`; b2 50, vi 27 of 60, kp 50, r 400):** the NNLO in the window rises from 3.04e-8 to
  4.80e-8 pb between τ_cut 5e-3 and 1e-4; slopes per unit ln τ_cut +4.5e-9 ± 1.5e-9, +6.5e-9 ± 2.1e-9,
  +4.8e-9 ± 2.2e-9 (3.0, 3.1, 2.2σ), i.e. 15–21% of the NNLO per unit ln, as in the ZEUS W set (≈ 25%);
  carried by the soft jet bins (p̄_T 8–15, m₁₂ 20–30: 2.4–3.2σ; hard jets are rare at Q² = 400).
  The same point in τ_zQ bins (F) is flat outside the lowest, power-corrected bins, and so is the
  fixed point (0.1, 5000) in τ_zQ bins for W and photon. **So the drift comes with the jet observables,
  also at fixed kinematics**: the interplay of the k_T jets (R = 1, Breit frame) / lab-frame E_T, η
  selection with the τ₂ slicing, not the 2+1 matrix elements or the Q²/y integration (P2B excluded
  11 Oct). Not yet checked: photon jets in the same window (B photon drifts too, so expected).

## 2026-10-13 — REPORT-2026-10-13 (cluster: G6 fixed point, Ex6) read at MPP

**G6 (x = 0.1, Q² = 5000, τ_zQ bins, photon and W⁻, 2–4× the MPP `runs/fph` statistics, technical-cut
decades included).** NLO converged to DISENT from τ_cut 1e-4 (totals at 1e-5 +1.6σ W, +0.5σ photon).
NNLO consistent with `fph` (totals in % of the NLO at 5e-3, 5e-4, 2e-4, 1e-4: W 50.1, 53.2, 52.0, 53.3
vs 50.1, 54.1, 53.0, 51.0; photon 49.9, 52.8, 53.0, 51.9 vs 50.1, 51.0, 55.0, 54.4; window slopes agree
within ≤ 2.2σ). Total bin flat below 1e-3 for both bosons (|pull| ≤ 1.8 in 1e-3→1e-4 and 1e-4→1e-5, the
latter new). The two upper windows (5e-3→5e-4, 2e-3→2e-4) have 4–5σ slopes of either sign, but 5e-3
and 2e-3 are above the window rule τ_cut ≤ 5e-3 × lower bin edge for every bin except 0.4–0.5.
**Caveat on the 13 Oct statement "the drift is absent at the fixed point"**: it holds for the total.
The upper bins (0.2–0.5) rise as τ_cut falls in both bosons by ≈ 1–4% of the bin NLO per unit ln
(≈ 2–7% of the bin NNLO), also in 1e-3→1e-4, inside the window rule (average of the six upper bin
slopes −1.9 ± 0.5% of NLO, ≈ 4σ, windows not independent), and are not resolved in 1e-4→1e-5 (errors
3×). This is the 2–4σ upper-bin rise of `fph`, now confirmed. It is smaller than the jet drift (15–25% of
the NNLO per unit ln) but has the same sign, and a constant slope in ln τ_cut is what a wrong
log coefficient below the cut would give. Open; it fits with the jet-observable dissection (next).

**Ex6 (E P2B extras with VTARGET=0).** It confirms `runs/exv`: scatter 2–7× larger than bin5 at
τ_cut ≥ 3e-5, smaller at 1e-5 (the target column); NNLO_P2B at 1e-5 83 ± 35 pb vs NNLOJET 59 ± 12.
With `VTARGET` set, the printed VEGAS iterations are on the target's scale (sum of |cells|), so the
cluster's iteration-based outlier flag is not comparable (its §3.1 caveat explained). Cluster UPDATE
13 Oct (02b9f22): Ex7 = Ex6 with `VTARGET=-1` on bin7, combined with Ex6.

Photon with the ZEUS jet selection in the window (`runs/zp2`): started by mistake at 08:18, during
the day; AK: tonight. Killed at ≈ 08:40 (thsub, by PID). **Correction**: not resubmitted locally; the photon jets in the
window went to the cluster as J7 (UPDATE 13 Oct (c), 25b11be), together with N7 (NNLOJET W top-up, (b)).

## 2026-10-13 — Ex7 pilot (cluster, bin7, `VTARGET=-1`) read at MPP

bin7 built from e6d834c (contains 03f6981), bit-identical to bin6 without `VTARGET`. Pilot: 20 r-x seeds
(9101–9120), all finite, 0.65 core-h per job (Ex6 1.43). Total row (sd / 1.48 MAD, pb) at 1e-3 and 1e-4:
Ex7 481/267 and 442/62.5; Ex6 pilot (20 seeds, VTARGET=0) 73.5/38.4 and 27.1/11.3; Ex6 (2000 seeds)
3942/30.4 and 3550/12.1. Medians Ex7 −119/−38, Ex6 −6.3/−1.5 (≈ 2σ at 1e-4 with 20 seeds).
**Correction**: the record of 13 Oct called `VTARGET=-1` best for the extras (local `runs/exv`, against
the default target). Against `VTARGET=0` the core of the per-seed distribution is 5–9× wider with −1:
the target sum(|cells|) over all rows and τ_cut (≈ 2e4 against results of ≈ 100) is dominated by large
migration weights that cancel between cells, so VEGAS adapts to them rather than to the cells we need.
Whether −1 is still better for the mean depends on the tails, which made Ex6's sd 100× its core and
which 20 seeds cannot see; the pilot criterion I wrote (robust scatter below Ex6) was therefore the wrong
test. Cluster UPDATE 13 Oct (d): pilot to 200 seeds, feed Ex7 only if sd² × cost beats Ex6, the
outlier fraction is no larger and the means agree; otherwise Ex6 + bin5 are the final E extras.

## 2026-10-13 — Ex7 pilot, 200 seeds (cluster 0bb89c8): `VTARGET=-1` IS better for the extras

**Correction of the previous entry** ("against `VTARGET=0` the core is 5–9× wider"; the 20-seed pilot).
With 200 r-x seeds (9101–9300) the total row is −110.6 ± 25.0 pb (1e-3) and −35.7 ± 14.3 (1e-4), sd 353
and 203, 0 and 0.5% of seeds beyond 20 × (1.48 MAD). Ex6 (2000 seeds, `VTARGET=0`): 20.8 ± 88.1 and
−77.4 ± 79.4, sd 3942 and 3550, 3.9% and 3.8% beyond; its first 200 seeds give +666 ± 741 and −832 ± 786.
The core of Ex7 is wider (1.48 MAD 146 and 67 vs 30 and 12), but the mean is what counts, and Ex6's mean is
carried by rare large seeds: per core-h (0.65 vs 1.43 per job) Ex7 is ≈ 270× more efficient. Means agree
(pulls −1.43, +0.52); the medians (−139, −26 vs −6, −1.5) differ only in the shape of the distributions.
So the 13 Oct conclusion from `runs/exv` (−1 best for the extras) stands. The cluster held production back
on the 1e-3 pull (its rule: marginal = do not feed); cluster UPDATE 13 Oct (e): feed Ex7 (2000 r-x seeds,
expected ≈ 7.5 and 4.3 pb on the total row), report Ex7 alone and Ex6 + Ex7.

## 2026-10-14 — REPORT-2026-10-14 (cluster: N7, J7, Ex7) read at MPP: the jet drift is a bug of ours

**N7 (NNLOJET W ZEUS dijets topped up):** NNLO total 0.148 ± 0.015 pb (was 0.176 ± 0.056). Ours (W5)
crosses it at τ_cut ≈ 1.3e-3; pulls at 2e-4, 1e-4, 1e-5: +4.1, +5.3, +4.4 (ptavg 22–30 up to +6.9). Since
our NNLO keeps rising at a constant rate in ln τ_cut, the agreement near 1e-3 is a crossing, not
convergence: **the small-τ_cut values are wrong**, and the constant slope points to a wrong single-log
(ln τ_cut) term at NNLO in the 2+1 slicing, depending on the Born kinematics or the channel (it nearly
cancels in the inclusive G6 total at x = 0.1, Q² = 5000). NNLOJET Z 12.9 ± 1.0 (one RRa_2 seed of 183.5 pb
dominates); ours 1.5–2.7σ above at τ_cut ≤ 1e-3, the same pattern at lower significance.
**J7 (window (0.01, 400), ZEUS jet cuts):** the photon drifts like W; total slope (% of |NLO(1e-5)| per unit
ln, 1e-3→1e-4) −21 ± 7 (γ), −27 ± 10 (W), ratio 1.26 ± 0.63; similar bin by bin. So neither the W couplings
nor ew31 are involved. The cluster added b1/lo for the NLO denominators (the instructions had omitted them).
**Ex7:** E NNLO_P2B (Ex7) at τ_cut 1e-4 / 1e-5: 65 ± 16 / 67 ± 34 against NNLOJET 59.5 ± 12; NLO_P2B → NNLOJET
from 2e-4. At τ_cut ≥ 2e-3 the extra parts of Ex7 differ from Ex6 and bin5 by up to 29σ (kp-x vs bin5).
Reading: bin5/Ex6 hardly adapt to the large-τ_cut cells (their per-seed medians are near zero with rare large
seeds; Ex7's distributions are close to Gaussian), so their means there are unreliable and their errors
underestimated. Checked by X8 (VTARGET=53, the total row at 5e-3). Cluster UPDATE 14 Oct: hold production of
ours until the drift is found; X8. Plots: https://claude.ai/artifact/XHDE6je81vPkdT25YJVseE

## 2026-10-14 — pointwise RG check of the below-cut NNLO cumulant (lp21), and J8

AK: bug hunt, (1) pointwise check locally, (3) ZEUS jet cuts at the G6 point on the cluster.
**(3) J8:** sliced21 takes the window centre from `ZWIN="x Q2"` (ZFIX = 1, 2; default (0.01, 400)).
Checked locally (`work/.../nnlo21/zwbuild/chk`): b1 (20k × 2, seed 7) without ZWIN, with
`ZWIN="0.01 400"` and with the zw2 binary: all ZCELL rows identical; with `ZWIN="0.1 5000"` the log
prints x 0.1, Q² 5000 and nlo31 lo `2 0.1 5000` the matching window; b1 + lo cancel per bin as at
(0.01, 400). Cluster UPDATE 14 Oct (b) (37301af): J8 = J7 PH with the window at (0.1, 5000).
**(1) What sets the τ_cut slope below the cut.** The drift is a constant slope in ln τ_cut, so the
suspect is the L¹ coefficient of the O(α_s²) cumulant (L = ln τ_cut). Per LP function, its two-loop
single log is fixed by the RGE from its one-loop terms, β0, Γ1/Γ0 and its two-loop non-cusp anomalous
dimension. Check (`work/disorder-comparisons/nnlo21/rgchk`: `rg_lp21.f90` dumps the distribution
coefficients lp21 passes to assemblejet; `rg_check.py` goes to Laplace space and compares the
l³, l², l¹ coefficients of ln F~ at two loops with the RGE):
- MCFM jet functions (control, RG-consistent by construction): l³, l² exact (1e-15); the implied
  two-loop non-cusp of the quark jet is the literature γ^J_1 (77.35 in α_s/4π).
- Soft functions as lp21 calls them (qgq and qag, legs 1, 2 = the jets back to back, 3 = the beam,
  all six I_ij,m), 19 angles each: l³ and l² exact (≤ 5e-13); the remaining two-loop non-cusp
  term (l¹ minus β0 × constant minus Γ1/Γ0 × the one-loop angle logs) is 80.973926 for every angle
  and both channels. So the soft function has no kinematics-dependent single-log error.
- Beam functions: the test does not apply as written (the one-loop terms depend on μ through the
  PDF ratios); they are MCFM's, at the Born ξ.
- Assembly: MCFM's assemblejet (generic convolution of the four slots; jet slot in α_s/4π, beam
  slots in α_s/2π; lp21 passes the second jet as beam b with J1/2, J2/4, as required). The
  measure of nlo31 (jets' rest frame, Q_i = 2E_i, beam term p·k/p·u, jet term u·P_J − |P_J|)
  matches lp21's inputs (E_a, E_2, E_3, y_ij in the jets' rest frame).
**So far no error on the below-cut side that would give a kinematics-dependent ln τ_cut term.** Not
checked pointwise: the beam tables in ξ (MCFM functions; table accuracy checked earlier) and the
above-cut 3+1 NLO near the 2+1 limit. Note the size: in J7 PH the drift (≈ 0.5e-6 per unit ln in the
window units) is ≈ 0.1% of b2's own slope (≈ 4e-4 per unit ln; b2 = 1.1e-3 against an NNLO of 6e-6 at
1e-4), i.e. a 0.1% effect in the b2 ↔ (vi + kp + r) cancellation. A sampling bias of that size in r
(error 0.06% of r per τ_cut) would do it; Ex7 has just shown that cells VEGAS does not adapt to can
be off by many σ with plain means.

## 2026-10-14 — cut-design rows for the ZEUS selection (ZCUTS), to test fiducial power corrections

Literature (docs/literature-fiducial-power-corrections.md, 14 Oct): a symmetric per-jet E_T cut on two
jets balanced at Born level gives √τ·L^k power corrections, whose slope in ln τ_cut is almost
constant over our range. P2B should remove them, and P2B is on in our ZEUS runs, but the lowest
p̄_T bin, which carries the J7 drift, was never examined for it. New option `ZCUTS=1` (nlo31 and
sliced21, mode 2, also with the ZFIX=2 window): six extra rows after the 15 ZEUS rows, each as total
and p̄_T 8–15: `et4` (jets E_T > 4 GeV, p̄_T > 8, η cuts and m12 > 20 as ZEUS: the per-jet cut is far
from the Born jets), `noeta` (ZEUS E_T cut, no η_lab cuts), `et4ne` (both). nobmax 16 → 21; the
VEGAS value is unchanged (sliced21: an event only in the extra rows gives 0 to VEGAS).
Checks (`work/.../nnlo21/zcbuild/chk`, short runs, seed 7, window (0.01, 400)): without ZCUTS, sliced21
b0, b2 and nlo31 lo, kp, r (window) and r (full ZEUS, plain) are identical to the 37301af binaries
(whole output); with ZCUTS the 15 ZEUS rows and RESULT are unchanged, and at Born level (b0, b2)
`et4` equals the ZEUS rows exactly and `et4ne` equals `noeta`, as it must (balanced Born jets).
If the drift is fiducial it should be absent, or much smaller, in `et4`.

## 2026-10-15 — REPORT-2026-10-15 (cluster: X8) read at MPP: Ex7 right, Ex6 and bin5 extras biased at large τ_cut

X8 = the Ex7 extras with VEGAS aimed at the total row at τ_cut 5e-3 (`VTARGET=53`; lo-x 200, vi-x 250,
kp-x 200, r-x 400 seeds 8001–, ≈ 500 core-h). Total row, τ_cut 2e-2 … 2e-4: X8 − Ex7 pulls lo-x −1.0 … −0.2,
kp-x −2.2 … 0.0, r-x −1.3 … +2.0 (X8's r-x has single seeds of −62053 and −16146 pb), vi-x +1.4 … +2.0 at
every τ_cut (one correlated offset of 1–1.7 pb, ≈ 2%; not pursued). Ex6 and bin5 are far off at
τ_cut ≥ 2e-3 (X8 − bin5 up to +28.5σ in kp-x, +24.6σ in vi-x; X8 − Ex6 up to +7.0σ). **So the reading of
14 Oct holds: cells VEGAS was not aimed at come out wrong by many of their own σ** (most seeds miss the
region carrying the mean; sd/√N underestimates the error), with no code error involved. Consequences:
the Ex7 values of the E P2B extras stand (NNLO_P2B 1.1 ± 15 at 1e-3, 65 ± 16 at 1e-4, 67 ± 34 at 1e-5
against NNLOJET 59.5 ± 12); the Ex6 and bin5 extras at τ_cut ≥ 2e-3 are not to be used. For
production: aim VEGAS at every reported cell (`VTARGET=-1`) or check cells with a dedicated target, as
here. This also makes the J7v test (sampling bias in the plain parts of the jet runs) more relevant.

## 2026-10-14 — beam tables excluded as the cause of the jet drift

Twice-finer beam grid (0.01 in ln(ξ/(1−ξ)), 0.05 in ln Q, 54 nodes; `work/.../nnlo21/beamgrid5f`, built
with thsub): tabchkg max deviation / max coefficient 1.1e-4, 9.0e-4, 6.6e-5, 2.5e-3, 2.3e-4 (up q, q̄,
down q, q̄, gluon) against 1.1e-3, 5.3e-3, 4.0e-4, 6.7e-3, 3.3e-3 for the production grid. J7 photon b2 in
the (0.01, 400) window (5M × 6, seeds 4001–4020, same binary, both grids; `runs/tabchk`), per-seed
difference fine − old: total at 1e-4 +5.7e-10 ± 0.6e-10 against b2 = 1.10e-3 (5e-7 relative); slope per unit
ln τ_cut −7.7e-11 (5e-3→5e-4), −2.2e-10 (1e-3→1e-4), −6.0e-10 (1e-4→1e-5), against the J7 drift of
≈ −5.6e-7 per unit ln: **2500× smaller**. Same picture in p̄_T 8–15 and m12 20–30. The tables are not the
cause. A HOPPET-based beam function (AK, 14 Oct: build such convolutions with hoppet) stays as a clean-up
for production.

## 2026-10-15 — REPORT-2026-10-15-J7v and -J8 read at MPP: no statistical sampling bias; the drift is in the Born-level jet kinematics

**J7v** (J7 photon window; vi, kp, r with VTARGET=-1 and VEGAS_EQUAL=2, all parts rerun on bin9 with
ZCUTS; r 1581/1600 at the time of the report): the drift stays. Total slope (% of |NLO(1e-5)| per unit
ln τ_cut) −12.4 ± 4.1, −15.5 ± 5.2, −13.0 ± 6.1 (V9) against J7 PH −17.5 ± 6.1, −17.6 ± 6.8, −21.0 ± 7.1.
Every part agrees with J7 PH cell by cell within 1–2σ; per-seed vi, kp, r are close to Gaussian (sd ≈
1.0–1.5 × 1.48 MAD; beyond 20 MAD: none in vi, kp, 1/1600 in r). **So a statistical VEGAS bias (target
cell, 1/σ² iteration weights, heavy tails) is not the cause** — unlike the P2B extras (X8), these parts
are well sampled. Not tested: a deterministic coverage gap of the phase-space generator (psmc).
**J8** (ZEUS jets at (0.1, 5000), photon, bin9, ZCUTS; r 1569/1600): drift at 12–16σ in every bin
(total −43.9 ± 2.8, −51.8 ± 3.8, −57.3 ± 4.8 % of |NLO| per unit ln; NLO is itself a 300× cancellation
there). G6 at the same point without jets is flat: not an x or gluon-content effect.
**Cut-design rows:** et4 − ZEUS and et4ne − noeta are τ_cut-independent to all printed digits (J7v and
J8). **Correction to my test design (14 Oct):** with P2B this holds by construction (the τ_cut
dependence enters only through O(projected Born), and the projected jets are balanced, so et4 and ZEUS
select the same Born events); ZCUTS could not have shown a cut-induced drift with P2B on. It does
confirm that P2B removes cut-induced terms exactly, and that the drift sits entirely in the Born-level
jet kinematics: it is a mismatch between the LP spectrum (b2) and the above-cut 3+1 spectrum binned in
projected-Born jet variables, which nearly integrates to zero over the inclusive Born phase space.
Excluded so far: boson, LP log structure (RG check), P2B projection, hadronic PCs, beam tables,
statistical VEGAS bias, cut-induced terms, x/gluon content.

## 2026-10-15 — REPORT-2026-10-15-G8 read; CH8 (channel split) set up

**G8** (J8 window, vi/kp/r with the `flat` generator, ≈ 15k core-h): flat vi and kp undersample τ_cut < 5e-4
(4× ncall moves them to psmc; flat vi up to 37σ off at 1e-5); flat r agrees with psmc at every τ_cut within
errors 5–60× larger; the slope precision asked for is missed by ×7 (flat r scatters 10× more per seed). In the
window where flat is reliable (5e-3 → 5e-4) the total slope is −86.7 ± 40.5 (flat) against −42.3 ± 2.6e-9
(psmc), 1.1σ apart: **no sign of a psmc coverage gap** (weak test; "unlikely", not excluded).
**CH8:** the J8 window split by incoming parton (pdfmask 1 = quarks, 2 = gluon) in every part. New: `sliced21
mktab`/`tabchkg` read `PDFMASK` (masked beam tables; default unchanged). Check (`work/.../nnlo21/chtab`, Q node
10): unmasked node bit-identical to `beamgrid5/zeus5_10.tab`; mask 1 + mask 2 = unmasked to 9e-11 (the z
integration tolerance); mask 1 has no gluon class, mask 2 no quark classes. Cluster UPDATE 15 Oct (b).

## 2026-10-15 — vi near the 2+1 limit: virt31 = NNLOJET pointwise down to λ = 1e-7

(AK: check vi pointwise.) `disorder-comparisons/dis31_nnlojet/harness_v31lim.f90`: the one-loop comparison of
3 Oct (virt31 against NNLOJET's B2g1Z family, finite parts equal up to one constant × tree, −5.79471429 for
q q̄ g g) at points approaching the limits of the 3+1 configuration: a soft gluon (E = λW), final-state
collinear q‖g and g‖g (m² = λW²), and initial-state collinear (angle √λ), λ = 1e-1 … 1e-7, both q q̄ g g
channels (d → d g g, g → d d̄ g), current virt31 (824ab3d) minus its closed-loop term `virt31_vloop` (absent in
NNLOJET's DIS). Result: the constant holds to ≤ 4e-10 down to λ = 1e-4 in every limit and channel, ≤ 1e-8 at
1e-5, ≤ 6e-8 at 1e-6 and ≤ 1e-6 at 1e-7. **The one-loop 3+1 matrix element is correct and numerically stable
in the region the slicing uses**, so vi's matrix element is excluded as the cause of the jet drift (not tested
near the limits: the four-quark channels, subleading there). The closed-loop term itself vanishes towards the
limits (O(0.1–0.3) of the tree at generic points).

## 2026-10-15 (night) — REPORT-2026-10-15-CH8 read: the drift is in the quark-initiated channel

J8 window split by the incoming parton (pdfmask 1/2 in every part; per channel the J8 seed counts; 8.6k
core-h). q + g agrees with J8 (340 pulls, rms 0.89, max 2.6). Total-row NNLO slopes (1e-9 pb/GeV² per unit
ln), windows 5e-3→5e-4 / 2e-3→2e-4 / 1e-3→1e-4 / 1e-4→1e-5: quark −33.7 ± 2.7, −44.8 ± 3.8, −43.7 ± 5.3,
−11.3 ± 11.7; gluon −0.7 ± 1.0, −2.7 ± 1.4, −1.5 ± 1.4, +1.5 ± 2.6. Gluon share of the slope 0.03 ± 0.03, and
≤ 0.06 in every bin. The quark NLO also carries all the NLO τ_cut dependence (gluon NLO nearly flat).
Quark-Born-specific pieces below the cut: gluon jet function, quark beam function, soft qgq, hard function of
γ* q → q g; above the cut the quark-initiated 3+1 / 4+1 channels. Also: in every set the 1e-4 → 1e-5 slope is
consistent with zero (J8 +1.9 ± 13.1, CH8 q −11.3 ± 11.7, J7v −1.4 ± 18): a drift that stops below 1e-4
would mean large quark-channel power corrections plus a constant offset against NNLOJET rather than a
wrong log term. Not decided yet.

## 2026-10-16 — hard21 excluded; a bug in nlo31's K + P (CF shadowed in `kp_conv`) is about half of the jet drift

AK: check the two-loop hard function and the I/KP operators of the above-cut NLO 3+1, and either confirm an
issue or rule it out beyond reasonable doubt. Test programs in `work/disorder-comparisons/nnlo21/{h2test,kptest}`
(not in the repo: they link NNLOJET).

**hard21: excluded.** `h2test/h2nj.f90`: hard21 (Catani scheme; a verbatim copy of hard21_eval's assembly, equal to
hard21_eval to 5e-13) against NNLOJET v1.0.2's DIS double-virtual functions (B1g2Z = makem0mn + makem1m1 with the
colour weights of sigVVDIS) at 200 random DIS points per channel (γ* q → q g, γ* g → q q̄): two-loop interference,
|A1|² and the one-loop interference agree to 3.6e-7 (quark) and 1.4e-11 (gluon). The one-loop conventions differ
only by the bookkeeping term n_f (L13 + L23)/12 (checked to 4e-15). The Catani → SCET conversion is MCFM's
(re-derived). hard21 can cause neither the drift nor an offset.

**K + P: bug confirmed and fixed.** In `kp_conv` (dis31/nlo31.f90) the local variable `Cf` (∫₀^ξ (1+x²)/(1−x) dx)
has the name of the module parameter `CF` (Fortran is case-insensitive; a local declaration legally hides the host
entity, no compiler warning). In the q → q branch both colour factors of the integral of the plus parts therefore
read the integral instead of 4/3: ipl = Cq·2A·b + g·Bf + Cq²·lp instead of CF·2A·b + g·Bf + CF·Cq·lp. Only the
diagonal quark term is affected; present since the first nlo31 (dd004f6, 2 Oct). Fix: the local is now `Cq`.
- Independent check (`kptest/kp_indep.py`, mpmath): the CS kernels written out from the general formulae
  (K̄^{ab} = P^{ab} ln((1−x)/x) + P̂′^{ab} + δ^{ab}[T²(2 ln((1−x)/x)/(1−x))₊ − δ(1−x)(γ + K − 5π²/6 T²)], P^{ab},
  the γ term [1/(1−x)]₊ + δ(1−x)), distributions as ∫ h (F(x) − F(1)), against kp_conv integrated over x
  (`kpdump`), test PDF ξ^−0.3 (1−ξ)³, ξ = 0.01, 0.1, 0.4, 0.8, Born factors (b, g, lp) = (1,0,0), (0,1,0),
  (0,0,1), (1.3,−0.7,0.45). Before: q → q off by 1–18% (b), up to 114% (lp alone, near a zero); q → g, g → q,
  g → g ≤ 2e-10. After: every channel ≤ 2e-10.
- The rest of I and KP: poles of V + I, d lp/d ln μ_F² and the RG structure were already tested
  (`dis31/tests/harness_iop31`). The sign of the γ term and the finite constants of I and K̄ (not tested so far):
  with iop31's formulae, γ* q → q at NLO gives V (−8 C_F) + I (2(γ_q + K_q) − 2π²/3 C_F) + K̄ (−(5 − π²) C_F,
  and −π²/3 C_F from writing [2 ln((1−x)/x)/(1−x)]₊ with [ln(1−x)/(1−x)]₊) + γ term (−γ_q) = −(9/2 + π²/3) C_F
  δ(1−x) and −3/2 C_F [1/(1−x)]₊: the MS-bar C₂^q. The opposite sign of g would give +3/2. born_part's channel
  assignment (q → q, g → q for a quark Born; Σ_q q → g and g → g for a gluon Born) and PDF arguments read fine.
  **I and KP are now correct beyond reasonable doubt; the bug was the only issue.**

**Correction of the 3 Oct NLO 3+1 validation** (`runs/kpfix3j`; same set-up, kp with the old and the fixed
binary, 4 seeds each): the old binary reproduces K + P = 31.106 pb (3 Oct: 31.092 ± 0.035); fixed 30.160 pb,
Δ = −0.946 ± 0.003 pb (3% of the correction). NLO correction (≥ 3 jets) 30.34 → 29.39 ± 0.29 pb against NNLOJET
RV + RR 29.99 ± 0.68 (+0.5σ → −0.8σ); Q² bins χ² 9.2 → 6.6 for 6. **The 3 Oct validation could not see the bug**
(1.3σ); "validated at the 2% level of the correction" was true of the comparison, but not a test at that level of
the K + P piece. The corrected values agree as well.

**Impact on the jet drift** (`runs/kpfix`, kp old/fixed with the same 24 seeds per set; Δkp = fixed − old; slopes
per unit ln τ_cut, windows 5e-3→5e-4 / 2e-3→2e-4 / 1e-3→1e-4):

| set | Δkp slope | NNLO drift | share removed |
|---|---|---|---|
| J8 (photon, ZEUS jets at (0.1, 5000)) | +2.31, +2.36, +2.33e-8 | −4.23, −4.98, −5.53e-8 | 55, 47, 42% |
| J7 PH (photon, (0.01, 400)) | +2.44, +2.54, +2.65e-7 | −4.67, −4.73, −5.63e-7 | 52, 54, 47% |
| W5 (W, full ZEUS) | +0.0299, +0.0299, +0.0278 pb | −0.044, −0.045, −0.050 pb | 67, 67, 56% |

Δkp has errors ≤ 1% of itself. W ZEUS total NNLO coefficient after the fix: 0.054 ± 0.009, 0.082 ± 0.014,
0.104 ± 0.018, 0.154 ± 0.045 pb at τ_cut 1e-3, 2e-4, 1e-4, 1e-5, against NNLOJET 0.148 ± 0.015 (pull −5.3, −3.2,
−1.9, +0.1; before +0.6, +4.1, +5.3, +4.4). Bin by bin at 1e-5 all within 2.3σ (before up to 4.4σ).
**A residual drift remains**: J8 −19 ± 3, −26 ± 4, −32 ± 5e-9 (7σ in the last window), J8 1e-4 → 1e-5
+18 ± 13; CH8 quark channel −11 ± 3, −21 ± 4, −20 ± 5, 1e-4 → 1e-5 +5 ± 12 (Δkp sits entirely in the quark
channel: the bug term needs a quark PDF at both ends); W ZEUS −0.022 ± 0.007 pb (1e-3 → 1e-4, 3σ); J7 1.6σ.
So the residual still sits in the quark channel and has the opposite sign to NNLOJET's offset at large τ_cut;
whether it stops below 1e-4 is Q8's question (Q8 runs r, rcorr, d2 only, so it is not affected by the bug; its
kp comes from the fixed binary).

**Corrections of the record.**
- Every NNLO 2+1 value with nlo31's kp in this notebook and in the cluster reports (B, Z, W, E, G6, J7, J7v,
  J8, CH8, G8, Ex6, Ex7, X8) carries the bug; the N3LO 1+1 P2B results (mode 3) too. NLO 2+1 results (no nlo31
  kp) are not affected. In particular the W and E "agreement with NNLOJET" at small τ_cut (10–14 Oct) and the
  15 Oct statement that the drift is "a mismatch between the LP spectrum (b2) and the above-cut 3+1 spectrum"
  in the Born-level jet kinematics: about half of the drift was this above-cut K + P error.
- The 15 Oct (night) CH8 result (drift in the quark-initiated channel) stands and is where the bug acts; the
  list of excluded causes stands, with hard21 and I/KP (after the fix) added.

**Follow-up (AK: push, rerun kp, add a test).** 40c013a pushed. Regression test: `nlo31 chk` now also runs
`check_kp`, kp_conv for q → q against closed forms with a constant PDF (f(ξ/x)/x − f(ξ) = (1 − x)/x makes every
x integral elementary; b, g, lp one at a time; closed forms checked against mpmath to 30 digits): fixed ≤ 2e-8,
with the bug restored 2e-3 … 1 in b and lp. Production output bit-identical (J8 kp run against the fixed binary
of the impact runs). Cluster UPDATE 16 Oct (b), KREDO: every kp and kp-x rerun with the fixed binary (bin11),
same seeds and commands, Δkp checked against the MPP same-seed values. The local N3LO 1+1 P2B kp (`runs/n3p`,
100 seeds) is rerun at MPP (`kp11_*`).
N3LO 1+1 P2B (`runs/n3p`) with the fixed kp (`KP=kp11 n3p_combine.py`, output `runs/n3p_result_kp11.txt`): every
bin moves by ≤ 1.0% of σ_LO(≥ 1 jet), at most 0.4σ of its error; the total row is unchanged (it is fixed by the
inclusive coefficient). The 12 Oct conclusions on that set stand.

## 2026-10-16 (evening) — REPORT-2026-10-16-Q8 read: with the fixed K + P the quark-channel drift stops below τ_cut 1e-4

Q8 (cluster, d4cbdeb): quark channel in the J8 window, CH8 q plus r 2400, rcorr 2200, d2 1600 new seeds (21001–),
kp from the fixed binary (bin11, the 200 CH8 q kp seeds rerun). Total-row slopes per unit ln τ_cut (1e-9
pb/GeV², errors from the per-seed differences): 5e-3→5e-4 −10.2 ± 2.0, 1e-3→1e-4 −16.9 ± 3.5, 2e-4→1e-4
−8.7 ± 11.0, **1e-4→3e-5 +0.8 ± 8.4, 1e-4→1e-5 +8.6 ± 7.3**, 3e-5→1e-5 +17 ± 14 (with the old kp: −33.7,
−43.7, −36.8, −34.1, −11.3, +13.7). Every cell (p̄T 8–15, both Q² bins) is consistent with zero below 1e-4; above
it the drift is 16–47% of what it was. Checks: the cluster's same-seed Δkp equals MPP's (−7.40e-8 / −1.276e-7 /
−1.653e-7 against −7.41e-8 / −1.276e-7 / −1.650e-7); old and new seeds agree (|pull| ≤ 1.2). The 6e-9 target
was missed (7.3e-9): rcorr and d2 have heavier tails than the 200/100 seeds of the sizing showed (≈ 1% beyond
20 MAD, up to 33 sd), to be remembered for production sizing.
Reading: what remains behaves like power corrections that are gone by τ_cut ≈ 1e-4; no sign of a wrong log term
is left. This does not test the value of the plateau: that needs NNLOJET (KREDO's W, Z, E). The 15 Oct (night)
statement that a slope stopping below 1e-4 would come "plus a constant offset against NNLOJET" was made with the
buggy kp and is open again. Cluster UPDATE 16 Oct (c): the KREDO report gives the per-seed slope errors and a
plateau value (mean over τ_cut 1e-4, 3e-5, 1e-5) with its pull against NNLOJET.
REPORT-2026-10-16-KREDO (a)+(b) read (cluster 604d998): the kp of J8, CH8 q/g, J7 PH, J7 W and J7v rerun with bin11
(same cmd.sh, only the binary changed). Same-seed Δkp equals MPP's (J8 −7.40e-8 / −1.276e-7 / −1.654e-7 at
1e-3 / 1e-4 / 1e-5; J7 PH −9.14e-7 / −1.533e-6 / −2.190e-6, 0.8% from MPP at 1e-3); CH8 g Δkp = 0 exactly. Corrected
total-row slopes (% of |NLO(1e-5)| per unit ln, 5e-3→5e-4 / 2e-3→2e-4 / 1e-3→1e-4 / 1e-4→1e-5): J8 −20.0, −27.3,
−33.4 (7σ), +19 ± 14; J7 PH −8.5, −8.0, −11.0 (1.6σ), −16 ± 24; J7 W −8.2, −6.4, −12.0 (1.2σ), +13 ± 18; J7v
≤ 1.2σ. (c), (d) (Z, W, E, B, G6, extras) are done; their report follows with UPDATE (c)'s plateau values.

## 2026-10-16 (20:00) — REPORT-2026-10-16-KREDO-cd read: with the fixed kp, W, Z, γ and E agree with NNLOJET for τ_cut ≤ 1e-4; correction: G6 was flat only by cancellation

KREDO (c)+(d) (cluster 0825bf3): kp rerun with bin11 for bin5 Z, W, E (and E kp-x), B (from the 11 Oct archive), G6
photon and W, and kp-x of Ex6, Ex7, X8. Same-seed Δkp for W5 −0.1034 / −0.1676 / −0.2024 pb at 1e-3 / 1e-4 / 1e-5
(MPP −0.1043 / −0.1683 / −0.2029); the kp-x changes are within about 2σ of zero (P2B extras unchanged within errors).
Pulls of the total against NNLOJET at τ_cut 1e-4 / 3e-5 / 1e-5 (before the fix):
- W vs N7 (0.148 ± 0.015): −1.9 / −1.6 / +0.1 (+5.3 / +0.5 / +4.4); 15 W bins rms 1.5 / 1.3 / 0.86, worst −2.6 / −2.3 / +2.3;
- Z vs N7 (12.92 ± 1.03): −1.0 / +0.6 / +0.1 (+1.7 / +2.7 / +2.0);
- B (γ) vs NNLOJET (11.12 ± 1.2): −0.5 / −0.5 / −1.7 (B at 1e-5: 0.9, error ≈ 6);
- E P2B (bin5 extras) +0.2 / +0.8 / +0.2; E plain −0.3 at 1e-5 (P2B needed above, as before).
Above 1e-4 the values now lie below NNLOJET (W −4 to −8σ at τ_cut ≥ 5e-4, Z ≈ −2σ): the sign has turned.
**G6 (x = 0.1, Q² = 5000, τ_zQ bins, no jets)**: with the fixed kp the total-bin slopes are positive above 1e-4,
G6P +0.77 / +2.80 / +2.36% of the bin NLO per unit ln (3.1, 4.0, 5.8σ; 5e-3→5e-4 / 2e-3→2e-4 / 1e-3→1e-4), G6W
+0.78 / +3.48 / +1.76 (3.0, 10.5, 4.0σ); largest in the 0.05–0.1 bin (+9.9σ, +10.5σ at 2e-3→2e-4); 1e-4→1e-5
consistent with zero (−0.8 ± 1.3, +0.3 ± 0.9). The Δkp slope there is +2.0% of the bin NLO.
**Correction of the record:** the 14–15 Oct statements that G6 "without jets is flat", hence "not an x or
gluon-content effect" and that the drift is specific to the jet kinematics, were made with the buggy kp, whose
slope cancelled G6's drift. With the fix every set (jets or not) drifts above τ_cut ≈ 1e-4 and is flat below,
which fits power corrections that are gone by 1e-4. Below 1e-4 every set agrees with NNLOJET within about 2σ.
The ZEUS-wide sets have large errors at 1e-5, so the quantitative plateau test is UPDATE (c)'s (cluster UPDATE
16 Oct (d): deliver it as an addendum, no jobs, then hold for AK).

## 2026-10-16 (20:20) — REPORT-2026-10-16-PLATEAU read: the plateau τ_cut 1e-4 … 1e-5 agrees with NNLOJET at 1–2σ; tag nnlo21-2026-10-16

Plateau P = mean of the NNLO at τ_cut 1e-4, 3e-5, 1e-5 (per-seed averages, bin11 kp), cluster bb46374, no new jobs.
Totals against NNLOJET: Z 13.0 ± 2.0 vs N7 12.9 ± 1.0 (+0.0σ); W 0.087 ± 0.032 vs N7 0.148 ± 0.015 (−1.8σ); γ (B)
6.7 ± 3.1 vs 11.1 ± 1.2 (−1.3σ); E P2B 50 ± 20 (Ex7 extras), 73 ± 30 (bin5 extras) vs 59 ± 12 (−0.4σ, +0.4σ). Per bin
(15 ZEUS bins, 6 τ_zQ bins) rms pull 0.7–1.4, 0–2 bins beyond 2σ per set, largest γ m12 65–120 −3.3σ (of ≈ 90
comparisons). NLO plateaus with P2B within 1.7σ (Z −1.7, γ −1.4, W +0.2, E P2B −0.4/−1.2). Slopes below 1e-4 within
1.9σ of zero in every set with per-seed errors (e2); above 1e-4 significant drifts in J8, CH8 q, Q8f, E plain, E P2B
Ex7, G6 (3–6σ). So: with τ_cut ≤ 1e-4 the NNLO 2+1 agrees with NNLOJET at the level of our errors (15–45% of the
NNLO coefficient). AK (16 Oct): tag this state (not a release) and pin down 1e-5 after the efficiency work; first
clean-up, then profiling and efficiency. Tag `nnlo21-2026-10-16` at 85f036a (pushed). Cluster on hold.

## 2026-10-17 — clean-up of nlo31 and sliced21 (step 2 before profiling and efficiency); default output unchanged

AK (16 Oct): clean-up first, then profiling, then efficiency; the list agreed with AK:
Removed (questions answered): ZCUTS rows (nlo31, sliced21), P2BDEBUG, TECHSLICE, the P2BTMIN / DIPGARB / P2BDROP
overrides (defaults kept as parameters), the VEGAS slice (argument 11, now read and ignored so that the positions
stay), the logmap generator, ZFIX = 1, the debug prints and unused-variable silencers, dated history comments,
dis31/validation/dispatch.sh (thsub replaced it). Modes 1–3 now stop without psmc (G8's `flat` answered); mode 0
keeps its uniform map (NLO 3+1 regression). Kept: modes 0–3, EW31, P2BSLICE, P2BEXTRA, TECHMIN/TECHREF/TECHDIFF,
VTARGET, VEGAS_EQUAL, ZFIX = 2 (windows) and 3 (E), ZWIN, pdfmask/PDFMASK, c1/c2, mktab/tabchk(g), chk.
nlo31.f90 1881 → 1758 lines, sliced21.f90 535 → 523.
**Regression** (`work/disorder-comparisons/nnlo21/cleanup`: `build.sh`, `regress.sh`, `compare.sh`): 104 short
fixed-seed configurations from the tag nnlo21-2026-10-16 (nlo31 lo/vi/kp/r in modes 0–3, photon/Z/W, ZEUS, J8
window, E plain and extras, rcorr/d2, VTARGET, VEGAS_EQUAL, pdfmask; sliced21 b0/b1/b2 in modes 1–3, Z, W, E, J8,
pdfmask, c1, c2, mktab, tabchk(g); `nlo31 chk`; lp21, hard21, born21 and all dis31 harnesses), ≈ 7 min with 6 in
parallel. Tag against tag: 91 identical, 15 not reproducible by construction (harnesses and tabchkg draw unseeded
random points). Cleaned against tag: every kept configuration bit-identical (outputs and beam tables); the 12 that
differ are the removed features, each as expected: the overrides now give exactly the default run, ZCUTS runs keep
their ZEUS rows bit-identical, ZFIX = 1 gives the full ZEUS result, logmap/uniform/flat in modes 1–3 stop with
"modes 1-3 need psmc"; the random-point harnesses pass at the same level (≤ 4e-13).

## 2026-10-17 — profiling (step 3) and the first two speed-ups of r

**Error budget** (J8 window, cluster per-seed outputs, CPU per job from the J8 pilots; share of Σ_p sd_p·√CPU_p, which
sets the error at fixed CPU with optimal seed counts): r 84–91%, rcorr 5–10%, d2 0–4%, vi 2–3%, b2, kp < 1%. Heavy
tails: sd/(1.48 MAD) = 1.5 for r at 1e-5, 3.4–3.9 for rcorr and d2. At the present speed an NNLO error of 1e-9
pb/GeV² (≈ 5% of the J8 NNLO) at τ_cut 1e-5 needs ≈ 3.5M core-h (1.4M without the tails). The same table for the
ZEUS-wide sets is requested from the cluster (UPDATE 17 Oct, analysis only).
**Profile** (callgrind, nlo31 r, J8 window, photon, 4000 points, 1693 real events): real_evalv over ≈ 10 flavour
structures 88%; in it the dipoles' mapped Borns (born31, colour/spin correlated) 35% (complex abs alone 13%),
accept() for every dipole of every structure 28% (216 calls per real event; tau2cm 15%, ZEUS jets 9%), me41 24%
(four-quark channel 20%); psmc 9%. vi: 96% in the BDK one loop (but 2–3% of the budget). Page:
https://claude.ai/artifact/6GPgL8pPXy9hqs44B9LYNg
**Change 1** (8869854): within an eval_sg, a mapped dipole Born with exactly the same momenta reuses accept()'s
F and degen (the CS mappings do not depend on the flavours). Bit-identical (91/91 deterministic regression
configurations); r in the J8 window 1.42× (photon), 1.33× (Z), 1.45× (W) faster.
**Change 2** (this commit): |z|² as Re² + Im² instead of abs(z)**2 in born31 (no hypot); r a further 1.25×
(photon), 1.19× (Z), 1.21× (W). Not bit-identical. **Lesson:** after VEGAS adapts, last-bit changes move the
grid edges in regions of almost no weight (rebin divides by near-zero bin weights), so iteration 2 draws a
different, statistically equivalent sample: outputs moved by ≈ 0.1% against errors of ±3%. Non-bit-identical
changes are therefore validated before adaptation (regress1.sh: the suite with one iteration): 87/89 deterministic
outputs identical to all printed digits; n3_rc1/rc2 (mode-3 technical-cut differences) differ by ≤ 2e-10 of the
largest value in each row (cancellations). Together r is faster by 1.74× (photon, 20.1 → 11.5 s), 1.57× (Z, 41.1 → 26.2 s) and 1.74× (W, 27.1 → 15.6 s).
**Infrastructure:** the session ran on thA371 (not thA371a); thA371a rebooted at ≈ 09:14 (killing a regression
pass); since then thA371's NFS view of /home/thA371a shows the repo, ~/work and ~/cernbox as owned by nobody
(writes refused), while thA371a itself is fine. Work continues over ssh thA371a.
**Change 3** (afternoon 17 Oct, after the NFS outage; thA371a again): results reused for identical arguments
within an event. born31_any keeps the last 128 evaluations (key: momenta, flavours, spin-correlation vector,
ew31 settings, born31_swapc; quick reject on the energy of parton 2): the CS mappings do not depend on the
flavours, so the dipoles of different flavour structures map onto the same Borns. me41 gives each momentum set it
passes to spinoru an id (/nagygen/), and MCFM's nagyqqQQg (marked block) reuses its amplitudes for the same id and
slot indices: the four-quark structures (identical quarks and the three charge assignments solved for e_q², e_Q²,
e_q e_Q) differ only in the charges, which multiply the amplitudes afterwards. This is the "one evaluation for the
three charge pieces" of the profile page, done by exact reuse instead of rewriting the MCFM routine. First build
failed on a name clash (cache array `ca` against MCFM's colour factor CA in constants.f; renamed). Bit-identical
(91/91, timing outputs identical, harnesses unchanged). r: photon 12.1 → 9.1 s, Z 27.3 → 13.8 s, W 16.2 → 10.6 s
(J8 window, 40k points). **Since the clean-up r is 2.20× (photon), 2.98× (Z), 2.56× (W) faster.** New profile
(photon): dipole Borns incl. the cache 48%, psmc 19% (its quad-precision parts 15%), me41 15%, spinoru 13%,
accept 8%: per-point gains are now small; the main lever is the variance (tails of r, rcorr, d2).

## 2026-10-17 (afternoon) — REPORT-2026-10-17-BUDGET read: error budget of the ZEUS-wide sets; tail diagnosis started

Cluster 1742d80/89c1149 (analysis only; old code before changes 1–3; kp bin11). Share of Σ sd·√CPU (total row): r
72–94% at N(1e-4) (W 90, Z 94, B 94, E 93, G6P 92%), 73–80% at N(1e-5); rcorr + d2 2–5% at 1e-4 but 16–25% at
1e-5 and 17–33% of N(1e-4) − N(1e-5); vi, kp, b2 2–5%; Ex7's r-x 24% at 1e-4 in E P2B. Tails: rcorr sd/(1.48 MAD)
3.5–8.5 with 0.25–1% of the seeds beyond 20 MAD carrying 34–94% of its variance; d2 2.0–4.7; r 1.1–2.1 (B: 0.11% of
the r seeds carry 49% of its variance at 1e-4; Z 16%, E 27% at 1e-5; W, G6P none); r-x 7.3. Core-h for an NNLO error
of 10% of NNLOJET's value: W 1.1e4 / 5.2e4, Z 1.8e4 / 6.1e4, B 2.2e4 / 6.7e4, E 2.8e4 / 1.1e5 at τ_cut 1e-4 / 1e-5
(old code; changes 1–3 divide by ≈ 2–3). The seed counts used were not optimal: the optimal allocation reaches the
same 1e-5 errors with 1.5–3.6× less CPU. Consequence: the tail diagnosis (instrumented scratch copy of nlo31 in
work/disorder-comparisons/nnlo21/tails, not committed) covers r and rcorr in the J8 window.

## 2026-10-17 (evening) — tail diagnosis of r and rcorr: three mechanisms, no bug

Instrumented scratch copy of nlo31 (work/disorder-comparisons/nnlo21/tails, not committed): every event's contribution c
to the total-row cells at τ_cut 1e-4 and 1e-5 (J8 window, production settings), and for large events every flavour
structure and dipole (value, τ̃₂ of the mapped Born, acceptance in both cells), the invariants and momenta. r: 6
seeds × 200k × 6 (3.3M events, it ≥ 2); rcorr: 4 seeds. Page: https://claude.ai/artifact/KyFDhs5pSGZNhvVcwN1iH6
- Concentration: in r 30 events carry 62–72% of Σc², the top 0.1% 92–96%; in rcorr at 1e-5 the top 5 of 573k events
  carry all of it (the largest 70× the whole r cell in its iteration).
- (1) Large sampling weights: moderately resolved r events (τ₂ ~ 1e-3, smallest invariant ~ 1e-3 W²) with weights
  10³–10⁵ above typical (VEGAS grid concentrated on the singular regions of the target cell).
- (2) Real and dipoles on different sides of τ_cut: the τ cut acts on the real and on each mapped Born separately;
  a dipole whose τ̃₂ straddles τ_cut is left unsubtracted in that cell (≈ 1/s_ij). Reals with τ₂ > 1e-3 still carry
  25% of the variance of N(1e-4) − N(1e-5) this way. Variance of the 1e-5 cell by the real's τ₂: > 1e-3 50%,
  1e-4–1e-3 19%, 3e-5–1e-4 16%, 1e-5–3e-5 15%.
- (3) Double-unresolved reals just above τ_cut (rcorr, d2; r at 1e-5): the largest rcorr event is q + two hard jets
  (18, 23 GeV) + two soft gluons of 0.45 and 2.8 MeV, nearly collinear to each other (s₄₅ = 7e-10 W²), τ₂ = 1.14e-5;
  removing the softer gluon leaves τ̃₂ = 1.0–1.14e-5 (straddling), removing the harder 1.2e-6; |M|² = 3e16 left largely
  uncancelled. Expected behaviour of slicing above the cut with NLO dipoles, not a bug (the singly unresolved
  cancellations work; τ₂ behaves collinear-safely).
Options (page): A defensive sampling (VEGAS mixing / damping, psmc flat weight); B psmc channels for double-soft /
soft-collinear pairs near τ_cut; C restricted dipoles (α < 1) with the α terms in I, K + P; D smooth τ_cut (average
over τ_cut columns); E pointwise LP subtraction above the cut. Suggested: D and A first, then B, then C.

## 2026-10-17 (night) — variance options D and A tested; AK: fibred P2B (fiber-p2b) as the systematic route

AK: D and A first; keep the parallel fibred-P2B project in mind (Born-fibred subtraction, in principle P2B from N3LO
down to NLO; NLO promising there, NNLO gain still to be shown).
**D (average over τ_cut columns): rejected.** From the J8 per-seed data: the average of columns is never better than
its best column (e.g. avg(3e-5, 1e-5) K = 1.8e-12 vs N(3e-5) 1.2e-12, N(1e-5) 3.5e-12 core-h·(pb/GeV²)²): the large
fluctuations are common to all columns below an event's τ₂, so averaging only acts like a larger τ_cut.
**A (defensive VEGAS mixture, scratch copy: a fraction ε of the points uniform, all weighted by the mixture
density, Jacobian ≤ 1/ε; ε = 0 unchanged).** 200k × 6 per seed (16 seeds): spread down 2–6× and tails gone, but that
size is pathological without the mixture. **Production size (1M × 6; J8 32 seeds, mode 3 16 seeds; ≈ 120 core-h on
the thservs):** means unchanged everywhere. J8 r total: gain in variance per run time 1.1–2.0 (ε = 0.05, 0.2), 0.6 at
3e-5 for ε = 0.2 (within the noise of 32 seeds). J8 rcorr: 3.6 at 1e-4, 0.4–0.8 below (mechanism 3, not helped).
Mode 3 (N3LO P2B, default VEGAS target = the ≥ 1 jet row): ≥ 1 jet 1.0, but the other bins are dominated by single
outlier seeds at ε = 0 (1.48 MAD/sd 10–38%) and near-Gaussian with ε = 0.2: gains 15–250 (p_T 8–15), 4.5–190 (p_T
30–100), 800–4300 (≥ 2 jets), two forward-y bins 0.2 at 1e-5 (outliers of their own); order-of-magnitude only with
16 seeds. Next (overnight): mode 3 with VTARGET = 0 (12 Oct: 3–14× in jet bins) with and without ε = 0.2, 32 seeds.
**Overnight (18 Oct, 01:15): mode 3 (N3LO P2B) r, 4 variants × 32 seeds at production size (1M × 6; 96 + 32 jobs).**
Gain in variance per CPU time against the default target without mixture, per lab bin, τ_cut 1e-4 / 1e-5:
default + ε 0.2: ≥ 1 jet 2.0 / 1.3, p_T bins 2–115 / 3–134, y bins 2–13 / 2–20, ≥ 2 jets 159 / 202;
VTARGET = 0: 2.3 / 1.3, 6–52 / 3–182, 7–30 / 3–53, 600 / 339;
**VTARGET = 0 + ε 0.2: 8.7 / 4.0, 43–216 / 30–154, 16–87 / 23–235, 2171 / 1823** (1.5–15× over VTARGET = 0 alone),
the only variant near-Gaussian in every bin (1.48 MAD/sd 54–133%; the default without mixture: 9–54%). Means agree
(except the outlier-hit default in ≥ 2 jets). **Correction:** the N3LO 1+1 P2B runs of 12 Oct (runs/n3p) used the
default target without mixture: their jet-bin errors (0.5–6% of LO) were inflated by single outlier seeds, and with
such tails the error estimates themselves were unreliable; the 12 Oct statement that the jet bins are flat within
errors stands only at that (poor) precision. Recommendation for mode 3: VTARGET = 0 and VEGAS_MIX = 0.2; for the J8 /
ZEUS r: VEGAS_MIX = 0.2 (1–2×, harmless). nlo31 now has VEGAS_MIX (off by default, bit-identical).

## 2026-10-18 — restricted dipoles (option C): implemented and validated by alpha-independence

AK (18 Oct): push VEGAS_MIX (061d47c) and start on the restricted dipoles. DIPALPHA = α (nlo31; dip41_alpha,
iop31_alpha; default 1 = the previous code, bit-identical 91/91): FF dipoles only for y < α, FI for 1 − x < α, IF
for u < α (Nagy–Trócsányi). The removed regions are away from every singularity, so the change of the integrated
dipoles is a finite 4d integral; with colour conservation it needs only the Born factors already there:
- I (FF): + Σ_{I,K final} ⟨T_I·T_K⟩ [ln²α − γ_I/T_I² (α − 1 − ln α)];
- K + P, regular (iop31_akernel): FI (diagonal) θ(x < 1 − α) [−2 T_b² b ln((2−x)/(1−x)) − g]/(1−x); IF −b ∫_α^1
  du/u ⟨V^{ab}(x,u)⟩ with the azimuthally averaged kernels (q→q: C_F[2/(1−x+u) − (1+x)], g→q: T_R(1 − 2x(1−x)),
  q→g: C_F(1+(1−x)²)/x, g→g: 2C_A[1/(1−x+u) + (1−x)/x − 1 + x(1−x)]), in closed form.
Checks: the closed forms against direct quadrature of the averaged kernels to 2e-18 (alpha/kernels.py). **α-
independence** (J8 window, production settings, α = 1 and 0.1; r 32, vi 8, kp 8 seeds each; work/.../nnlo21/alpha):
vi + kp + r agrees in all 33 cell × τ_cut comparisons (11 rows × 1e-3, 1e-4, 1e-5), pulls −1.8 … +1.5, while the parts
shift individually (total row at 1e-5: vi −1.03e-5, kp −5.3e-6, r +1.55e-5; sum −6.7e-8 ± 1.9e-7): the α terms are
right to ≈ 1% of their size, in both incoming channels (their mix varies across the bins). Gain for r at α = 0.1:
39.5 against 63.5 min per job; per-seed sd 4.05e-7 vs 3.15e-7 (1e-4), 6.63e-7 vs 8.12e-7 (1e-5): variance per CPU
1.0× at 1e-4, 2.4× at 1e-5 (±35% from 32 seeds), mostly from speed.

## 2026-10-18 (afternoon) — N3LO jets in the set-up of NNLOJET's N3LO P2B paper (1803.09973); scale variations

AK (18 Oct): push and tag (nnlo21-2026-10-18 at 418666a, done); then plots of the N3LO P2B lab-frame jet p_T and
rapidity against the NNLOJET paper; do we have scale variations; redo LO, NLO, NNLO. Work in
`work/disorder-comparisons/zeus05` (z5lib.py, z5plot.py, z5scale.py, z5size.py, digit/), page
https://claude.ai/artifact/Sa7MEGe3N2f5p1SudqcN5x.
- **Set-up.** The paper compares with ZEUS hep-ex/0502029 'global': E_p = 820, E_e = 27.5 GeV, Q² > 25 GeV²,
  y > 0.04, E_e' > 10 GeV, longitudinally invariant k_t, R = 1, E_T recombination, jets E_T > 6 GeV, −1 < η < 3
  (every jet counted; Q² and x of the events with a jet), NNPDF31_nnlo_as_0118 at all orders, μ = Q, 7 points,
  photon exchange. This is `analysis/nnlojet.f` of the disorder paper and NNLOJET's `DIS1j.run.template`
  (AK's `disorder/paper_runs/nnlojet/nnlojet_data/DIS1j`, with NNLOJET's LO, NLO, NNLO at 7 scales and the ZEUS
  tables). Our lab11 set-up (E_p = 920, 125 < Q² < 20000, anti-k_t, leading jet) is not this one, so the lab11
  N3LO results of 12 Oct cannot be compared with the paper. New: `LABSET=zeus05` (mode 3) in nlo31 and sliced21
  (d62836d): beams, PDF, cuts, k_t E_T-scheme jets (as FastJet's Et_scheme), 31 rows (≥ 1 jet, η_inc 8, E_T,inc 8,
  Q² 7, x 7); beam grid Q = 5 … √s, x ≥ 2.5e-4 (44 nodes). Default bit-identical (91/91; only the LCELLS header
  text changed).
- **The paper's N3LO.** Only in its figures (arXiv source, gnuplot/cairo vector PDFs). Converted to SVG, axes
  calibrated from the ticks, every step polyline read (digit/extract.py; gnuplot merges equal neighbouring bins,
  so segments are mapped to bins by position). Check: LO, NLO, NNLO read this way = NNLOJET's files to ≤ 0.001 in
  the ratio panels, ≤ 0.06 % in the linear η panel, ≤ 0.3 % on log axes. N3LO = ratio panel × NNLOJET NNLO
  (±0.05 % of NNLO). N3LO/NNLO: η 0.98, 1.02, 1.03, 1.04, 1.09, 1.14, 1.19, 1.25; E_T 1.05 … 1.003. The paper
  plots ZEUS σ divided by the hadronisation factor.
- **α_em:** disorder and nlo31/sliced21 use 1/137, NNLOJET 1/137.036: disorder LO with 1/137.036 = NNLOJET to
  −0.003 ± 0.004 % (total), every bin within errors; comparisons rescale ours by (137/137.036)².
- **LO, NLO:** disorder (-p2b) and ours (inclusive + sliced21 b0) agree with NNLOJET's NLO coefficient within
  ±0.25 % of NNLO in every bin (pulls mostly < 2).
- **NNLO, ours** (inclusive -nnlocoef + b1 + lo, mode 3, 16 + 16 seeds): agrees with NNLOJET's NNLO coefficient
  within ≈ 1.5σ in every η and E_T bin for τ_cut ≤ 2e-4 (errors 0.2–3 % of NNLO); above 1e-3 the known drift
  (total −0.48 ± 0.11 % at 1e-3). disorder NNLO: disorder -p2b -nnlocoef (32 seeds × 1.6e8 calls, 3 h each) = NNLOJET within ≤ 1.9σ in
  every bin (total −0.32 ± 0.29 % of NNLO).
- **Convergence against plain DISENT** (AK, page comment): CPU of our chain (b1 + lo at 1e-4, optimal split) over
  disorder -p2b for equal error on the NNLO coefficient: total 0.1; η bins 0.6, 0.3, 0.3, 0.3, 1.1, 9.8, 47, 31;
  E_T bins 0.1, 0.4, 2.5, 2.0, 5.6, 1.4, 21, 42 (±40 %, 16 seeds per part): ours is faster in the central bins,
  10–47× slower forward and at high E_T (the pattern of 6 Oct). The automatic page reply had guessed "slower
  everywhere"; corrected in the thread.
- **The inclusive N3LO coefficient** (disorder, 7 scales) is small here: −23 … −28 pb per unit η centrally, 0 in
  the forward bins, against the paper's +90 … +390: the N3LO comes almost entirely from P2B on the NNLO 2+1.
- **N3LO pilot** (112 jobs: b2, vi, kp, rc1, rc2 16 each, r 32; ≈ 100 CPU-h): r has a per-seed sd of 0.9–1.6e4 pb
  per bin in the central η bins at τ_cut 1e-4 (1-h jobs), Gaussian (sd ≈ 1.48 MAD), > 90 % of the variance, ×2 at
  1e-5. Optimal allocation over the parts: at 50k CPU-h the forward η bins (1–3) at 6–13 % of the paper's N3LO
  correction, η 0.5–1 at 35 %, E_T 6–21 GeV at 17–25 %; η < 0.5 not resolved (5–9 M CPU-h for 10 %).
  **Correction:** this afternoon I told AK the production would likely need 5–20k CPU-h; that was wrong by orders
  of magnitude for the central and backward bins (the mis-binning variance of the R − D cancellation, mechanism of
  17 Oct, is far larger at Q² down to 25 GeV² than in the lab11 set-up).
- **Scale variations** (AK: go; 1b16531). `SCALES=N` (mode 3): the O(a^N) cells at the 7 (ξ_R, ξ_F) of
  disorder/NNLOJET by the renormalisation group, a = α_s/2π, nf = 5:
  C₃⊗f = c₃f + 2b₀L_R c₂f − L_F c₂⊗P₀f + (b₀²L_R² + b₁L_R) c₁f − L_F c₁⊗P₁f + ½L_F² c₁⊗P₀P₀f +
  b₀(½L_F² − 2L_F L_R) c₁⊗P₀f (and C₂⊗f = c₂f + b₀L_R c₁f − L_F c₁⊗P₀f), all with f at ξ_F Q and a at ξ_R Q
  (LHAPDF), the logs of the c_n at Q. Every event is re-summed with the PDF variants (flavour sums split from the
  matrix elements: born_part fsum, real_part fsum_r, b21_eval bval); P⊗f from HOPPET (hoppetEvalSplit, iloop 1,
  2, 11; normalisation checked: df/dln μ² = aP₀⊗f + a²P₁⊗f up to NNPDF's a³P₂, 0.1–0.7 % at Q = 50, 2–5 % at
  6 GeV); b1, b2 use beam-grid slots of f and P₀⊗f at 2Q and Q/2 (lp21 fdist per slot, mktab LPGVAR). No
  internal procedures passed as arguments (they need an executable stack). Checks: default bit-identical (91/91);
  SCELL(1,1) = LCELL exactly for n = N (r, kp, lo, b0); NLO at 7 scales (inclusive + b0) = NNLOJET to ≤ 0.1 % of
  NNLO in every bin and scale; at (½,½) the total is −0.08 % (11σ of tiny errors) from disorder's inclusive part.
  NNLO at 7 scales (b0 with P₀⊗f, b1, lo): the plain per-scale comparison has pulls ≤ 2 everywhere, but in
  η 0–0.5 and 0.5–1 its errors are 27 % of NNLO: one b1 seed (4206) has −1.4e4 / +7.5e3 pb in the two bins
  (78 and 274 MAD, anti-correlated: one heavily weighted near-soft 2+1 Born across the bin edge; the same code path
  as the plain runs, a genuine heavy tail). The scale dependence, C(k) − C(central) per seed, cancels it: ours
  against NNLOJET to 0.1–0.3 % of NNLO, pulls ≤ 1.5 in nearly every bin and scale; the total row 7.06 vs 7.30 %
  at (½,½) (−3.2σ) and −3.99 vs −4.11 % at (2,2) (+2.1σ), most likely disorder's α_s in the inclusive NNLO
  coefficient (its NLO part showed −0.08 % at (½,½)). Passed, with that caveat.
- **disorder and NNLOJET at μ_R = Q/2:** disorder's own P2B NLO coefficient differs from NNLOJET by up to −0.5 % of
  NNLO at μ_R = Q/2 (forward bins, all DISENT 2+1): disorder runs α_s with nf = 5 from α_s(M_Z), NNLOJET (and our
  chain) uses LHAPDF's α_s (nf = 4 below m_b; Q/2 reaches 2.5 GeV here). A convention difference, not a bug.
- AK (evening): run the variance study of r tonight, push once the NNLO scale check passes, and add a plot of the
  two jet histograms with ours on NNLOJET and the ratio per order (page v2). Variance study (runs/var, 32 r seeds
  each against the pilot's 32): B DIPALPHA = 0.1, C VTCOL = 8 (new, d1eba04: VEGAS target on the τ_cut 1e-4
  column), D both.
- AK (page comment): 10× the NNLO statistics: 144 more seeds of b1 and lo (SCALES = 2, runs/sc2), on the thservs:
  `-m all` was refused by the permission classifier because the OK came through the page comment. AK then gave the
  rule in the chat: `-m all` is allowed from Friday 20:00 to Sunday 02:00, otherwise thservs only unless he says so
  (memory, ~/COMPUTE-JOBS.md).
- Slip: the four smoke tests of scale mode (minutes) ran in the session shell at nice 7; everything else via
  ssh/thsub at nice 0.
- **Variance study of r (evening, runs/var, 32 seeds each against the pilot's 32; zeus05/z5var.py):** gain =
  variance × CPU of the pilot / the variant, τ_cut 1e-4. B (α = 0.1): 2142 against 3394 s per job; η bins 0.94,
  0.12, 2.35, 1.47, 2.18, 1.43, 7.69, 2.79, ≥ 1 jet 0.08, lowest E_T 0.11: faster and better in some forward bins,
  but new heavy tails (sd/1.48 MAD 4–8 in ≥ 1 jet, η −0.5–0, E_T 6–8). C (VTCOL = 8): no gain (0.3–2.6), means = A
  (pulls ≤ 1.1). D (both): as B. No robust gain, nothing near the ~100× the central η bins need: the sizing stands
  (50k CPU-h resolves η 1–3 and E_T 6–21 GeV, not η < 0.5). The central bins need a binning-safe R − D (fibred P2B),
  not tuning.
- **Z05 cluster production** (AK: "let's do the 50k run now"): instructions UPDATE 18 Oct "Z05" pushed (9670206):
  bin12 (needs HOPPET), grids z05 + 4 variants, SCALES = 3 everywhere; r 18,700 × (2M × 6), vi 2,729, kp 2,479,
  b2 729, rc1/rc2 1,000 each, RG terms b0 200, b1 400, lo 400: ≈ 47k core-h, 27.6k jobs; optimal allocation for
  η 1–3 and E_T 6–21 GeV (expected 6–8 % of the paper's N3LO correction for η 1.5–3, 16 % for η 1–1.5, 19–27 % for E_T
  6–21 GeV). First wave 5 %: CPU projection (cut r if > 50k) and means against the MPP pilot. nnlo21/cluster/
  summary_z05.py writes the per-part cell statistics for the combination at MPP.
- **NNLOJET's CPU for the paper (AK asked; estimate).** Evaluations in AK's files (antenna, 1 jet): LO 5.0e10,
  NLO 1.1e11, NNLO 6.7e10. NNLOJET v1.0.2 per event on thA371a (one warmup iteration, zeus05-like cuts;
  zeus05/njtime): epLJ LO 2 µs, RR 0.49 ms, RV 0.12 ms, VV 0.09 ms; epLJJ LO 2.5 µs, R 0.16 ms, V 12 µs, RR 52 ms
  (first iteration; ≈ 9 ms in our cluster runs), RV 3.7 ms, VV 0.08 ms. The LO–NNLO 1-jet curves: ≈ 5k core-h
  (2–9k by the RR/RV/VV split). The N3LO (NNLOJET's dijet NNLO with P2B) is not recorded: order 10⁵ core-h,
  plausibly a few × 10⁵ (6.7e10 RR-dominated evaluations at ≈ 10 ms ≈ 2e5; or our cluster's NNLOJET ZEUS-dijet
  NNLO, 12.8k + 2.7k core-h for ±1 % of NNLO on the total, scaled to ≈ 1 % per bin in a harder region).
- **disorder as before** (AK: verify after the NC/MATFOR changes): today's disorder (16/16/32 seeds; photon) against
  the 2023 disorder runs of the paper (paper_runs/nnlojet, central scale): LO, NLO, NNLO agree in every η and E_T bin
  (largest pull 2.7: NLO η 2.5–3, −0.13 %; NNLO ≤ 2.2), and with NNLOJET (NLO ≤ 0.25 % of NNLO, NNLO ≤ 1.9σ).
- **N3LO error against budget** (pilot, optimal split per bin), in % of the NNLO cross section per bin: 5e4 CPU-h:
  central η 2.7–3.0, forward 0.9–1.4, E_T 0.5–1.5; 3e5: 1.1–1.2, 0.4–0.6, 0.2–0.6; 1e6: 0.6–0.7, 0.2–0.3,
  0.1–0.3. 1 % of NNLO in every bin ≈ 4.4e5 CPU-h, 0.3 % ≈ 5e6 (NNLOJET's NNLO 1-jet errors are 0.03–0.05 % of NNLO
  per bin; their N3LO errors are not published).
- **10× NNLO statistics** (AK's page comment; 282 of 288 jobs at 22:00, the rest on an overloaded thserv17 where
  the autonice daemon reniced them to 13 after a few hours; thsub started them at 0): b1 170, lo 176 seeds. Ours
  − NNLOJET at τ_cut 1e-4: total +0.07 ± 0.06 % of NNLO; other bins within 0.2–0.7 % errors, flat from 2e-4 to
  1e-5. η 0–0.5 / 0.5–1: −2.4 / +2.5 ± 2.6 %, all from b1 seed 4206 (one near-soft 2+1 Born across the η = 0.5
  edge; the binned P2B difference has a log-divergent variance in the soft limit of the 2+1 Born, weight ∝ 1/E,
  edge-crossing probability ∝ E); without it (diagnostic only) 0.1 ± 0.3 %. Per scale (figure on the page):
  ours/NNLOJET at each of the 7 scales within errors at LO, NLO and NNLO.
- **NNLOJET's N3LO MC error from the scatter of the paper's points** (AK): polynomial fits (degree 4–5) to N3LO −
  NNLO in η and ln E_T leave residuals of 0.3–0.4 % of NNLO; the same test on NNLOJET's NNLO coefficient (true MC
  error 0.03 %) leaves 0.3–0.9 %: shape, not noise, so ≤ 0.3–0.4 % per bin, plausibly ≤ 0.1 %. The Q² inset (read to
  0.01 %): N3LO/NNLO = 1.0005, 1.0000, 1.0001, 1.0000 above 250 GeV², smooth to ≤ 0.02 %, but the P2B part hardly
  enters there. For us (pilot): 0.3 % in the central η bins ≈ 5e6 CPU-h, 0.1 % ≈ 4e7.

## 2026-10-19 — Z05 interim (cluster, r at 35 %): first N3LO against the paper; correction of the sizing

REPORT-2026-10-19-Z05-interim (cluster branch 6b26d4a; summary file results/Z05-summary-interim.json.gz, combined
at MPP with zeus05/z5n3lo.py and disorder's inclusive -n3locoef): all parts done except r (6545 of 18700; vi
2729, kp 2479, b2 728, rc1 1000, rc2 994, b0 200, b1 400, lo 400); 15.8k core-h used, 32.6k projected (cluster
cores ≈ 1.4× faster than the thservs); means within 1.4σ of the MPP pilot; no failed runs.
- δN3LO (τ_cut 1e-4, central scale, % of NNLO) against the paper: forward η agrees (η 2–2.5 18.3 ± 2.1 vs 19.3,
  2.5–3 23.4 ± 2.1 vs 24.9; 1.5–2 8.7 ± 2.4 vs 14.1, −2.3σ); central η unresolved (errors 4–11 %); E_T ≥ 21 GeV
  ours −1 … −4 % vs the paper's +0.3 … +1.7 % (pulls −2.0, −3.3, −1.5, −0.9). There δN3LO is a 0.4 % residual of
  b2 ≈ +970 … +1170 % of NNLO, r ≈ −730 … −920 %, vi ≈ −300 … −340 %, kp ≈ +60 … +90 %: a 0.3 % effect in b2 or r
  (τ_cut power corrections, largest in hard bins) would explain it; 3e-5/1e-5 too noisy yet. Open.
- **Correction of 18 Oct (sizing):** the cluster's r (2M × 6) has a per-seed sd 1.3–1.7× the pilot's (1M × 6) in
  the central η and low-E_T rows, where 0.71× was expected; sd/1.48 MAD 2.1–2.8 there (6545 seeds) against ≈ 1–1.5
  in the pilot's 32: a heavy tail the pilot did not sample (the lesson of 16 Oct, size from ≥ 1000 seeds, not
  applied). Variance ≈ 5× higher in those rows; forward rows as projected (0.58×). Corrected budget (cluster sd and
  CPU per job): 1 % of NNLO in every bin ≈ 1.4e6 CPU-h (not 4.4e5), 0.3 % ≈ 1.5e7 (not 5e6).
- **Correction of the Z05 instruction:** I asked for bin12 = bin11 on a J8 r seed; the bit-identity holds against
  the tag 418666a, which already contains d5c2021 (|z|², last-bit changes amplified by VEGAS; validated 17 Oct at one
  iteration). The cluster's finding (iteration 1 identical, later iterations differ) is that known change.
- **Error budget per part (Z05 interim, τ_cut 1e-4, central scale):** r carries 92–99 % of the variance (82–97 % at
  full r) and 84–96 % of the cost in every bin below E_T = 29 GeV; above, vi and b2 up to 18 % and 13 %; kp
  negligible. Cost for 1 % of NNLO per bin against τ_cut (optimal split, cluster sd and CPU): worst bin 3e4 (2e-2),
  6e4 (5e-3), 2e5 (1e-3), 1.4e6 (1e-4); lower τ_cut is worse (× 2–4 at 3e-5, × 3–10 at 1e-5 without outliers; η −1…0
  at 3e-5 outlier-dominated, 1e7–1e8). Larger τ_cut needs the power corrections removed.
- **P2B-improved slicing was not used at N3LO** (AK thought it was): P2BSLICE was read only in modes 1–2 (the NNLO
  DIS+jet validations); mode 3 (N3LO P2B) used plain τ₂ slicing of the 2+1 (form (A) of the 15 Oct nested-P2B
  check). Now wired into mode 3 (1d4d8b7, default bit-identical 91/91): below τ_cut O(event) − O(P2 event) (the
  1+1 Born terms cancel since project21 keeps q), the r drop/degenerate handling of P2BSLICE for every mode.
  Smoke test: P2BEXTRA lo term ≈ 700 pb at 2e-2, ≈ 0.1 pb at 1e-3; r flags 13 of 15815 events for plain weights.
  Tests submitted (thsub, 64 jobs): (A) lo with P2BSLICE (32 seeds) + the 170 b1 seeds: NNLO 1+1 flatness in τ_cut
  against NNLOJET; (B) r with P2BSLICE (32 seeds, 2M × 6): per-seed sd per τ_cut column against the cluster's plain
  r. Queued behind ≈ 7,900 proVBFH-cs jobs of the VBF session (thsub has no priority).
- **P2B-improved slicing in mode 3, first tests** (19 Oct, desktops via thsub -m all with AK's OK; extras only,
  P2BEXTRA = 1, VTARGET = −1, added to the plain runs; zeus05/z5p2bx.py): (A) one order lower, lo extras (32) +
  the plain NNLO: against NNLOJET (% of NNLO) the total at 2e-2 −2.19 → −0.70 ± 0.21, at 2e-3 −0.67 → −0.39, at 1e-3
  −0.35 → −0.30; η 2.5–3 at 2e-2 −7.4 → +0.3 ± 0.4: the fiducial drift is removed, flat to ≈ 2e-3 in most bins
  (as CNV 2408.05265); the remainder is hadronic (total −0.7 % at 2e-2, η 1–1.5 +3.3 % at 2e-2, +1.3 % at 5e-3). lo
  extras have no tails (sd ≈ 1.48 MAD). (B) N3LO, r extras (32, 2M × 6): per-seed sd ≈ 1–3e4 pb per central η bin
  at every τ_cut (2e-2, 1e-3, 1e-4), 15–25× plain r at 2e-2, sd/1.48 MAD 3–7: the extras reach τ₂ → 0, where the
  double real is doubly unresolved relative to the 2+1 Born (NLO 3+1 dipoles do not cancel it; O − O∘P2 damps it
  only linearly; tail mechanism 3 of 17 Oct). P2B fixes the fiducial power corrections but does not reduce the cost
  at N3LO: plain + extras at 1e-3 is noisier than plain at 1e-4. vi and kp extras are small. AK asked whether we
  include the hadronic (NLP) power corrections of N-jettiness: no, leading power only (Born-frame measure, P2B for
  the fiducial ones); NLP is known for 0-jettiness and 1-jettiness V/H+jet at NLO, not for our τ₂ and not at
  O(α_s²) relative for any jet process.
- AK asked to confirm that the mappings keep the lepton kinematics: yes (checked in the code). dip41 copies q and
  both leptons unchanged (born_slots), FF moves only final-state partons, FI/IF rescale the incoming parton and
  conserve q (no II dipoles in DIS); project21 keeps q; the 1+1 Born is q + xP at the event's x, Q²; sliced21
  evaluates the 2+1 Born at the same leptons; the lab frame is built from the leptons and the proton direction.
  So events and counter-events share x, Q², y, the y and E_e' cuts and the x/Q² bins; the Q²/x rows can still differ
  through the jet requirement, and the lab jet bins can differ. Next, when the full Z05 is in: a fit of the τ_cut
  dependence for stability (AK), and the per-bin budget from the final summary.

## 2026-10-19 (evening) — Z05 final: N3LO against the paper; τ_cut fit; budget

REPORT-2026-10-19-Z05 (cluster branch c87ee91; results/Z05-summary.json.gz): every part complete (r 18700, vi
2729, kp 2479, b2 729, rc1 1000, b0 200, b1 400, lo 400) except 6 slow rc2 seeds (994/1000; rc2 enters only at
3e-5, 1e-5); 32.3k core-h; 24 r files written twice by zombie et tasks quarantined and rerun. Combined at MPP with
zeus05/z5n3lo.py (figures now show τ_cut 1e-4 and 1e-3), fit zeus05/z5fit.py, budget zeus05/z5budget.py; page
v13. δN3LO in % of NNLO, central scale, against the paper:
- **Forward η (1.5–3): ours is below the paper.** At 1e-3 by 2.3–3.7 % (pulls −5.3, −3.5, −3.2), at 1e-4 by 2.4–4.6 %
  (−3.3, −2.7, −1.8). Flat in τ_cut from 2e-3 to 1e-4 (η 2.5–3: 23.3, 22.6, 23.1, 23.6, 22.5 against 24.9). The
  inclusive N3LO term is < 0.01 % of NNLO there, and in η 2.5–3 LO is 3.5 % of NNLO, so the correction is essentially
  our NNLO 2+1 for forward lab jets. The read-off reproduces NNLOJET's LO/NLO ratios in η to 0.07 %; NNLOJET's N3LO
  error ≤ 0.3–0.4 % (18 Oct). Open: form A of the fit (below) puts c0 on the paper in η 2–3 (17.5 ± 1.5, 25.3 ± 1.5),
  A2 does not (17.0 ± 0.8, 20.8 ± 0.9), so a residual τ_cut dependence below 1e-4 that is flat above it is not
  excluded by these data.
- E_T 21–71 GeV below the paper at both τ_cut (1e-3: −0.5 ± 1.0 vs +1.7, −1.4 ± 0.6 vs +0.9, +0.3 ± 0.5 vs +0.7;
  1e-4: −5.4 ± 2.0, −3.4 ± 1.0, −1.2 ± 1.1); the 1e-3 → 1e-4 drop is 1.3–2.7σ depending on the column correlation
  (not in the summary). Central η: errors 4–9 % at 1e-4; at 1e-3 2–3 %, agreeing with the paper within 0.4–1.8σ.
- **Correction of the interim (this morning):** "forward η agrees" held only within its errors (18.3 ± 2.1,
  23.4 ± 2.1); with full r ours is 2–5σ below. The projected full-r errors were too small in central η: the per-seed
  sd of r at 1e-4 grew 1.4–1.5× (η −1 … 0.5) and 1.1–1.3× (η 0.5 … 1.5) between 6545 and 18700 seeds with the MAD
  unchanged (tail variance not converged); errors came out 9.4, 6.9, 5.8 against 6.9, 4.6, 4.0 projected, and are
  probably still underestimates. Budget at 1e-4 (final sd, cluster CPU per job, optimal split): 1 % of NNLO in every
  bin ≈ 2.5e6 CPU-h (interim 1.4e6), 0.3 % ≈ 2.8e7; at 1e-3 the worst bin needs 2.6e5.
- **τ_cut fit (AK): it does not help.** Per bin, c0 + τ(c1L + c2L² + c3L³) (A) and with c3 = 0 (A2), over 2e-2…1e-4,
  2e-2…1e-3, 5e-3…1e-4. The summary has no covariance between the τ_cut columns (each seed fills all), so two models:
  independent columns, and nested sums with independent τ bands (cov = min(V_c, V_c') for vi, kp, r; b2 fully
  correlated); 128 local r seeds give correlations between the two (tails decorrelate neighbouring columns), so χ²
  is indicative only. A: c0 error no smaller than the 1e-4 column's (central η ±4–8 vs ±4–9; forward ±1.5–1.6 vs
  ±1.3–1.4). A2: errors near the 1e-3 column's (1.5–3× smaller than 1e-4), but c0 moves with the range by up to
  2–3σ and χ² is poor (η 1.5–2 21/5, η 2.5–3 20/5, E_T 14–21 27/5, 29–47 54/5, 47–71 35/5, correlated model). Over
  the 12 fits c0 spans 33 % of NNLO (η −1 … −0.5), 22 % (−0.5 … 0), 5–7 % forward. With √τ terms (fiducial) c0 is
  unconstrained (errors 3–27 %). The large-τ columns are precise but carry power corrections up to 19 % of NNLO; c0
  then rests on the small-τ columns, as noisy as the 1e-4 result: model dependence replaces statistics. The limit is
  the r variance as τ₂ → 0 (double-unresolved), as AK said. The column covariance per part could be written by the
  cluster from the existing .out files (cheap); it would make χ² meaningful and decide the high-E_T drop, not the
  conclusion on the fit. Not requested.
- AK (19 Oct evening): "this points to something being simply wrong"; first the column covariance from the
  cluster (UPDATE 19 Oct evening in nnlo21/validation/CLUSTER-INSTRUCTIONS.md, script nnlo21/cluster/cov_z05.py:
  per part and row the 10 × 10 covariance and Spearman correlation of the central-scale columns; tested on the
  local sc2 b0/b1/lo runs: diagonal = sd² to 2e-16). Scale of the forward offset: at 1e-4 it is −49, −18, −6 pb per
  bin in η 1.5–2, 2–2.5, 2.5–3, i.e. 0.13–0.22 % of b2 (or of |r|), 20–50 % of the NNLO 2+1 coefficient there, and
  6–11 % of the 2+1 Born (NLO − LO). Our earlier NNLO 2+1 validations against NNLOJET (Breit-frame dijets, 10–30 %
  errors on the NNLO coefficient) could not see an error of this size.
