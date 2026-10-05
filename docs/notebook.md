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
