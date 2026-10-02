# disorder release notes

Changes since version 1.0.0 (SciPost Phys. Codebases 32, arXiv:2401.16964).
The 1.0.0 manual (`docs/disorder-1.0.0-manual.tex`) is left unchanged;
where it no longer applies, this file says so.

## 2.2.1

### Physics

- **DISENT no longer drops the O(αs²) part of a fraction of its events.**
  With 1 − X < CUTOFF for the collinear X of the K and P terms, `VIRTHR`
  ended the event after the three-parton Born had been handed over, so the
  O(αs²) contributions of about CUTOFF^(1/npow2)/2 of the events (0.5% with
  the defaults `-npow2 4`, `-cutoff 1e-8`) were lost; in expectation this
  removed the pieces that do not depend on that X (the virtual and the reals
  with a final-state spectator). Now only the excluded sliver of the x
  integral is dropped. With `-p2b`, O(αs²) distributions change by about half
  a percent of their O(αs²) 2+1 part; total cross sections are unchanged.

## 2.2.0

### Physics

- **p2b for all NC and CC processes, up to NNLO.** DISENT's matrix elements
  (`src/libdisent.f`) now include Z and W exchange: NC with γ, Z or γ/Z
  (`-includeZ`, `-Zonly`, `-intonly`), CC, charged-lepton and neutrino beams
  of either charge (`-positron`, `-neutrino`). `-p2b` can therefore be
  combined with `-includeZ`, `-CC` and NC neutrino beams, giving
  differential predictions (jets, event shapes) at O(αs) and O(αs²) for all
  these processes. The electroweak couplings are the same as in the
  inclusive (structure-function) mode. For photon exchange the results are
  bitwise identical to before.
  - O(αs): the three-parton tree `MATTHR` (derivation in
    `derivations/matthr/`); validated against NNLOJET 1.0.2 and
    POWHEG-BOX-RES for all processes each code supports.
  - O(αs²): the four-parton trees (`MATFOR`, including for W exchange the
    identical-final-state interferences), the spin-correlated dipoles
    (`CONTHR3`, `SUBFOR`) and the finite one-loop part (`VIRTHR`, with the
    Bern–Dixon–Kosower amplitudes ported from MCFM); derivation and
    validation in `derivations/o2/`. Neglected, as in the structure
    functions: contributions proportional to the sum of the quark axial
    couplings (Z coupling to a second quark line or to a quark loop), which
    vanish for complete generations. T-odd one-loop terms (sin φ
    correlations from absorptive parts) are dropped; they integrate to zero
    for reflection-symmetric observables.

### Analyses

- **`analysis/cmp_nnlojet_powheg.f`** (build with
  `-DANALYSIS=cmp_nnlojet_powheg.f`; it has its own anti-kT and needs no
  FastJet): the analysis used for the comparisons with NNLOJET and
  POWHEG-BOX-RES. It fills Q², x, y, the leading lab-frame jet (anti-kT,
  R = 1, p_T > 5 GeV) and the Breit-frame event shapes τ_zE, B_zE and ρ_E,
  with `*_Ec` versions that require an energy E_cur > Q/10 in the current
  hemisphere (for τ_zE and B_zE this makes them infrared safe at O(αs²),
  as NNLOJET's `dis_eventshapes`). The observables are defined in
  `analysis/cmp_obs_core.f` and the cuts in `analysis/cmp_obs_cuts.h`,
  shared with the POWHEG-BOX analysis.

### Testing and validation

- **Unit tests.** `ctest` now runs unit tests of the Born phase-space
  generator, the lab ↔ Breit frame transformations, the matrix element (at LO
  against the quark-parton-model expressions for NC γ/Z/γZ exchange, CC, and
  neutrino beams; beyond LO the combination of structure functions), the
  interface to DISENT, DISENT's three-parton matrix element for every
  process (`test_matthr`), the consistency of the four-parton matrix element
  with the dipole subtraction terms in all single-unresolved limits for every
  process (`test_subtraction`), and the refusal of unsupported option
  combinations.
  They take about a second: `ctest -L unit`.
- **Validation through ctest.** Each configuration of the validation matrix
  (now listed in `validation/configurations.txt`) is a ctest test, compared
  with `validation/ref_runs/` by `validation/run_validation.py`. The
  comparison allows a relative difference of 1e-5 (configurable) instead of
  requiring identical text, so it no longer fails on a different compiler or
  platform, and it ignores the position of library banners in the logs. The
  validation tests are registered when configuring with
  `-DNEEDS_FASTJET=ON -DANALYSIS=exclusive_lab_frame_analysis.f`.
  ```
  cmake -S . -B build -DNEEDS_FASTJET=ON -DANALYSIS=exclusive_lab_frame_analysis.f
  cmake --build build -j
  ctest --test-dir build -j 8               # everything, a few minutes
  ctest --test-dir build -j 8 -LE slow      # without the PDF-uncertainty run
  ```
- A NaN or infinity anywhere in the output of a validation run is a
  failure, also when generating references.
- The validation matrix covers p2b for every NC/CC mode and beam: NC γ,
  γ/Z, Z only and interference only; CC; e±, ν and ν̄; NC+CC; at NLO and
  NNLO; with and without scale variations.
- `validation/validate_or_generate.sh validate|quick|generate` still works and
  now uses ctest. It no longer needs GNU parallel (this replaces the
  corresponding paragraph of section "Validating the code" of the 1.0.0
  manual). Extra cmake flags can be given in `EXTRA_CMAKEFLAGS`.
  `validation/ref_runs_quick/` has been removed (it duplicated files in
  `ref_runs/`).
- **CI** runs the unit tests and the full validation matrix through ctest,
  and keeps its scheduled (weekday) runs from being disabled by GitHub
  after 60 days without repository activity.

### Fixes

- HOPPET's grid now always extends to x = 0.1, also when a very large
  `-xmin` is requested.
- `-scale-choice 4`: the printed header and `-help` text said the central
  scale was Q(1−x)/x. The scale used is, as intended, μ² = Q²(1−x)/x, i.e.
  μ = Q·sqrt((1−x)/x); the texts now say so. Results are unchanged. (The
  1.0.0 manual numbers the scale choices differently from the code: in the
  code 3 is μ² = Q²(1−y) and 4 is μ² = Q²(1−x)/x.)

### Internal

- The DISENT callbacks moved to `src/mod_disent_interface.f90` and the PDF
  and structure-function set-up to `src/mod_pdf_setup.f90`; the physics code
  is built as a CMake object library shared by `disorder` and the tests.
  The results of the validation matrix are bitwise identical to before.
- CMake 3.12 or newer is required.
