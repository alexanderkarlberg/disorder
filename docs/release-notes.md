# disorder release notes

Changes since version 1.0.0 (SciPost Phys. Codebases 32, arXiv:2401.16964).
The 1.0.0 manual (`docs/disorder-1.0.0-manual.tex`) is left unchanged;
where it no longer applies, this file says so.

## Unreleased

### Testing and validation

- **Unit tests.** `ctest` now runs unit tests of the Born phase-space
  generator, the lab ↔ Breit frame transformations, the matrix element (at LO
  against the quark-parton-model expressions for NC γ/Z/γZ exchange, CC, and
  neutrino beams; beyond LO the combination of structure functions), the
  interface to DISENT, and the refusal of unsupported option combinations.
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
