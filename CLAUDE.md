# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`disorder` is a Fortran program that computes fixed-order (LO through N3LO) predictions for
neutral-current (NC) and charged-current (CC) deep inelastic scattering (DIS). Physics results
must match published references (SciPost Phys. Codebases 32, arXiv:2401.16964), so numerical
changes to the perturbative calculation are scientifically sensitive — see the Validation
section below before touching anything under `src/`.

## Build

Depends on [Hoppet](https://github.com/hoppet-code/hoppet) (v2.1.0+) and
[LHAPDF6](http://lhapdf.hepforge.org/), with [FastJet](https://fastjet.fr/) optional.

```
mkdir build && cd build
cmake ..
make [-j]
```

Useful cmake flags:
- `-DHOPPET_CONFIG=/path/to/hoppet-config -DLHAPDF_CONFIG=/path/to/lhapdf-config` if the
  `*-config` scripts aren't on `$PATH`.
- `-DNEEDS_FASTJET=ON [-DFASTJET_CONFIG=/path/to/fastjet-config]` to link FastJet.
- `-DANALYSIS=my_analysis.f` to compile a specific analysis from `analysis/` instead of the
  default `analysis/simple_analysis.f` (see Analysis system below).

This produces three executables: `disorder` (main program), `mergedata` and `getpdfuncert`
(auxiliary tools built from `aux/`).

## Running

```
./disorder -pdf MSTW2008nlo68cl -n3lo -Q 12.0 -x 0.01 -scaleuncert
./disorder -help   # full list of command-line options
```

All run parameters are command-line switches parsed in `set_parameters` (`src/mod_parameters.f90`)
— that subroutine is the authoritative reference for every flag (`-lo/-nlo/-nnlo/-n3lo`, `-NC/-CC`,
`-p2b`, `-pdf`, `-scaleuncert`, `-pdfuncert`, `-alphasuncert`, `-neutrino`, `-positron`,
`-includeZ`, etc.), rather than duplicating them here.

A run produces `xsct_[...].dat` (total cross section, MC error, scale/PDF uncertainties) and,
if an analysis is linked in, `disorder_[...].dat` histogram output.

## Validation (the test suite)

There is no unit-test framework; correctness is checked by running a fixed matrix of physics
configurations and diffing the output against committed reference results.

```
cd validation
./validate_or_generate.sh validate   # full matrix: build + run + diff against ref_runs/
./validate_or_generate.sh quick      # smaller/faster subset, diff against ref_runs_quick/
```

- `validate_or_generate.sh generate` re-runs the full matrix and overwrites `ref_runs/` — only
  do this deliberately, when a physics/numerical change is intended and reviewed, since it
  redefines "correct".
- CI (`.github/workflows/cmake-single-platform.yml`, driven by `.github_CI.sh`) builds then runs
  `validate_or_generate.sh validate` (the full matrix), not the quick one.
- All three modes hardcode command-line configurations (`prefix_full`/`cmdline_full` for
  `validate`/`generate`, `prefix_quick`/`cmdline_quick` for `quick`) paired positionally with
  files in `ref_runs`/`ref_runs_quick`; diffs ignore volatile lines (timing, version banners) via
  `grep -v`.
- `ctest` in the CMake build itself defines no tests (`disorder`'s CMakeLists.txt has no
  `enable_testing()`/`add_test`) — the CI's `ctest` step is a no-op for this repo; the real check
  is the validation script.

When changing physics code (matrix elements, phase space, structure functions, DGLAP setup),
always run the validation script and treat any diff as something to explain, not silence.

## Architecture

Entry point `src/disorder.f90` drives two very different computational modes, selected by
`-p2b` / `inclusive` (default):

- **Inclusive mode** (default): cross sections are built directly from DIS structure functions
  (via Hoppet's `structure_functions` module) evaluated at fixed order, without generating
  explicit final-state radiation. Loops over PDF members here when `-pdfuncert` is requested.
- **P2B ("projection to Born") mode** (`-p2b`): adds real-radiation corrections through
  `DISENTFULL` (`src/libdisent.f`, ~3100 lines, adapted from Mike Seymour's `disent`/`dispatch`),
  which generates 2→3 and 2→4 real-emission kinematics and dipole-subtracted matrix elements on
  top of the Born. Currently NC-only, and capped below N3LO (see the guard clauses in
  `set_parameters`).

Both modes funnel differential contributions through `dsigma` (`src/mod_dsigma.f90`), the
function passed to the VEGAS integrator (`src/integration.f`). Per-call flow:

1. `gen_phsp_born` (`src/mod_phase_space.f90`) generates Born-level phase space (x, y, Q²) from
   VEGAS random numbers, with lab-frame ↔ Breit-frame boosts (`mlab2breit`/`mbreit2lab`).
2. `eval_matrix_element_new` (`src/mod_matrix_element.f90`) evaluates NC/CC structure functions
   (F1, F2, F3) at the requested perturbative order via Hoppet, applies electroweak couplings
   (`compute_sigmas_new`), and handles `xR`/`xF` scale variations (`muR_muF`).
2. Results are accumulated per scale-variation point (`sigma_all_scales`, indexed by the
   7-point `scales_mur`/`scales_muf` tables in `mod_parameters.f90`) and, in P2B mode, per
   real-emission multiplicity from disent.
3. If an analysis is linked in, `mod_analysis.f90`'s `analysis` routine hands Born/real/double-real
   momenta (lab and Breit frame) to the user-supplied `user_analysis`.

`mod_parameters.f90` is the central hub: all physical constants, EW parameters, and run-mode
flags live here as module-level state read once in `set_parameters` and used throughout
`src/*.f90` via `use mod_parameters`.

`src/toy_pdfs.f90` provides simple analytic starting PDFs used for validation/regression runs
independent of LHAPDF, evolved through Hoppet from a low scale — this is what the validation
suite's `-toyQ0` runs exercise.

### Analysis system

Analyses are separate Fortran files in `analysis/`, compiled in via the `-DANALYSIS=<file>`
cmake flag (default `simple_analysis.f`, a no-op skeleton). Each analysis implements
`user_analysis` and typically books/fills histograms via the POWHEG-derived booking library
(`analysis/pwhg_bookhist-multi.f`/`.h`, adapted from the POWHEG-BOX). FastJet-based analyses
(e.g. `exclusive_lab_frame_analysis.f`, `vbf.f`) need `-DNEEDS_FASTJET=ON`. Directories like
`panscales_sliceNLO_runs/`, `paper_runs/`, `vbf_analysis/`, and the various `test-*`/`test_*`
directories are run scripts and outputs for specific physics studies/papers, not part of the
core library — check their own scripts before assuming a shared convention.

## Third-party code boundaries

Do not casually "clean up" or restyle these — they are vendored/adapted, not house style:
- `src/libdisent.f`, `src/integration.f` — from `disent`/`dispatch` (Mike Seymour, GPLv3).
- `src/io_utils.f90`, `src/lcl_dec.f90` — command-line/IO utilities by Gavin Salam (GPLv3).
- `analysis/pwhg_bookhist-multi.*`, `aux/mergedata.f` — adapted from the POWHEG-BOX (GPLv2).
- Some `src/`/`analysis/` code is adapted from [proVBFH](https://github.com/fdreyer/proVBFH/).
