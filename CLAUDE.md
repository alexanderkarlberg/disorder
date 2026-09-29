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

## Tests

Everything runs through `ctest` (`enable_testing()` in the top-level CMakeLists.txt,
`-DDISORDER_TESTS=OFF` to skip). There are two kinds of test:

- **Unit tests** (`tests/`, label `unit`, ~1 s in total): small Fortran programs linked against
  the physics code (`disorder_core` object library) plus `analysis/simple_analysis.f`, using the
  asserts in `tests/test_utils.f90` (a failed check exits non-zero):
  - `test_boosts`, `test_phase_space`: `mlab2breit`/`mbreit2lab`, `gen_phsp_born`, `p2bmomenta`
    (kinematics, frame properties, Jacobian integrated against the phase-space volume).
  - `test_matrix_element`: run with disorder's own flags (one ctest entry per process: γ/Z/γZ,
    e±, ν/ν̄, NC/CC). At LO it compares `eval_matrix_element_new` with parton-model formulas
    written out independently in the test, and also the per-parton couplings `parton_couplings`
    used by DISENT; beyond LO (photon only) it checks the F2/FL combination. Also checks
    `muR_muF` and the scale-variation labels.
  - `test_matthr`: DISENT's three-parton matrix element for each process (same flag sets):
    bitwise identical to the original DISENT expression for photon exchange, decomposition into
    parity-conserving/violating parts, collinear limits.
  - `test_subtraction`: DISENT's four-parton matrix element (MATFOR) against the sum of its
    dipoles (SUBFOR) in collinear, soft and initial-state limits, for each process (symmetrised
    over the labels of the outgoing partons, dipoles with an unresolved Born left out).
  - `test_disent_interface`: drives the real `DISENTFULL` with a checking callback and verifies
    the momentum mapping/boosts used by `user` (`src/mod_disent_interface.f90`).
  - `guard_*`: invalid flag combinations must be refused. The guards use `stop`, which exits
    with status 0, so these tests match the message.
- **Validation matrix** (`validation/`, labels `validation` + `fast`/`slow`): one test per line of
  `validation/configurations.txt` (`labels | prefix | args`). `validation/run_validation.py`
  runs disorder and compares every output file with `validation/ref_runs/<prefix>*` with a
  relative tolerance (default 1e-5), after dropping volatile lines (timings, banners,
  citations). Only registered when configured with
  `-DNEEDS_FASTJET=ON -DANALYSIS=exclusive_lab_frame_analysis.f` (the reference histograms
  come from that analysis); runs whose LHAPDF set is missing are disabled.

```
cmake -S . -B build -DNEEDS_FASTJET=ON -DANALYSIS=exclusive_lab_frame_analysis.f
cmake --build build -j
ctest --test-dir build -j 8 --output-on-failure        # everything (~5 min, dominated by one run)
ctest --test-dir build -L unit                         # unit tests only
ctest --test-dir build -LE slow -j 8                   # skip the slow MSHT20 PDF-uncertainty run
ctest --test-dir build -R validation_inclusive_cc_Q_10 --output-on-failure   # a single test
```

`validation/validate_or_generate.sh validate|quick|generate` wraps the same thing in a fresh
`validation/build` (extra cmake flags via `EXTRA_CMAKEFLAGS`). `generate` overwrites
`ref_runs/`: only do this deliberately, when a physics/numerical change is intended and
reviewed, since it redefines "correct". CI (`.github/workflows/cmake-single-platform.yml`,
dependencies from `.github_CI.sh`) configures with FastJet and the exclusive analysis and runs
the full `ctest`. A scheduled `keepalive` job re-enables the workflow through the GitHub API so
that the 60-day inactivity rule doesn't switch off the scheduled runs.

The log files in `ref_runs/` hold a run's stdout followed by its stderr (as GNU parallel,
used for the original references, wrote them), and `run_validation.py` writes them the same way.
The FastJet and LHAPDF banners come from C++ streams, so their position is not stable and they
are ignored.

When changing physics code (matrix elements, phase space, structure functions, DGLAP setup),
always run the full `ctest` and treat any failure as something to explain, not silence (e.g. by
loosening a tolerance).

## Notebook and release notes

- `docs/notebook.md`: dated log of development work (what was checked, how, findings, and
  explicit corrections of earlier conclusions). Append an entry for each piece of work.
- `NEWS.md`: the release notes of every version (the GitHub release texts), newest
  first; add the entry for a new version there when releasing.
- `docs/release-notes.md`: user-facing changes since 1.0.0. The 1.0.0 manual
  (`docs/disorder-1.0.0-manual.tex`) matches the published paper and is not edited; document
  changes in the release notes instead.

## Architecture

Entry point `src/disorder.f90` drives two very different computational modes, selected by
`-p2b` / `inclusive` (default):

- **Inclusive mode** (default): cross sections are built directly from DIS structure functions
  (via Hoppet's `structure_functions` module) evaluated at fixed order, without generating
  explicit final-state radiation. Loops over PDF members here when `-pdfuncert` is requested.
- **P2B ("projection to Born") mode** (`-p2b`): adds real-radiation corrections through
  `DISENTFULL` (`src/libdisent.f`, ~3100 lines, adapted from Mike Seymour's `disent`/`dispatch`),
  which generates 2→3 and 2→4 real-emission kinematics and dipole-subtracted matrix elements on
  top of the Born. All NC/CC processes are supported up to NNLO (N3LO is not available; see
  the guard clauses in `set_parameters`). The electroweak couplings of DISENT's matrix
  elements come from `parton_couplings(_split)` (`src/mod_matrix_element.f90`, mirroring
  `eval_matrix_element_new`) via `disent_couplings(4)` (`src/mod_disent_interface.f90`);
  photon exchange is bit-identical to the original DISENT. The parity-violating and CC pieces
  are in `src/disent_o2_trees.f` (generated by `derivations/o2/make_trees.py` from FORM; do not
  edit by hand) and `src/disent_virt3.f` (one loop, BDK amplitudes ported from MCFM).
  Derivations and validation: `derivations/matthr/`, `derivations/o2/`.

The disent callbacks (`user`, `dis_cuts`, `disent_muf`) live in `src/mod_disent_interface.f90`,
and HOPPET/PDF/structure-function set-up (`read_PDF`, `setup_structure_functions`) in
`src/mod_pdf_setup.f90`, so that the tests can use them.

Both modes funnel differential contributions through `dsigma` (`src/mod_dsigma.f90`), the
function passed to the VEGAS integrator (`src/integration.f`). Per-call flow:

1. `gen_phsp_born` (`src/mod_phase_space.f90`) generates Born-level phase space (x, y, Q²) from
   VEGAS random numbers, with lab-frame ↔ Breit-frame boosts (`mlab2breit`/`mbreit2lab`).
2. `eval_matrix_element_new` (`src/mod_matrix_element.f90`) evaluates NC/CC structure functions
   (F1, F2, F3) at the requested perturbative order via Hoppet, applies electroweak couplings
   (`compute_sigmas_new`), and handles `xR`/`xF` scale variations (`muR_muF`).
3. Results are accumulated per scale-variation point (`sigma_all_scales`, indexed by the
   7-point `scales_mur`/`scales_muf` tables in `mod_parameters.f90`) and, in P2B mode, per
   real-emission multiplicity from disent.
4. If an analysis is linked in, `mod_analysis.f90`'s `analysis` routine hands Born/real/double-real
   momenta (lab and Breit frame) to the user-supplied `user_analysis`.

`mod_parameters.f90` is the central hub: all physical constants, EW parameters, and run-mode
flags live here as module-level state read once in `set_parameters` and used throughout
`src/*.f90` via `use mod_parameters`.

The user manual source is `docs/disorder-1.0.0-manual.tex`, and `docs/disent_manual.pdf` is the manual for the original disent.

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
- `src/disent_virt3.f` — ported from MCFM 10.3 (BDK one-loop amplitudes; GPLv3 or later).
- `src/disent_o2_trees.f` — generated (`derivations/o2/make_trees.py`); regenerate rather than edit.
