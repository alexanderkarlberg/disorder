# disorder news

Release notes for each version, newest first. The detailed list of changes since 1.0.0, with pointers into the code, is in [docs/release-notes.md](docs/release-notes.md).

## v2.2.0 (2026-09-30)

This release extends the projection-to-Born (P2B) mode, `-p2b`, from photon exchange to all neutral- and charged-current processes. DISENT's matrix elements now include Z and W exchange up to O(αs²). Fully differential predictions (jets, event shapes) at NLO and NNLO are therefore available for NC γ/Z exchange (also Z only or interference only) and CC, for charged-lepton and neutrino beams of either charge. For photon exchange the results are bitwise identical to v2.1.1.

The release also adds unit tests, runs the validation matrix through `ctest` with tolerance-based comparisons, and extends CI accordingly.

### Features

* **P2B for all NC and CC processes, up to NNLO.** `-p2b` can now be combined with `-includeZ` (and `-Zonly`, `-intonly`), `-CC`, `-neutrino` and `-positron`. The electroweak couplings are the same as in the inclusive (structure-function) mode.
  * O(αs): the three-parton tree-level matrix element with γ/Z/W exchange.
  * O(αs²): the four-parton tree-level matrix elements, including the identical-quark and W-attachment interferences; the spin-correlated dipoles; and the finite one-loop part. The one-loop parity-violating part uses the Bern–Dixon–Kosower amplitudes, ported from MCFM.
  * Validated against NNLOJET 1.0.2 and POWHEG-BOX-RES at O(αs) for every process these codes support, against NNLOJET at O(αs²) for NC γ, NC Z and CC e±, and point by point against an independent implementation of the matrix elements for all processes and beams.
  * As in the inclusive coefficient functions, contributions proportional to the sum of the quark axial couplings are neglected. They vanish for complete generations, so with 5 flavours only the b quark contributes. The T-odd parts of the one-loop amplitudes are dropped; they integrate to zero for observables that are symmetric under reflection.
  * `-p2b` is still not available at N3LO, with a variable flavour number, or with `-pdfuncert`.
  * When running `-p2b` at NNLO with Z or W exchange, please also cite Z. Bern, L.J. Dixon, D.A. Kosower, Nucl. Phys. B 513 (1998) 3 [hep-ph/9708239] (see the citation policy in the README).
* **New analyses.**
  * `analysis/cmp_nnlojet_powheg.f`: the analysis used for the comparisons with NNLOJET and POWHEG-BOX-RES. It needs no FastJet. It fills Q², x, y, the leading lab-frame jet and the Breit-frame event shapes τ_zE, B_zE and ρ_E, also with a current-hemisphere energy cut E_cur > Q/10 that makes τ_zE and B_zE infrared safe at O(αs²).
  * `analysis/inclusive_panscales_nlo_analysis.f`: used for checks in the DIS PanScales matching paper.
* **Tests and validation.**
  * `ctest` runs unit tests (about a second, `ctest -L unit`) of:
    * the phase space and frame transformations;
    * the matrix element against quark-parton-model expressions;
    * the interface to DISENT;
    * DISENT's three-parton matrix element for every process;
    * the consistency of the four-parton matrix element with its dipole subtraction terms in all single-unresolved limits for every process;
    * the refusal of unsupported option combinations.
  * Each entry of the validation matrix (`validation/configurations.txt`) is a `ctest` test. It is compared with the stored references to a relative tolerance (default 1e-5), so it no longer fails on a different compiler or platform. A NaN or infinity in any output is a failure.
  * The matrix now covers P2B at NLO and NNLO for NC γ/Z, Z only, interference only and CC, for e±, ν and ν̄ beams, and NC+CC combined, with and without scale variations.
  * `validation/validate_or_generate.sh validate|quick|generate` uses `ctest` and no longer needs GNU parallel.
* **CI** runs the unit tests and the full validation matrix, and keeps its scheduled runs from being disabled after 60 days without repository activity.

### Bug fixes

* HOPPET's grid now also covers x up to 0.1 when a very large `-xmin` is requested.
* `-scale-choice 4`: the printed header and the `-help` text gave the central scale as Q(1−x)/x. The scale used is, as intended, μ² = Q²(1−x)/x. Results are unchanged.
* The validation scripts have been fixed, and `quick_validation.sh` is now the `quick` mode of `validate_or_generate.sh`.

### Other changes

* The reference runs have been updated to HOPPET v2.3.0.
* The DISENT callbacks and the PDF/structure-function set-up moved to their own modules (`src/mod_disent_interface.f90`, `src/mod_pdf_setup.f90`). The physics code is built as a CMake object library shared by `disorder` and the tests. CMake 3.12 or newer is required.
* The README now describes the two modes, the tests, and the third-party code (including the MCFM port).

**Full Changelog**: https://github.com/alexanderkarlberg/disorder/compare/v2.1.1...v2.2.0

## v2.1.1 (2026-01-22)

This is a minor release which mainly rectifies that the v2.1.0 was given v2.0.1 as a tag. Now the tag and version are consistent again as v2.1.1. Additionally reference runs ahve been updated to use hoppet v2.2.0.

### Features
The code now compiles with Library Time Optimisation (LTO). This seems to speed up certain runs by ~10%.

### Bugs fixes
There are no bug fixes.

**Full Changelog**: https://github.com/alexanderkarlberg/disorder/compare/v2.0.1...v2.1.1

## v2.1.0 (2025-10-13, tagged v2.0.1)

This release contains a few modifications to align with v2 of HOPPET. In particular disorder now uses the streamlined hoppetSetCoupling to set up the coupling. It also uses the faster interpolation orders in hoppet by default. They can be changed on the command line.

**Full Changelog**: https://github.com/alexanderkarlberg/disorder/compare/v2.0.0...v2.0.1

## v2.0.0 (2025-08-12)

This update breaks backwards compatibility because it now requires hoppet v2.0.0 or newer, and will not compile if that is not present. Part of that change now means that there is support for using the N3LO evolution in hoppet, but otherwise this is to be considered a small update.

**Full Changelog**: https://github.com/alexanderkarlberg/disorder/compare/v1.1.0...v2.0.0

## v1.1.0 (2024-09-17)

This is mainly a bug release but with a few new features. Given the accumulation of various minor new features it seems timely to increase the minor version number!

### Bugs
* toyPDFs did not work properly with p2b and NO scale variations
* mergedata would cut long command lines upon combination

### Features
* New analysis (caesar.f) that computes the broadening and thrust, relevant for extracting logs for these event shapes. The event shapes are taken from EvtShpLib in CAESAR.
* Updated mergedata to the version that can accumulate histograms (in addition to many other features)
* Put in protection against NaN in the analysis. NaN can happen when running with a very small cutoff and large values of npow. The program prints the offending event whenever this happens.
* Updated help and README.md

**Full Changelog**: https://github.com/alexanderkarlberg/disorder/compare/v1.0.2...v1.1.0

## v1.0.2 (2024-08-20)

This is mainly a bug release with a few minor new features:

### Bugs (all fixed)
* alphaEM was not set correctly in disent when NOT running with -scaleuncert (reported by Silvia FR)
* -p2b together with  -order-min/-order-max did not have the expected behaviour (reported by Silvia FR)
* If the user specified both an order (e.g. -nnlo) and -order-min/-order-max the program would ignore the latter (reported by Melissa v B)
* The program did not run/compile if the user had installed hoppet or lhapdf outside of their path (reported by Melissa v B)
* analysis.f had a bug (reproted by Melissa v B)

### Features
* Added an nf=5 toy PDF initial condition and support on the command line for picking it (-pdf toyHERALHC || toyNF5) and setting the alphas value at Q0 (-toy-alphas-Q0).
* Quark masses can now be set on the command line when using a toy PDF (-mc, -mb, -mt)
* Added support on the command line for flags that specify a specific perturbative coefficient (-nlocoef, -nnlocoef, -n3locoef)

**Full Changelog**: https://github.com/alexanderkarlberg/disorder/compare/v1.0.1...v1.0.2

## v1.0.1 (2024-07-30)

This is a bug fix release. 

1. The Variable Flavour Number Scheme was not implemented correctly, which has now been fixed. At the same time the quark masses are now read directly from the PDF (if using a toy PDF the internal HOPPET values are used by default).

2. The header is now printed with more digits.

3. Validation runs updated.

**Full Changelog**: https://github.com/alexanderkarlberg/disorder/compare/v1.0.0...v1.0.1

## v1.0.0 (2024-06-17)

This first stable release is published with the acceptance of the paper in SciPost. No changes from the final beta release.

## v1.0.0-beta.4 (2024-06-05)

This release corresponds to the resubmission to SciPost after the first editor decision.  

The main changes are

*) The documentation has been updated to reflect more clearly the limitations of p2b.
*) The code support neutrino beams in the inclusive mode (with the flag -neutrino or -neutrino -positron for anti-neutrino).

**Full Changelog**: https://github.com/alexanderkarlberg/disorder/compare/v1.0.0-beta.3...v1.0.0-beta.4

## v1.0.0-beta.3 (2024-04-12)

This is a bug fix beta release. 

### Bug fixes
* The alpha_em set on the command line was not passed to disent at all, and therefore in p2b mode results were not correct if the commandline had a different value from 1/137.
* fastjet was acutally not included correctly in CMakeLists.txt if the user tried to compile having installed fastjet from a non-standard path.

### Updates
* help output has been updated and made more useful
* gev2pb has been updated to most up2date value from PDG
* command-line now supports Q2min and Q2max arguments in addition to Qmin and Qmax.
* disorder now uses hoppet αS running everywhere and not a mix of lhapdf and hoppet. This has an impact if one is using a PDF with a different order than what one is calling disorder with as in that case the hoppet and LHAPDF running would differ.

**Full Changelog**: https://github.com/alexanderkarlberg/disorder/compare/v1.0.0-beta2...v1.0.0-beta.3

## v1.0.0-beta.2 (2024-02-01, tagged v1.0.0-beta2)

Added correct arXiv identifier and fixed a typo in the header.

## v1.0.0-beta.1 (2024-01-30, tagged v1.0.0-beta1)

This is the release created at the time of submission to arXiv.
