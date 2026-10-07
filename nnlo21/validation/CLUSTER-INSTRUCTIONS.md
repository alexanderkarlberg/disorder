# NNLO DIS 2+1 by τ₂ slicing: cluster production — instructions for Claude

**Written 6 Oct 2026 on thA371a (draft of the night, updated in the evening
with projected slicing and the new technical cut).** For the Claude session that AK starts on the cluster (Slurm, partition
`alma`, 24 h wall time, at most 10k running and 30k queued jobs). Read this
file, `CLAUDE.md`, `docs/nnlo21-plan.md` and the nnlo21 entries of
`docs/notebook.md` (3–6 Oct) first. The page with the current status:
https://claude.ai/artifact/ChKctwxDdEau7H9DGMxwyW

## 0. First: ask AK (one message, numbered)

1. Which disorder checkout: the branch `2026-10-nnlo21` was pushed on 5 Oct
   (afternoon); everything below needs the commits of 6 Oct (local on
   thA371a, up to at least 347c248: projected slicing, the r fixes, the
   technical cut `TECHMIN`). They must be pushed first (only with AK's OK).
2. Work/scratch directory and quota; account, if any.
3. Modules: gfortran (≥ 9), LHAPDF 6 with the Fortran interface and the set
   `NNPDF30_nlo_as_0118`; hoppet (≥ 2.1.0, for disorder); MCFM-10.3 sources
   (the NNLO pieces are compiled from them, not stored in the repo); NNLOJET
   v1.0.2 (public, GPL; installed, or build it — needed for the reference runs
   and for `libnnlojet_core.so`, which `hard21` links). Does `nnlojet-run`
   support Slurm here, or write job arrays yourself?
4. Where to put the results: a new branch (e.g. `2026-10-nnlo21-cluster`) with
   the combined numbers, the scripts and a notes file; no per-job outputs.
5. The target precisions below, confirmed after the pilots, with the CPU
   estimates.
6. **Event shapes**: which ones (proposal in 1.C) and whether our side should
   have them before the production starts (they have to be added to
   `nlo31`/`sliced21` mode 2 first; that is done on thA371a, not on the
   cluster).

Rules: push only with AK's OK; report at the checkpoints (pilots and CPU
estimates before the full submissions; results before pushing); log the work
(notes file, section 5); when a result contradicts an earlier conclusion, say
so and correct the record. Nothing from the private dis-1jet code may be used
or committed. Kill only your own jobs, by job id.

## 1. Physics goals

All runs: photon exchange, α = 1/137, μ_R = μ_F = Q, NNPDF30_nlo_as_0118,
E_e = 27.5 GeV, E_p = 920 GeV.

### A. ZEUS-like dijets: high-quality NNLOJET references at LO, NLO and NNLO

Runcard `nnlo21/validation/nnlojet_epLJJ_zeus2j.run` (cuts of NNLOJET's ZEUS
dijet example: 125 < Q² < 20000 GeV², 0.2 < y < 0.6, Breit-frame k_t jets
R = 1 with E-scheme recombination, E_T = √(p_T² + m²) > 8 GeV, lab −1 < η <
2.5, ≥ 2 jets, m₁₂ > 20 GeV; histograms Q², p̄_T, m₁₂). Channels LO; R, V
(NLO coefficient); RR, RV, VV (NNLO coefficient).

Status on the MPP machines (5 Oct): LO 103.29 ± 0.05 pb (4 seeds), NLO
coefficient 10.49 ± 0.07 pb (R 20, V 10 seeds). Warmups: LO, V, R, VV done;
RV and RR slow (> 6 ms and > 3 ms per point on a thserv; the first 2M/4M-point
warmup iterations did not finish in 3 h).

Targets (per bin, to be confirmed with AK):
- LO: 1e-4 relative (fast).
- NLO coefficient: ≤ 0.3% of itself in the total, ≤ 1% per bin.
- NNLO coefficient: ≤ 2% of itself in the total, a few % per bin.

### B. Our side in the same set-up (τ₂ slicing, projected)

Terminology: **projected slicing** = P2B-improved τ₂ slicing of the 2+1
process itself (Campbell–Neumann–Vita 2408.05265, eq. 2.14): below τ_cut
the 3+1 events contribute O(event) − O(projected 2+1 Born) instead of
nothing (env `P2BSLICE=1` in `nlo31`, projection `project21`). It is not the
P2B step to N3LO 1+1 (that comes later, on top of this NNLO 2+1).

Mandatory settings for every `nlo31` run of parts lo, vi, kp, r in mode 2:
- `P2BSLICE=1`;
- psmc edge 1e-12 (12th argument), technical-cut argument `1d-9`;
- leave `P2BDROP` at its default (1).

**Technical cut: run production with the default cut and add the
correction as a separate difference run.** Do not run production with the
lowered cut `TECHMIN` directly. The default cut, s_min < 1e-9 W² (W² the
partonic invariant mass), biases τ_cut ≤ 3e-5 low. The correct cut,
min(1e-10 W², 1e-8 Q²), is validated. But a direct run with it is
contaminated by rare numerically garbage events: single seeds give r = 0 or
+2000 instead of −3260 (6 Oct, `runs/tcut11`). The difference runs below are
clean (mean ≈ median ≈ trimmed mean), and so is the default production. So:
- r production: default cut, no `TECHMIN`;
- r correction: `TECHDIFF=1 TECHMIN="1d-10 1d-8" VEGAS_EQUAL=2` (same
  arguments). This integrates only the events between the two cuts. Add its
  seed average to r, with errors in quadrature. About 60 seeds gave ±2.9 pb
  at 1e-5 (6 Oct).
- For every part, check that mean, median and trimmed mean of the seeds
  agree. If not, look for outliers before combining.
`sliced21` (b0, b1, b2) is unchanged.

Status on the MPP machines (6 Oct):
- **NLO:** projected slicing agrees with NNLOJET (10.49 ± 0.07 pb) in every
  bin from τ_cut ≤ 5e-4. Plain slicing drifts like √τ (fiducial power
  corrections) and needs about 3e-5.
- **NNLO:** projected slicing is flat from τ_cut = 2e-4 to 1e-5. Total
  40.1 ± 2.7, 41.7 ± 3.4, 42.4 ± 5.0, 43.2 ± 6.1 pb at 2e-4, 1e-4, 3e-5,
  1e-5. These used the default cut plus a measured `TECHDIFF` correction:
  the procedure to repeat. Plain slicing has errors 3–8× larger. There
  is no NNLOJET NNLO reference yet: that is the main goal of A.
- **Remaining checks for this production:**
  - the convergence of the technical cut at 1e-5: a residual up to
    ≈ 6 ± 4 pb is allowed;
  - the psmc edge, which shifts vi by about +2 ± 0.9 pb at small τ_cut.

Targets (to confirm with AK): NNLO coefficient to ≈ 1 pb in the total at
τ_cut = 1e-4 and 3e-5 (both on the plateau), a few % per bin.

### C. Event shapes (AK, 5 Oct)

NNLOJET's DIS process has `dis_thrust` (τ w.r.t. the boson axis; our τ_zQ),
`dis_thrust_c`, `dis_JB`, `dis_JB_c`, `dis_JM2`, `dis_C`, `y23`, `tau1`. At
NNLO in the 2+1 region (shape > a minimum) they come from the dijet process
(epLJJ) with the shape selected away from 0, as in NNLOJET's published NNLO DIS
event shapes. Proposal: τ_zQ (dis_thrust) in 0.05 … 0.5 as in our fixed-point
study, plus dis_C and y23, in the inclusive DIS cuts (125 < Q² < 20000,
0.2 < y < 0.6), each with its own minimum. Needs AK's choice and our mode-2
implementation first.

### D. One order lower: NNLO 1+1 by P2B on our τ₂-sliced 2+1

Lab-frame jets (`analysis/lab11_analysis.f`: anti-k_t R = 1, p_T > 5 GeV, −1 <
y < 2.5; leading-jet p_T and y). Ours: disorder inclusive (`-nnlocoef`) +
`sliced21 b1` + `nlo31 lo` in mode 3 (P2B). Reference: `disorder -p2b
-nnlocoef -cutoff 1d-10` with the same analysis. On the MPP machines: χ² ≈ 15
bins for τ_cut ≤ 2e-4, errors ≤ 0.2% of the jet rate (ours) and ≤ 0.13%
(disorder). Target: per mille of the jet rate per bin on both sides (AK);
the NLO step already agrees to < 0.01%.

## 2. Building

1. disorder (with the lab-frame analysis for D):
   `cmake -S . -B build-lab11 -DANALYSIS=lab11_analysis.f -DDISORDER_TESTS=OFF && cmake --build build-lab11 -j`.
2. `sliced21` and `nlo31`: `nnlo21/build_sliced21.sh <MCFM-10.3> <dir of libnnlojet_core.so> <builddir>`.
3. Tables (once, shared by all jobs):
   - beam grid: `sliced21 mktab i i <prefix>` for i = 0 … 28 (29 jobs of
     ~1.5 min), Q from √125 to √20000 GeV;
   - soft table: built on first use (`<softtable>` argument; 80 s, 19 MB).
4. Check, before any production: modes 0–2 of the new binaries reproduce the
   results in `docs/notebook.md` (e.g. ZEUS LO 103.29 pb with `b0`).

## 3. Runs

### NNLOJET (A)

Copy the runcard per channel (`CHANNEL`, `SEED`, `WARMUP`, `PRODUCTION`).
Pilot first: warmups LO 1M[5], R and V 2M[5], VV 1M[5], RV and RR as large as
fits in 24 h on one core (NNLOJET can share a warmup grid between productions:
check `nnlojet-run`); then 20 production jobs per channel to measure time and
error per point; report and size the production with AK.

### Ours (B, D)

Job arrays of independent seeds. Combine seeds with equal weights;
inverse-variance weights of VEGAS errors are biased for R − D. Every `nlo31`
job of B needs the environment

    export P2BSLICE=1

and the r correction jobs in addition `TECHDIFF=1 TECHMIN="1d-10 1d-8"
VEGAS_EQUAL=2` (section 1.B).

Commands, with MPP timings on a thserv in brackets:
- `sliced21 b0 7500000 6 SEED zeus -` (~12 min)
- `sliced21 b1 3000000 6 SEED zeus <beamgrid prefix> <softtable>` (~25 min)
- `sliced21 b2 5000000 6 SEED zeus <beamgrid prefix> <softtable>` (~30 min)
- `nlo31 lo 3000000 6 SEED 1d-9 2 0 0 psmc 0 0 1d-12` (~25 min)
- `nlo31 vi 1000000 6 SEED 1d-9 2 0 0 psmc 0 0 1d-12` (~50 min)
- `nlo31 kp 2000000 6 SEED 1d-9 2 0 0 psmc 0 0 1d-12` (~5 min)
- `nlo31 r 1000000 6 SEED 1d-9 2 0 0 psmc 0 0 1d-12` (~45 min). This part
  dominates the error: about 300 seeds give ±3.4 pb at τ_cut 1e-4.
- r correction (`TECHDIFF`, see above): about 100 seeds.
- Cross-checks, smaller sets (about 50 r seeds each):
  - cut convergence: `TECHDIFF=1 TECHMIN="1d-11 1d-9" TECHREF="1d-10
    1d-8" VEGAS_EQUAL=2`, must be ≈ 0 at 3e-5 and 1e-4 (6 Oct: 1.0 ± 3.8
    and −0.1 ± 0.3; 6.1 ± 3.6 at 1e-5);
  - vi and r with psmc edge 1e-13 (edge);
  - plain slicing, i.e. the same commands without `P2BSLICE`, with its own
    `TECHDIFF` correction (comparison).
- The fixed-point τ_cut test (mode 1, x = 0.01, Q² = 400, `nlo31 … 1d-9 1
  0.01 400 psmc 0 0 1d-12`, `sliced21 … 0.01 400`) with high statistics,
  plus its `TECHDIFF` correction. MPP status (6 Oct): plateau 8.5 ± 1.2
  pb/GeV² from 5e-4 to 1e-4; corrected 13 ± 3 at 3e-5 and 15.5 ± 5.5 at
  1e-5, i.e. 1.3–1.5σ above the plateau. Resolving this is part of the job.
- D: replace `zeus` by `p2b` (sliced21) and the 6th argument `2` by `3`
  (nlo31); there `P2BSLICE` does not apply (mode 3 is already the P2B to
  1+1). Add a `TECHDIFF` correction run as well. disorder: `disorder -p2b -nnlocoef -cutoff 1d-10 -pdf
  NNPDF30_nlo_as_0118 -Q2min 125 -Q2max 20000 -ymin 0.2 -ymax 0.6 -Ehad 920
  -ncall1 20000 -ncall2 40000000 -iseed SEED -prefix c2_` (~26 min), and the
  inclusive part without `-p2b` (seconds).
- Diagnostics, if something looks off. Validate each on a known answer first
  (see notebook 6 Oct):
  - `TECHDIFF=1 TECHREF=…`: the difference between two technical cuts,
    same events;
  - `P2BEXTRA=1` with `VTARGET=1 VEGAS_EQUAL=2`: the projected-slicing term
    alone;
  - `P2BDEBUG=1`.

Before production: check that each new binary reproduces the MPP numbers,
e.g. ZEUS LO 103.29 pb (b0), the NLO with projected slicing, and a short r
run against the notebook.

## 4. Comparison

- Ours: `nnlo21/combine_zeus.py` (ZEUS), `slicing/p2b11_plot.py` (D, with
  the plot), `nnlo21/plot_tcut.py zeusnlo`, `nnlo21/plot_zeus_dist.py`.
- NNLOJET: average the per-seed histogram files with equal weights (or
  `nnlojet-run`'s combination), convert fb per unit to pb per bin.
- Report per bin: values, errors, pulls, χ², against τ_cut for ours.

## 5. Notes and push

A notes file in the results branch: what ran (commands, seeds, CPU), the
results, and any correction of earlier statements. Push only with AK's OK.
