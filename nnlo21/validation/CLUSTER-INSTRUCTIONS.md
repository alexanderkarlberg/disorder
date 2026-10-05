# NNLO DIS 2+1 by τ₂ slicing: cluster production — instructions for Claude

**DRAFT (written 6 Oct 2026, night, on thA371a) — to be reviewed by AK before
use.** For the Claude session that AK starts on the cluster (Slurm, partition
`alma`, 24 h wall time, at most 10k running and 30k queued jobs). Read this
file, `CLAUDE.md`, `docs/nnlo21-plan.md` and the nnlo21 entries of
`docs/notebook.md` (3–6 Oct) first. The page with the current status:
https://claude.ai/artifact/ChKctwxDdEau7H9DGMxwyW

## 0. First: ask AK (one message, numbered)

1. Which disorder checkout: the branch `2026-10-nnlo21` was pushed on 5 Oct
   (afternoon); later commits are local on thA371a. Push them first (only with
   AK's OK).
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

### B. Our side in the same set-up (τ₂ slicing)

- LO: `sliced21 b0 … zeus` (fast).
- NLO coefficient: `sliced21 b1` + `nlo31 lo`, mode 2, τ_cut scan built in
  (10 columns, 2e-2 … 1e-5). On the MPP machines: converges to NNLOJET like
  √τ_cut (fiducial power corrections next to the cuts), all bins within 1.5σ
  at 1e-5.
- NNLO coefficient: `sliced21 b2` + `nlo31 vi, kp, r`, mode 2, psmc edge
  1e-12, technical cut 1e-9. **Open issue (6 Oct):** in this set-up the NNLO
  coefficient shows no plateau down to 1e-5 (total 11 → 86 ± 44 pb from 2e-3
  to 1e-5), concentrated in the bins next to the cuts. Either huge fiducial
  power corrections of the recoil-free projection or a problem in mode 2 at
  NNLO; NNLOJET's NNLO (A) decides. Do not over-invest in (B)-NNLO statistics
  before that is understood; a moderate set is enough to compare.

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

Job arrays of independent seeds; equal-weight seed averages (inverse-variance
weights of VEGAS errors are biased for R − D). Commands (MPP timings on a
thserv in brackets):
- `sliced21 b0 7500000 6 SEED zeus -` (~12 min)
- `sliced21 b1 3000000 6 SEED zeus <beamgrid prefix> <softtable>` (~25 min)
- `sliced21 b2 5000000 6 SEED zeus <beamgrid prefix> <softtable>` (~30 min)
- `nlo31 lo 20000000 6 SEED 1d-9 2 0 0 psmc 0 0 1d-12` (~20 min)
- `nlo31 vi 1000000 6 SEED 1d-9 2 0 0 psmc 0 0 1d-12` (~50 min)
- `nlo31 kp 2000000 6 SEED 1d-9 2 0 0 psmc 0 0 1d-12` (~5 min)
- `nlo31 r 1000000 6 SEED 1d-9 2 0 0 psmc 0 0 1d-12` (~40 min)
- D: replace `zeus` by `p2b` (sliced21) and the 6th argument `2` by `3`
  (nlo31); disorder: `disorder -p2b -nnlocoef -cutoff 1d-10 -pdf
  NNPDF30_nlo_as_0118 -Q2min 125 -Q2max 20000 -ymin 0.2 -ymax 0.6 -Ehad 920
  -ncall1 20000 -ncall2 40000000 -iseed SEED -prefix c2_` (~26 min), and the
  inclusive part without `-p2b` (seconds).
- Optional for D: correlated sampling `sliced21 c1 … p2b …` with
  `C1_NEMIT=-1` (stratified; up to ~3× more efficient at τ_cut = 1e-5).

## 4. Comparison

- Ours: `nnlo21/combine_zeus.py` (ZEUS), `slicing/p2b11_plot.py` (D, with
  the plot), `nnlo21/plot_tcut.py zeusnlo`, `nnlo21/plot_zeus_dist.py`.
- NNLOJET: average the per-seed histogram files with equal weights (or
  `nnlojet-run`'s combination), convert fb per unit to pb per bin.
- Report per bin: values, errors, pulls, χ², against τ_cut for ours.

## 5. Notes and push

A notes file in the results branch: what ran (commands, seeds, CPU), the
results, and any correction of earlier statements. Push only with AK's OK.
