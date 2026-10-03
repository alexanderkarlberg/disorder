# NLO DIS 3+1 validation on the cluster: instructions for Claude

**Done without the cluster (3 Oct 2026, MPP machines):** NNLOJET RV + RR =
29.99 ± 0.68 pb against nlo31's 30.34 ± 0.29 pb (0.5σ), six Q² bins
χ² = 9.2/6; see `docs/notebook.md` (3 Oct, "NLO 3+1 validated"). The
instructions below are kept for a more precise or extended run (Z/W,
distributions).

Written 3 Oct 2026 (night, on thA371a) for the Claude session that AK starts
on the cluster (Slurm, partition `alma`, 24 h wall time, at most 10k running
and 30k queued jobs). Read this file, `CLAUDE.md`, `docs/dis31-plan.md` and
the dis31 entries of `docs/notebook.md` (2–3 Oct) first.

## 0. First: ask AK (one message, numbered)

1. Which disorder checkout to use. The work is on the branch
   `2026-10-dis31`; on 3 Oct it was local only. Has it been pushed, or
   should you push it (only with AK's OK)?
2. The work/scratch directory and its quota; account, if any.
3. Compilers/modules: gfortran (>= 9), LHAPDF 6 (with its Fortran
   interface), the set `NNPDF30_nlo_as_0118`. NNLOJET v1.0.2: is it
   installed, or should you build it from the public tarball
   (nnlojet.hepforge.org; built on thA371a with gfortran, see
   `~/cernbox/disorder-comparisons/NOTES.md` there)? Its `nnlojet-run`
   workflow tool: does it support Slurm on this cluster, or should you
   write job arrays yourself?
4. Where to put results: a new branch, e.g. `2026-10-dis31-validation`,
   with the combined numbers, the scripts and a notes file; no per-job
   outputs.
5. The target precision below, confirmed after the pilot, with its CPU
   estimate.

Rules: push only with AK's OK; report at the checkpoints (pilot and CPU
estimate before the full submission; results before pushing); log the work
(notes file, section 5); when a result contradicts an earlier conclusion,
say so and correct the record. Nothing from the private dis-1jet code may be
used or committed. Kill only your own jobs, by job id.

## 1. Physics goal

Validate the integrated NLO 3+1 of disorder's `dis31` (DIS 3+1 at NLO with
Catani–Seymour subtraction; photon exchange) against NNLOJET's DIS dijet
NNLO, where the ≥ 3-jet cross section at O(α_s³) is exactly the NLO 3-jet
correction (RV + RR; VV and the double-unresolved terms vanish with
njets ≥ 3).

Set-up (identical in `dis31/nlo31.f90` and
`dis31/validation/nnlojet_epLJJ_3j.run`): e⁻ p, 27.5 × 920 GeV; photon only;
α = 1/137; 150 < Q² < 15000 GeV², 0.1 < y < 0.9; inclusive kt jets, R = 1,
E-scheme, in the Breit frame, p_T > 5 GeV, ≥ 3 jets; μ_R = μ_F = Q;
NNPDF30_nlo_as_0118 (α_s from the set); n_f = 5. Observables: σ(≥ 3 jets)
and six Q² bins (edges 150, 200, 300, 500, 1000, 3000, 15000).

Status on thA371a (3 Oct):

| | nlo31 [pb] | NNLOJET [pb] |
|---|---|---|
| LO (3+1 tree) | 84.852 ± 0.092 | 84.739 ± 0.099 (channel R) |
| V + I | 28.881 ± 0.072 | |
| K + P | 31.092 ± 0.035 | |
| R − dipoles | −30.025 ± 0.455 | |
| NLO correction | 29.95 ± 0.46 | to do: RV + RR |

LO agrees in the total (0.8σ) and in all six Q² bins (≤ 0.5%). NNLOJET's
RV alone is ≈ 56 pb (its antenna terms distribute differently; only RV +
RR as a whole compares with ours). NNLOJET's RR is slow (about 0.1 s per
point on a loaded thA371a core, several sub-channels RRa, ...), hence the
cluster.

Target: the NLO correction to ±0.3 pb or better (1%) from each code, the
Q² bins to about 2%.

## 2. Building

- nlo31: `dis31/build_nlo31.sh <dir> [path to lhapdf-config]` (dis31
  sources and LHAPDF only). Usage: `nlo31 part ncall itmx seed [techcut]`,
  part = lo, vi, kp, r. Each job does its own VEGAS adaptation; the first
  iteration is not used in the result.
- The unit tests of dis31 need a disorder build (`dis31/tests/build.sh`);
  optional on the cluster, but run `harness_dip41` and `harness_virt31` once
  if you build disorder there.
- NNLOJET: runcard template `dis31/validation/nnlojet_epLJJ_3j.run`
  (replace SEED, WARMUP, PRODUCTION, CHANNEL). Pitfalls (see the runcard
  header): no `dis_frame = BREIT` (zero or NaN in this build; Breit is the
  default for DIS), and `collider = ep` instead of `beam1 = ...`.

## 3. Runs

- nlo31 cost (thA371a core, per VEGAS point): vi about 2 ms (finite part
  only; 2.8 times faster than with the poles, same result), r about
  0.6–1 ms, lo and kp negligible. Measure in the pilot. On thA371a 4 × 6 M points of r gave
  ±0.45 pb (seed scatter consistent with the VEGAS errors). For ±0.15 pb on
  r: about 40 jobs of 6 M points (`nlo31 r 1000000 6 <seed>`), about
  2 CPU-h each. vi: 8 jobs of `nlo31 vi 200000 6 <seed>`. kp: 2 jobs of
  `nlo31 kp 1000000 8 <seed>`. Use distinct seeds; combine with
  `dis31/combine_nlo31.py` (inverse-variance per part, then vi + kp + r).
- Technical-cut check: r with `techcut` 1e-7 and 1e-11 (default 1e-9)
  must agree (thA371a, 3 Oct: −29.35 ± 0.50, −30.02 ± 0.46, −29.06 ± 0.48 pb
  for 1e-7, 1e-9, 1e-11; consistent).
- NNLOJET: a pilot of RV and RR (e.g. 20 jobs each, warmup 100000[4],
  production 1000000[1]) to measure the time per point and the error per
  point; then size the production for ±0.3 pb on RV + RR and report to AK.
  NNLOJET reports fb (and dσ/dQ² in fb/GeV² for the histogram).
- Optional, cheap: NNLOJET channel R at high statistics for a sharper LO
  check.

### MPP machines (instead of, or in addition to, the cluster)

AK (3 Oct): the institute machines may be used, and should be scanned
first:
- thservs (thserv05–24, 64 logical cores each): nice 10, at most 30 of
  your jobs per machine, and keep the machine's total load at or below ~32
  (beyond that hyperthreading slows every job);
- other th machines (desktops; fewer cores but much faster): nice 19,
  never more than half of the logical cores, and not more than the cores
  that are actually free (several run other people's jobs);
- AK's own thA371a: never niced.
A dispatcher that respects these limits and starts queued jobs every 10
minutes: `dis31/validation/dispatch.sh <queue dir> <command>` (log next to it) on
thA371a (start it from an interactive session so that ssh has Kerberos
credentials; it skips hosts without the NFS mount of /home/thA371a, where
`ssh -f` would otherwise report success for a job that never started).
After a launch, check that the jobs really run (state R, not T: the
desktops' overheatd stops hot jobs) and audit running processes against the
queue.

## 4. Comparison

- σ(≥ 3 jets): nlo31 vi + kp + r against NNLOJET RV + RR; also LO + NLO.
- Q² bins: nlo31 prints σ per bin in pb; NNLOJET dσ/dQ² × bin width/1000.
- Agreement within the combined errors is the success criterion. A
  discrepancy: check first the set-up (runcard pitfalls), then the
  technical cut, then the per-part pieces we can compare (LO). Report it to
  AK before chasing it further.

## 5. Notes and push

Notes file `dis31/validation/NOTES-cluster.md` (what was run, CPU used,
results, comparison). Push (after AK agrees) the notes, the combined
numbers, the job scripts; not the per-job outputs (give AK their location).
