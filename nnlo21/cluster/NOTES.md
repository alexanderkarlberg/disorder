# nnlo21 cluster production (MPCDF t2, 6–7 Oct 2026)

Claude session on mppui1 for AK, following
`nnlo21/validation/CLUSTER-INSTRUCTIONS.md`. Branch `2026-10-nnlo21-cluster`
(from `2026-10-nnlo21` at d4b12b6), committed locally, NOT pushed.

Work directory (raw outputs, not committed): `/ptmp/mpp/akarlber/disorder-nnlo21/`
(`P` below).

## AK's answers to section 0 (6 Oct, ~22:30)
1. Checkout `2026-10-nnlo21` at d4b12b6 (contains 347c248).
2. Work/scratch: `$P`; default Slurm account.
3. Install what is needed (NNPDF30_nlo_as_0118, MCFM-10.3, NNLOJET v1.0.2).
4. Results in this branch; no per-job outputs.
5. Targets as in the instructions; size and submit after the pilots; reduce and
   flag any item that would need > 300k core-hours.
6. Event shapes: τ_zQ (dis_thrust) 0.05–0.5, dis_C, y23 in the inclusive DIS
   cuts. Our mode-2 implementation is done on thA371a; here only A, B, D, and
   the NNLOJET event-shape runcards/pilot if cheap.
7. Shared cluster: never touch the proVBFH jobs (`hxswg136-*`, `p1506-*`,
   `vbfnlo-*`, `h3jpw-*`). All jobs here are named `dis-*`, submitted with
   `--nice=100`.

## Software (what was built and where)
- `NNPDF30_nlo_as_0118` from lhapdfsets.web.cern.ch into `~/.local/share/LHAPDF`.
- MCFM-10.3 public tarball (mcfm.fnal.gov) unpacked in `$P/src/MCFM-10.3`.
- NNLOJET v1.0.2 (https://nnlojet.hepforge.org/nnlojet-v1.0.2.tar.gz):
  - `$P/src/nnlojet-v1.0.2/install` (serial build; its `lib64/libnnlojet_core.so`
    is the one `sliced21` links);
  - `$P/src/omp/nnlojet-v1.0.2/install` (same source, `-DOPENMP=ON`), linked as
    `$P/nnlojet`; used for all NNLOJET runs (warmups multithreaded, production
    with one thread; with one thread it is bit-identical to the serial build,
    checked on a small LO production).
  - Built on a compute node (16 cores, 6.5 min).
  - LHAPDF 6.5.6's `lhapdf-config` advertises `--incdir` but only accepts
    `--includedir`; NNLOJET's and disorder's CMake call `--incdir`. Wrapper
    `$P/src/lhapdf-wrap/lhapdf-config` maps it (used for both builds).
  - NNLOJET needs `ulimit -s unlimited` and `OMP_STACKSIZE=1024M` (the OpenMP
    build segfaults otherwise, even with one thread), and the sweep sizes as
    plain integers (`1M[5]` gives "parse_run_io: unknown option").
- disorder with the lab-frame analysis: `$P/builds/build-lab11`
  (`cmake -DANALYSIS=lab11_analysis.f -DDISORDER_TESTS=OFF`, LHAPDF via the
  wrapper).
- `sliced21`, `nlo31`: `nnlo21/build_sliced21.sh $P/src/MCFM-10.3
  $P/src/nnlojet-v1.0.2/install/lib64 $P/builds/sliced21`.
- Frozen copies used by all jobs: `$P/bin/{sliced21,nlo31,disorder}` (md5 in
  `$P/bin/MD5SUMS`).
- Tables: beam grid `$P/tables/beamgrid/zeus_{0..28}.tab` (29 jobs `sliced21
  mktab i i`, < 2 min each), soft table `$P/tables/soft.tab` (19 MB, built by
  a tiny b1 run on the login node, 4 min).

## Job machinery (`nnlo21/cluster/`)
- `mkjobs.py LIST BASE SEEDS [--env K=V] [--pre CMD] [--pattern RE] -- cmd`:
  one directory per seed with `cmd.sh` and the success pattern; appends to a
  jobs.list. Idempotent.
- `run_task.sh`: array task i runs line i of the list under `/usr/bin/time -v`
  (`run.log`, `time.log`), marker `done` only if exit 0 and the pattern is in
  the output; skips finished directories; logs host and attempt in `attempts`.
- `submit.py NAME LIST TIME MEM [--cpus N]`: submits every line that is neither
  done nor pending/running under NAME (so it is both the first submission and
  the resubmission of NODE_FAIL/TIMEOUT losses), within the caps (user total
  < 19,000, own `dis-*` < 11,500), `--nice=100`, `--exclude` from
  `/ptmp/mpp/akarlber/cs-production/bad_nodes` (+ `EXTRA_EXCLUDE`). Lines with
  ≥ 4 attempts are left alone and reported.
- `feeder.sh`: hourly self-renewing Slurm job (`dis-feeder`, never on et nodes)
  that runs `submit.py` for every group in a groups file; stops when all lists
  are done, when `<groups>.stop` exists, or after the deadline (14 Oct).
- `nnlojet_card.py`: runcard per part from `nnlo21/validation/nnlojet_epLJJ_zeus2j.run`.
- `combine_nnlojet.py`: equal-weight seed averages per part, parts summed per
  coefficient, output in the format of `combine_zeus.py --nnlojet`.

## Log

### 6 Oct, 22:15–23:05: builds, tables, checks, pilots
- Builds and tables as above. Cluster state at the start: 22 et nodes down*,
  the proVBFH production with ~1,650 jobs.
- Pilot directories `$P/runs/pilot/<group>/s<seed>` (seeds 1001…), all with
  the frozen binaries; B (ZEUS, mode 2) with `P2BSLICE=1 TECHMIN="1d-10 1d-8"`
  and edge 1e-12 for every nlo31 part; D (mode 3) with `TECHMIN` only.
- **Check 1, ZEUS LO (b0, 8 seeds × 7.5M × 6): 103.28 ± 0.005 pb**, against
  MPP b0 103.29 ± 0.03 and NNLOJET 103.29 ± 0.05. Reproduced.
- **Check 2, ZEUS NLO coefficient with projected slicing (b1 32 + lo 32
  seeds, the MPP statistics)**, total [pb] at τ_cut 2e-3, 1e-3, 5e-4, 2e-4,
  1e-4, 3e-5, 1e-5: 11.03, 10.88, 10.78, 10.70, 10.73, 10.68, 10.49
  (± 0.14 … 0.23), against NNLOJET 10.49 ± 0.07 and the MPP P2B values
  (10.49 × (1 + %LO)): 11.09, 10.85, 10.70, 10.48, –, –, 10.59. Bins: m12
  30–45 6.81 ± 0.12 (NNLOJET 6.76 ± 0.06), p̄T 15–22 5.15 ± 0.07 (5.07 ±
  0.02), Q² 125–250 7.64 ± 0.15 (7.46 ± 0.06) at 1e-4. Reproduced.
- **Check D (pilot: b1 16, lo 16 × 20M × 6, disorder p2b 8 × 40M, inclusive
  4):** χ²/15 bins against disorder = 244, 151, 79, 28.5, 13.8, 10.8, 15.4,
  14.8, 14.7, 13.5 at τ_cut 2e-2 … 1e-5, i.e. agreement for τ_cut ≤ 2e-4 and
  power corrections above, as on the MPP machines (χ² 9.3–15 there). Total
  (inclusive, exact): −42.119 ± 0.019 against disorder −42.138 ± 0.001.
  Pilot errors at 1e-4: ≤ 6.0 ‰ of the LO jet rate (1347 pb) per bin (ours),
  ≤ 4.1 ‰ (disorder, 8 seeds). Combination: `nnlo21/cluster/combine_p2b11.py`.
- Timings on the cluster (wall per job; these are much shorter than the MPP
  thserv timings for the mode-2/3 3+1 parts): b0 0.7 min (MPP ~12), lo
  (3M × 6) 2.9 min (~25), kp 2.6 min (~5), b1 20–25 min (~25), D lo (20M ×
  6) 16 min, D b1 15–18 min, disorder p2b (40M) 16–17 min (~26), disorder
  inclusive 1.3 s. Max RSS ≈ 10 MB (nlo31/sliced21), 170 MB (disorder).
- NNLOJET: per channel/part warmups (OpenMP), see below.

### 7 Oct, 01:20–01:40 (after a pause of the session from ~23:25)
- **D production submitted 6 Oct 23:04** (`$P/runs/prod/D`, groups file
  `$P/runs/prod/groups.txt`, feeder `dis-feeder` hourly): sliced21 b1 mode 3
  seeds 2001–4384 (array 48708057), nlo31 lo mode 3 2001–3984 (48708059),
  disorder -p2b -nnlocoef 2001–2592 (48708061), disorder inclusive -nnlocoef
  2001–2012 (48708062), disorder -lo 2001–2008 (48708063). With the pilots
  (seeds 1001…): b1 2400, lo 2000, disorder p2b 600. Sizing from the pilot
  (optimal split for the worst bin, y 0–0.5): ≈ 0.7 ‰ of the LO jet rate per
  bin at τ_cut 3e-5, ≈ 0.5 ‰ at 1e-4; ≈ 1,400 core-h.
- Node failures overnight: nearly all et nodes went down*/completing again;
  TIMEOUT (2 h limit, 16-min jobs) on et09/11/20/21/19 (hung but up),
  NODE_FAIL on et27/39/40 and a few gt/kt. Lost: b1 538, lo 326, dis 286 tasks;
  resubmitted by the feeder. Since 01:25 all `dis-*` groups exclude
  `et[07-40]` (6th column of the groups file).
- **B: r does NOT reproduce with the mandated settings** (P2BSLICE=1,
  TECHMIN="1d-10 1d-8", P2BDROP default 1, edge 1e-12; 40 seeds × 1M × 6):
  - 10 of 40 seeds collapse (r ≈ 0 at every τ_cut, ~1M of 3.2M events with
    plain weights against ~10k in healthy seeds); e.g. seed 1001: iteration 2
    = −1.26e9 ± 1.26e9, after which VEGAS's grid is destroyed. Others have
    χ²/it 30–690.
  - Healthy seeds are close to the notebook (≈ −2950 at τ_cut 1e-4 against
    −2995 ± 35), but the 40-seed mean is −2870 ± 763 (1e-4), −491 ± 47
    (1e-3; notebook −666 ± 15), 28.5 ± 3.5 (2e-2; notebook 43 ± 5).
  - The MPP production r (`runs/znnlop2b_r2`, which reproduced) used the old
    cut 1e-9 W² and P2BDROP = 0; the TECHMIN runs there were TECHDIFF
    differences, so this combination had not been run as a full r before.
  - B production is on hold. Diagnosis running (`$P/runs/rdiag`, same seeds
    1001–1040): v1 TECHMIN + P2BDROP=0; v2 old cut + P2BDROP=1; v3 old cut +
    P2BDROP=0 (the notebook's configuration).
  - The other B parts look healthy (b2 8, vi 16, kp 8 seeds; small smooth
    errors). Timings: r 1.4–1.9 h per job and vi 1.3–1.5 h (MPP ~45, ~50 min:
    about 2× slower here), b2 34 min (~30).
- **NNLOJET production pilot** (20 jobs per part, ~45 min each, `$P/runs/nnlojet/prod`):
  LO 103.283 ± 0.003 pb (ours 103.28 ± 0.005); NLO coefficient (R + V)
  10.678 ± 0.046 pb, 2.2σ above the MPP NNLOJET value 10.49 ± 0.07 (whole
  channels there, different warmups): to be watched with the production. Seed
  scatter / NNLOJET's quoted error = 0.6–1.2 per part.
- NNLOJET warmups (`$P/runs/nnlojet/warmup/<part>/s1`): LO 1M[5]; R, V 2M[5];
  VV_1–7 1M[5] (8 threads); RV_1–7, RRa_1–5, RRb_1–5 10M[6] (32 threads).
  Per-point CPU: LO 7.6 µs, V 0.10 ms, R 0.49 ms, VV 0.15–2 ms, RV 0.6–3 ms,
  RRb ≈ 4 ms, RRa ≈ 27–35 ms (one 10M iteration takes 2.1–2.8 h on 32
  threads; 6 iterations would need 13–17 h). RR iterations still fluctuate
  strongly after 6 iterations (e.g. RRb_3: −1492, 2433, 8768, 3479, 3192,
  4924 fb).
- Event shapes: runcards prepared in `nnlo21/cluster/eventshapes/` (not run).
  NNLOJET's DIS shapes are evaluated in the Breit-frame current hemisphere and
  normalised to Σ|p| there: `dis_thrust` = 1 − Σ|p_z|/Σ|p| (H1's τ_z type),
  **not** our τ_zQ = 1 − 2Σp_z/Q, contrary to the instructions (section 1.C).
  `y23` aborts (SIGABRT) at the first event in epLJJ. `dis_thrust` and
  `dis_C` run (tiny LO tests). Definitions must be matched before a pilot.

### 7 Oct, 02:40–07:30 (session paused ~02:45–07:30)
- r diagnosis (`$P/runs/rdiag`, seeds 1001–1040): with TECHMIN + P2BDROP=0
  (v1) all seeds are stable in the first iterations; the blow-ups (17 of 40
  seeds with an iteration of 1e5–1e9 under the default P2BDROP=1) come from
  P2BDROP=1 (plain-slicing weights for events with τ₂(real) < techcut or a
  degenerate dipole) together with the lower TECHMIN cut. Old cut 1e-9 W² (v2,
  v3) is stable with either drop rule. Final numbers below.
- **A (NNLOJET) production submitted 02:40** (`$P/runs/nnlojet/prod`, list
  `prod.list`, array name `dis-nj-prod`, first chunk array 48715220, rest by the
  feeder): per part (seeds incl. the 20-job pilots) LO 385, R 2066, V 100,
  VV_1–7 100 each, RV_1–7 200 each, RRb_1 700, RRb_2 700, RRb_3 2000, RRb_4 700;
  gated on their warmups (groups-file NEEDS column): RRb_5 700 (`dis-nj-RRb5`,
  released), RRa_1–5 500 each (`dis-nj-RRa<k>`, ncall 120k ≈ 1 h; waiting for
  the warmups, 10M[6] at ~2.1–2.8 h per iteration on 32 threads).
  ncall per part in `$P/runs/nnlojet/prod_ncall.txt` (~45–65 min per job).
  Sizing (`nnlo21/cluster/size_nnlojet.py`, targets LO 1e-4, NLO 0.3 % total /
  1 % per bin, NNLO 2 % of an assumed 40 pb / 5 % per bin): ~12k core-h
  without RRa. RRb_3 is heavy-tailed: one pilot seed of 20 gave 99,983 ±
  97,824 fb (others ±1–8k), so its per-job scatter (22 pb) is dominated by
  rare huge weights; NNLOJET's own combination (dokan) trims such outliers
  (MAD 3.5, ≤ 1 %); ours uses plain equal weights (rule), trimmed values only
  as a diagnostic.
- D production complete (07:07): b1 2383/2384, everything else done.

### 7 Oct, 07:30–08:00: AK's correction of the instructions (fefa1bc) and what changed
**The correction.** AK pushed fefa1bc (merged here as 014c275): production must
use the DEFAULT technical cut (1e-9 W², no `TECHMIN`); the lowered cut
min(1e-10 W², 1e-8 Q²) enters only through a separate difference run
`TECHDIFF=1 TECHMIN="1d-10 1d-8" VEGAS_EQUAL=2` (same arguments) whose seed
average is added with errors in quadrature. Reason (MPP, `runs/tcut11`, mode
1): direct TECHMIN runs are contaminated by rare garbage events. Also for D
("add a TECHDIFF correction run"). Cross-checks redefined (cut convergence via
TECHDIFF + TECHREF; plain slicing with its own TECHDIFF; fixed-point test with
high statistics + TECHDIFF). `P2BDROP` stays at its default.

**How it relates to what I found (01:30–07:30).** Independently, my B pilot r
with direct TECHMIN (default P2BDROP) did not reproduce: 17 of 40 seeds with
iterations of 1e5–1e9, 10 seeds collapsed to r ≈ 0. Diagnosis on the same 40
seeds (`$P/runs/rdiag`), total r [pb]:

| variant | 2e-2 | 1e-3 | 1e-4 | 3e-5 | 1e-5 |
|---|---|---|---|---|---|
| pilot: TECHMIN direct, P2BDROP 1 | 29 ± 3 | −491 ± 47 | −2869 ± 763 | −5802 ± 1396 | −9855 ± 2810 |
| v1: TECHMIN direct, P2BDROP 0 | 33 ± 2 | −657 ± 4 | −2951 ± 7 | −5336 ± 9 | −8501 ± 11 |
| v2: default cut, P2BDROP 1 (= corrected production) | 35 ± 2 | −664 ± 5 | −2966 ± 8 | −5353 ± 10 | −8548 ± 12 |
| v3: default cut, P2BDROP 0 (MPP r production) | 37 ± 2 | −662 ± 5 | −2966 ± 9 | −5350 ± 10 | −8545 ± 13 |
| notebook (300 seeds, MPP) | 43 ± 5 | −666 ± 15 | −2995 ± 35 | – | −8560 ± 44 |

- **Check 3 (short r run against the notebook) passes** with the default cut
  (v2, v3); the per-seed scatter here is ~10× smaller than the MPP 300-seed
  error suggests (the MPP set probably contained tail seeds).
- Paired v1 − v3 (TECHMIN effect): +15.8 ± 9.6 (1e-4), +13.5 ± 12.3 (3e-5),
  +44.7 ± 14.8 (1e-5), consistent with the MPP P′ correction (≈ +11 and ≈ +47).
  Paired v2 − v3 (drop rule, default cut): ≤ 4.6 ± 2.1 everywhere.
- So with P2B (mode 2) the blow-ups of direct TECHMIN come through
  P2BDROP = 1 (plain weights for events with τ₂(real) < techcut or a
  degenerate dipole); with P2BDROP = 0 the 40 seeds were clean. This does not
  contradict AK's mode-1 finding (no P2B there, so a different path); the
  corrected procedure (default cut + TECHDIFF) avoids both.
- **Correction of my 01:40 plan** ("run B with TECHMIN and P2BDROP = 0"): dropped
  in favour of AK's procedure. Nothing of B had been submitted.

**Jobs checked against the corrected text.**
- Nothing pending used the old settings (the pending dis-D-b1 task is sliced21,
  unaffected). Nothing cancelled.
- D: `dis-D-lo` (nlo31 lo mode 3, 2000 seeds incl. pilot) ran with TECHMIN
  directly, all complete. Outlier check: max |mean − 1 %-trimmed mean| = 0.77σ
  over all bins at τ_cut 1e-4…1e-5 (clean). **Kept as a cross-check only, not
  as the result**; its CPU (2000 × ~16 min ≈ 540 core-h) counts as spent on the
  superseded procedure. Retired from the groups file (copy of the old file:
  `groups.txt.v1-0740`). sliced21 b1 and the disorder runs are unaffected.
- B pilots lo/vi/kp/r (TECHMIN direct) are checks only; the direct-TECHMIN r
  pilot (40 × ~1.6 h ≈ 64 core-h) is unusable. v2 (40 seeds) = the corrected
  r production setting and is included in it. Pilot b0/b1/b2 (sliced21) are
  included in the production.

**Submitted 07:55 (corrected; feeder continues hourly), `$P/runs/prod`:**
- D: `lodef` nlo31 lo mode 3, default cut, the same seeds as the TECHMIN set
  (1001–1016, 2001–3984; paired comparison) — 2000; `locorr` TECHDIFF — 400.
- B (P2BSLICE=1, default cut): b1 +168 (200 with pilot), b2 +92 (100), lo 200,
  vi 250, kp 200, r 3460 (+40 v2 = 3500), rcorr (TECHDIFF) 400. Sizing
  (`size_b.py`, pilot per-seed scatter, optimal split): ±1 pb on the total at
  τ_cut 3e-5 needs b2 67, vi 193, kp 168, r 3430 ≈ 5.5k core-h (±1 pb at 1e-4:
  ≈ 3k core-h). The instructions' "about 300 r seeds" gives ≈ ±3.4 pb; AK's
  target is ≈ 1 pb, hence 3500.
- B cross-checks (Bx): cut convergence TECHDIFF TECHMIN="1d-11 1d-9"
  TECHREF="1d-10 1d-8" r 200; psmc edge 1e-13 r 300, vi 64; plain slicing
  (no P2BSLICE) lo 64, vi 64, kp 16, r 300, rcorr 100.
- Fixed point (F; mode 1, x = 0.01, Q² = 400): b1 64, b2 64, lo 64, vi 120,
  kp 32, r 1200, rcorr 300 (≈ 10× the MPP r statistics).
- Every part: mean, median and trimmed mean are compared before combining.
- 08:00: RV_3's 180 production jobs had a malformed command (no ncall: RV_3 had
  no pilot entry in `prod_ncall.txt`, so the seed went into the ncall slot);
  they failed in 0.04 s each. Fixed in place (ncall 430000 from the warmup's
  6.3 ms per point); they are requeued (attempt counter 1 of 4). All other
  NNLOJET commands checked (6 words each).
- Interim A (07:50, ~7,000 jobs): LO 103.2810 ± 0.0005 pb; NLO coefficient
  10.579 ± 0.006 pb (R 1318 + V 100 jobs), 1.3σ from the MPP NNLOJET 10.49 ±
  0.07; seed scatter = NNLOJET's quoted errors (ratio 0.8–1.1) for every part.
- `nnlo21/cluster/check_parts.py`: mean/median/trimmed check per bin. D b1
  (sliced21 mode 3, not affected by the cut question) is flagged: rare seeds
  with z = −100 … −630 in single bins (e.g. seed 2459, y −0.5…0, τ_cut 1e-4:
  −7592 against a median of 422); mean vs trimmed mean up to 1.8σ. D lo (both
  cuts) and v1/v2 r are clean.
- Cluster saturated at 07:50 (all et nodes down*/completing, ~52 ct/gt/kt
  nodes usable); B/F/D-lodef wait behind the NNLOJET jobs (FIFO, same nice).
- 08:20: RRa warmups are slower than estimated: 3.2 h per 10M iteration
  (RRa_1/4/5) and 4.8 h (RRa_2/3) on 32 threads (37–55 ms per point), so
  6 iterations would take 19–29 h (RRa_2/3 beyond the 24 h limit, and then no
  `done` marker, so their production would never be released). Added
  `nnlo21/cluster/freeze_warmup.sh`, run hourly by the feeder through
  `groups.txt.hook`: once a warmup has written 4 iterations (or its job ended
  with ≥ 1 grid) it cancels that warmup task by id (48707779_8 … _12), waits
  until it has left the queue, and writes `done` + `frozen`. Expected release:
  RRa_1/4/5 ≈ 11:30, RRa_2/3 ≈ 17:30. RRa production jobs (120k points) will
  take ~1.2–1.8 h each.

### 7 Oct, 08:20–11:30 (session paused ~09:05–11:20)
- **Cluster-wide failure 10:18–10:20:** every running job of mine was killed
  (state FAILED, exit 1, after ~43 min for the B r tasks), and every task that
  started in those minutes FAILED at launch with exit 0:53 (no output). Other
  users were hit the same way (sqzhang 18,260 FAILED, jmhenn 13; proVBFH
  hxswg136 128). Mine: B r 3,458, vi 250, rcorr 400, F 1,800+, Bx ~1,000,
  D lodef/locorr ~340, the 5 RRa warmups. The feeder resubmitted everything
  at its next run; all are running again.
- The RRa warmups (48707779_8…_12) died in that event after 2–3 iterations
  (RRa_3: 2, others 3). The hook then froze them (their grids were last
  written at 06:11–08:20, well before the kill, so they are intact) and
  released the RRa production (arrays 48755142, 48755244–48755247, 500 each).
  Caveat: fewer warmup iterations than planned (4), so the RRa grids are less
  adapted and the per-job variance is presumably larger.
- Done by 11:20: D lodef 1714/2000, locorr 348/400; B b1 165/168, b2 92, kp 200,
  lo 194/200; F lo 64, kp 31; NNLOJET RRb_5 592/700, all other A parts
  complete except RV_3 (rerunning).

### 7 Oct, 11:20–12:00
- **TECHMIN/TECHDIFF act only on nlo31 part r.** The technical cut s_min <
  min(cw W², cq Q²) and the TECHDIFF weight live in `real_part` only
  (`dis31/nlo31.f90` l. 576–590); lo, vi, kp ignore both. Checked: D lo with
  the default cut (`lodef`) is bit-identical to the direct-TECHMIN D lo for
  all 1,714 common seeds (1,698 production + 16 pilot). Consequences:
  - **Correction of my 07:50 entry:** the TECHMIN-direct D lo set was not
    "superseded" and no CPU was wasted on it; it *is* the default-cut result
    and is the D result. The `lodef` set (≈ 540 core-h) is a duplicate
    (kept, it confirms the identity).
  - The D "TECHDIFF correction" is zero by construction (D's nlo31 part is lo
    only). My `locorr` runs integrated the full lo integrand with equal-weight
    iterations, not a difference: invalid. 348 finished (≈ 95 core-h wasted);
    the remaining 52 running tasks were cancelled by job id (array 48754221),
    and the group is retired in `groups.txt`.
  - Likewise the B pilots lo/vi/kp with TECHMIN equal the default-cut runs;
    only the r pilot differs (and was the unstable one).
- Outage bookkeeping: 3,037 attempt entries of tasks killed 10:15–10:25
  removed from the `attempts` files (copies in `attempts.outage`; list of
  task ids in `$P/runs/outage_tasks.txt`); no line has more than 2 attempts.
- Crash: B r seed 4711 (48754427_2711, gt37) died after iteration 1 with
  "double free or corruption (out)" (the only one of 3,322 B r logs).
  Rerun with the same seed on gt37 (48759484) and elsewhere (48759485),
  `$P/tests/crash4711`.
- **D final** (`$P/runs/D_final.txt`, `.json`): b1 2400, lo 2000, disorder p2b
  600, inclusive 16, LO jet rate 1346.36 pb. χ²/15 against disorder: 24.2,
  18.0, 15.8, 16.3 at τ_cut 2e-4, 1e-4, 3e-5, 1e-5 (largest pull 2.8σ, y
  1.5–2.5 at 1e-4); above 2e-4 power corrections (χ² 287 at 1e-3). Inclusive
  total −42.1305 ± 0.0067 against −42.1368 ± 0.0001. Errors (ours) ≤ 0.44 ‰
  of the LO jet rate per bin at 1e-4 (0.50 ‰ at 3e-5), except the two y bins
  −1…−0.5 and −0.5…0 (2.5 ‰), which are dominated by a few b1 outlier seeds.
- **b1 trimmed vs plain (AK asked; diagnostic only, no seeds dropped):**
  1 % trimmed mean (24 of 2400 seeds per side per cell) minus plain mean,
  ‰ of the LO jet rate, at τ_cut 1e-4: y −1…−0.5 −2.80 (−1.1σ), y −0.5…0
  +2.83 (+1.1σ), p_T 20–30 −0.21 (−0.6σ), p_T 30–50 +0.10 (+1.0σ), all other
  bins ≤ 0.11 ‰ (≤ 0.3σ); same at 3e-5 and 1e-5. The outliers move weight
  between the two neighbouring y bins (events whose O(event) and O(Born) fall
  on either side of y = −0.5). With trimming the largest error drops to
  0.32 ‰ (1e-4); χ²/15 16.9, 12.8, 13.6 at 1e-4, 3e-5, 1e-5.
- The disorder reference is also heavy-tailed (per-seed scatter 44 pb in p_T
  8–15 against a robust 12 pb; trimmed = plain within 1σ) and reached only
  1.34 ‰ in its worst bin with 600 seeds. Added 1,000 seeds (2593–3592,
  array 48759492, ≈ 280 core-h) → ≈ 0.8 ‰ expected.
- **AK's decision (7 Oct ~11:50):** plain equal-weight means are the result for
  every part (A, B, D, F); trimmed means and medians are reported only as
  diagnostics alongside, until the comparison is ready. No seeds are dropped.
- **Crash of B r seed 4711 (reproducible, diagnosed; not fixed during the
  production).** Same seed reruns on gt37 and ct30 crash identically after
  iteration 1 (−10044.6050 ± 1374.8). A `-O2 -g -fbacktrace` build
  (`$P/builds/nlo31-g`) gives: abort in `free` returning from `kt_jets`,
  called from `zeus_jets` ← `zeus_bins` ← `accept_zeus` ← `real_evalv`
  (nlo31.f90:539, a mapped dipole configuration) ← `real_part`. Cause:
  in `kt_jets`, if a momentum is NaN (numerically degenerate mapped Born),
  every comparison with `dmin` is false, so `ii` stays 0 with `beam = .true.`
  and `jets(:,nj) = q(:,act(0)); act(0) = act(m)` writes before the automatic
  array `act` → heap-header corruption → "double free or corruption (out)".
  Proposed fix (for AK, on thA371a): in `kt_jets` (and `zeus_jets`), treat
  `ii == 0` / non-finite momenta as a failed event (pass = .false. or
  sg = NaN, which `real_part` already drops). Corrupting the malloc header
  aborts the program, so other jobs' results are not silently affected; the
  seed is left out (the feeder stops after 4 attempts). Effect of losing one
  of 3,500 r seeds: negligible.
- Also found with a bounds-checked build (`$P/builds/nlo31-bc`): line 625
  `hcacc = hcacc + sg*w*wgt` adds `sg(nv)` (nv = 150 in mode 2) to
  `hcacc(ncell = 160)` at every point. In the optimised build this reads 10
  values beyond `sg` and writes only the unused cells 151–160: harmless for the
  results, but a latent bug (and it prevents bounds-checked runs in mode 2).

## 7 Oct, 16:20–17:00: final state and results
All production lists are complete except 5 B r seeds, 1 B rcorr, 8 Bx conv,
2–3 Bx edge-r seeds that crash deterministically ("double free", the
kt_jets bug above; ≥ 4 attempts). Feeder stopped (`groups.txt.stop`, pending
feeder 48766335 cancelled). Combined numbers in `nnlo21/cluster/results/`
(raw outputs stay in `$P/runs`). **Plain equal-weight means are the results
(AK); everything marked DIAG is a diagnostic only.**

### A: NNLOJET (`results/nnlojet/`, `parts.txt`)
11,931 production jobs (+ warmups). LO 103.2810 ± 0.0005 pb; NLO coefficient
10.5797 ± 0.0048 pb; **NNLO coefficient 11.12 ± 1.19 pb** (target 2 % not
reached: ±1.19 is 11 %; dominated by RRa_3 ±0.80 and RRb_3 ±0.63 pb, heavy-
tailed RR parts; RRa warmups frozen after 2–3 iterations). Seed scatter =
NNLOJET's quoted errors for all parts (ratio 0.8–1.1).

### B: ours, ZEUS, P2B slicing, default cut + TECHDIFF (`results/B_*`)
- NLO (b1 200 + lo 200): 10.48 ± 0.08, 10.45 ± 0.09, 10.48 ± 0.10 pb at τ_cut
  1e-4, 3e-5, 1e-5; NNLOJET 10.580 ± 0.005 (pulls −1.3, −1.4, −0.9).
- **NNLO, plain means: not usable below τ_cut 5e-4.** r (3,495 seeds) and its
  TECHDIFF correction rcorr (399) contain blow-up seeds: r has 41 seeds with a
  VEGAS iteration > 1e6 (median 9e3) and seeds with −inf cells at 3e-5/1e-5;
  rcorr has 100 seeds with an iteration > 1e4 (median 1e2), values up to 1e163.
  Plain NNLO total: 41.8 ± 2.1 (5e-4), then 1.0e4, 4.3e5, −inf, −inf.
  (The 40-seed tests — v2, and AK's 60-seed TECHDIFF runs — were too small to
  show this.)
- DIAG (r: 64 of 3495 seeds, rcorr: 100 of 399 left out for non-finite cells or
  an iteration > 1e5 / 1e4): NNLO total 18.1 ± 0.6 (2e-3), 24.0 ± 0.8 (1e-3),
  33.2 ± 1.7 (5e-4), 40.7 ± 1.8 (2e-4), 48.1 ± 2.0 (1e-4), 49.3 ± 5.8 (3e-5),
  50.1 ± 6.0 (1e-5) pb.
- **Against NNLOJET (11.1 ± 1.2 pb): a factor ~4, 6–16σ at τ_cut ≤ 5e-4**
  (DIAG numbers; the plain ones are undefined there). Our values still rise
  between 2e-4 and 1e-4. **This contradicts the MPP statement (6 Oct) that
  projected slicing is flat from 2e-4 (40.1 … 43.2 pb)**, and it puts the
  whole NNLO 2+1 result in mode 2 into question (our side or the comparison).
  Only at τ_cut 2e-2 … 1e-2 (11.3, 11.8) is ours close to NNLOJET, which there
  is presumably a coincidence of power corrections.
- Cross-checks (paired seeds where possible; DIAG = blow-up seeds left out):
  - psmc edge 1e-13 − 1e-12, r: +10.8 ± 3.4 (1e-3), +22.0 ± 6.7 (1e-4),
    +27 ± 9 (3e-5), +44 ± 11 (1e-5) (DIAG; plain: +9.5 ± 5.1 at 1e-3, +20.9 ±
    18.7 at 1e-4). **Contradicts the MPP finding "r unchanged with edge
    1e-13"**; vi: +0.23 ± 0.25 (1e-4), unchanged.
  - cut convergence (TECHDIFF 1e-11/1e-9 against 1e-10/1e-8): 111 of 192 seeds
    blow up; DIAG rest +4 ± 5 (1e-4), +1 ± 21 (3e-5), +20 ± 30 (1e-5).
  - plain slicing: r(plain) − r(P2B), paired, DIAG: −8.7 ± 2.7 (1e-3), −9.7 ±
    4.5 (1e-4); plain rcorr DIAG +3.4 ± 1.5 (1e-4), +44.5 ± 4.9 (1e-5).

### F: fixed point x = 0.01, Q² = 400 (`results/F_*`)
- NLO: 16.68 ± 0.05, 16.82 ± 0.07, 16.80 ± 0.08 at 1e-4, 3e-5, 1e-5 against
  DISENT 16.885 ± 0.04.
- NNLO plain: 11.91 ± 0.26 (2e-3), 11.05 ± 0.34 (1e-3), 10.46 ± 0.44 (5e-4),
  10.33 ± 0.53 (2e-4), 11.47 ± 0.75 (1e-4); −74 ± 66 (3e-5) and −2608 ± 3100
  (1e-5) from rcorr blow-ups. DIAG (3 r, 9 rcorr seeds left out): 11.19 ± 0.74
  (1e-4), 11.7 ± 5.5 (3e-5), 7.1 ± 6.4 (1e-5); without the correction 11.1 ±
  0.8 (3e-5), 7.3 ± 1.0 (1e-5). The plateau 5e-4 … 1e-4 is ≈ 10.3–11.5,
  **≈ 2 above the MPP plateau 8.5 ± 1.2** (1.5σ of the MPP error); the
  TECHDIFF correction is too noisy to resolve 3e-5 and 1e-5.

### D: NNLO 1+1 by P2B (`results/D_*`)
b1 2400, lo 2000, disorder p2b 1600, inclusive 16. χ²/15 against disorder:
33.0, 25.3, 22.6, 22.4 at τ_cut 2e-4, 1e-4, 3e-5, 1e-5 (largest pull 2.8σ),
power corrections above. Inclusive total −42.1305 ± 0.0067 vs −42.1368 ±
0.0001. Errors: ours ≤ 0.44 ‰ of the LO jet rate (1346.4 pb) per bin at 1e-4
(0.50 ‰ at 3e-5) except y −1…−0.5 and −0.5…0 (2.5 ‰, b1 outlier seeds);
disorder ≤ 0.61 ‰. b1 trimmed-vs-plain (DIAG) as in the 11:20 entry.

### CPU (user time from the time.log files of the last attempts)
A: NNLOJET production 12,806 core-h, warmups 782 (+ ≈ 1,900 for the five RRa
warmups killed at 10:18, 32 threads × 11.7 h, not in time.log); B 5,570;
F 1,606; Bx 1,135; D 2,556 (incl. the duplicate lodef ≈ 540 and the invalid
locorr ≈ 100); pilots 113; r diagnosis 166. Attempts lost to the 10:18 outage
and to node failures are not included (several hundred core-h).

## 7 Oct, 21:00–: UPDATE 7 Oct evening (1a46f6a): r with equal iteration weights
- Merged origin/2026-10-nnlo21 (b4ca068, 1682b6d, 1a46f6a; fast-forward, my
  commits and e1884e7 included). Cause of the factor 4 (thA371a): nlo31's
  1/σ²-weighting of iterations biases r (+49.8 ± 7.1 pb at τ_cut 1e-5, MPP).
- New frozen binaries `$P/bin2/{sliced21,nlo31}` (1a46f6a; md5 in
  `$P/bin2/MD5SUMS`); old `$P/bin` kept. Check: b0, lo and r (mode 2, small
  runs) bit-identical between old and new binaries.
- **Step 2, bias from the existing outputs** (`nnlo21/cluster/iter_bias.py`,
  `$P/runs/iter_bias.txt`): per seed, equal-weight average of iterations
  2–6 (target cell) minus the reported cell. Plain means are undefined (garbage
  seeds), so DIAG statistics: B r reported − equal = +39.6 ± 24.5 pb (3,440
  unflagged seeds; trimmed +42.4, median +16.5); edge r +46 ± 13; plain r
  +5.8 ± 21 (trimmed +29); F r +9.45 ± 0.62 pb/GeV² (mode 1, significant).
  Consistent with AK's +50 pb.
- **Step 3, reruns with VEGAS_EQUAL=2 and bin2, same seeds** (paired with the
  old runs), groups file `$P/runs/prod/groups2.txt`, feeder started 21:35
  (48777071): B r-eq 3,500 (array 48777072), F r-eq 1,200 (48777074), Bx
  edge-r-eq 300 (48777075), Bx plain-r-eq 300 (48777076). TECHDIFF runs
  (rcorr, plain-rcorr, F rcorr, conv) kept as they are (already
  VEGAS_EQUAL=2); b0/b1/b2/lo/vi/kp and NNLOJET not rerun.
- Step 4 replays with P2BDEBUG=1 (bin2, old weighting, same settings):
  B r seed 3197 (iteration 2 = −1.2e162, then 0) job 48777079, seed 4637
  (iteration 4 = −9.4e102, then ≈ 0) job 48777080; `$P/runs/dbg/s*/dbg.log`.
- 8 Oct 02:30: reruns complete (B r-eq 3500/3500, F r-eq 1199/1200, edge-r-eq
  300, plain-r-eq 300); seed 4711 completed with bin2 (kt_jets fix). Early
  look (662 seeds, 22:56): with equal weights the same seeds still blow up
  (2042, 2382, 2592, 2754, 3197: iterations up to 1e162), so the plain r mean
  stays undefined at 3e-5/1e-5; median/trimmed shift as expected (−42/−60 pb at
  1e-5 against the old runs).
- 8 Oct 03:00: results of the equal-weight reruns in `REPORT-2026-10-08.md`
  and `results/` (B_eq_*, F_eq_*, B_eq_bins.txt, outliers_r.txt,
  iter_bias.txt, dbg/). Summary: with VEGAS_EQUAL=2 the ZEUS factor 4 is gone
  (plain NNLO 11.8–14.5 pb from τ_cut 2e-2 to 2e-3 against NNLOJET 11.1 ± 1.2;
  DIAG within 2.1σ down to 1e-5, errors 3–25 pb below 5e-4 because of the r
  tail). **Correction of my 7 Oct statements:** the factor 4, the "edge
  shift" (+22 pb, now median +8) and the fixed-point plateau (10.3–11.5) were
  artefacts of the 1/σ² iteration weighting; with equal weights the fixed
  point falls from 10.1 (1e-3) to 6.1 (1e-4) pb/GeV² (DIAG). Plain B means
  stay undefined below 1e-3 (59 r and 102 rcorr seeds with blow-ups; DBG
  replays: mapped-Born partons of ~5e-4 GeV give dipoles of 1e120–1e179).
  Feeder for groups2 stopped (`groups2.txt.stop`).

## 8 Oct, ~03:30: UPDATE 8 Oct (d638917, DIPGARB) — prepared, NOT submitted (AK wants the CPU estimate first)
- Merged origin/2026-10-nnlo21 at d638917. Built `$P/bin3/{sliced21,nlo31}`
  (bin and bin2 kept). Check (mode 2, small runs, VEGAS_EQUAL=2 for nlo31):
  b0 and lo bit-identical to bin2; r identical apart from the new line
  "r: events dropped (garbage dipole) 0 of 2776".
- Job lists (bin3, VEGAS_EQUAL=2, the same seeds and otherwise identical
  commands; each cmd.sh compared with its predecessor: 0 differences besides
  the binary path): `$P/runs/prod/{B/r-g 3500, F/r-g 1200, B/rcorr-g 400,
  F/rcorr-g 300, Bx/plain-rcorr-g 100, Bx/conv-g 200, Bx/edge-r-g 300,
  Bx/plain-r-g 300}`; groups file `$P/runs/prod/groups3.txt` (feeder not started).
- CPU estimate from the measured per-job CPU of the same sets: B r 1.42 h/job
  → 4,990 core-h; F r 1.09 → 1,310; B rcorr 0.50 → 200; F rcorr 0.32 → 100;
  plain rcorr 0.35 → 35; conv 0.60 → 120; edge r 1.24 → 370; plain r 1.13 →
  340. **Total ≈ 7,460 core-h**, 6,300 jobs (longest single job 2.2 h).
  Wall clock: ≈ 4–5 h with the ~1,900 concurrent slots of the last round,
  ≈ 3 h if ~4,000 slots are free (4,121 idle CPUs on alma at 03:30).
- To start after the go-ahead:
  `GROUPFILE=$P/runs/prod/groups3.txt nnlo21/cluster/feeder.sh`
- **8 Oct 08:49: submitted after AK's go-ahead** (feeder on groups3.txt,
  48788281): B r-g 48788282, F r-g 48788350, B rcorr-g 48788351, F rcorr-g
  48788352, plain rcorr-g 48788353, conv-g 48788354, edge r-g 48788355,
  plain r-g 48788356 (6,300 jobs, ≈ 7,460 core-h).

## 8 Oct, 14:40-20:30: results of the DIPGARB reruns (bin3, VEGAS_EQUAL=2) -- `REPORT-2026-10-08b.md`
- All 8 groups complete (B r-g 3500, F r-g 1200, B rcorr-g 400, F rcorr-g 300, plain rcorr-g 100,
  conv-g 200, edge r-g 300, plain r-g 300), no failed seeds. Feeder 48799109 cancelled (nothing left).
- Tools: `garb_report.py` (before=bin2 dirs, after=bin3 dirs; `results/garb_report_8b.txt`), `combine_b.py`
  (`results/{B,F}_g_{plain,diag,trim}.*`; DIAG limits r 1.79e5, rcorr 1.99e3 (B), 6.79e4/318 (F)).
- Seeds with an iteration > 1e30: B r 9 -> 0, B rcorr 7 -> 0 (F: 0). Garbage-dipole drops ~5e-4 of events (B r mean
  1677/seed). 20x-rule flagged: B r 59 -> 52, B rcorr 102 -> 96, F r 3 -> 3, F rcorr 13 -> 13. 35 of 3500 B r seeds
  change their max iteration at all.
- **Contradicts the expectation in the 8 Oct instructions** that the heavy tail would go away: the guard removes the
  1e60-1e179 garbage but the 1e5-1e6 tail (r) and the rcorr tail remain; plain means below 1e-3 are still
  dominated by single seeds (B 1e-4: 4e5 +- 4e5; 1e-5: -1e7 +- 1.3e7). Robust/sample sigma essentially unchanged.
- Replays with bin3, P2BDEBUG=1, VEGAS_EQUAL=2 (seeds 2042, 2049; reproduce the same iterations 3 / 2; stopped
  after the DBG lines were in; `results/dbg/g*`): the big iterations come from events with real m4 and t2 ~ 1e-8..1e-11
  GeV^2 (smin/W2 ~ 1e-9) where real and a dipole are 1e13-1e18 with residual 1e10-1e12 after subtraction (incomplete
  cancellation, dipole values far below the 1e30 guard). Not garbage dipoles; a guard on |D|>1e30 cannot catch them.
- Tables: see REPORT-2026-10-08b.md. B plain pulls vs NNLOJET +0.4..+1.8 down to 5e-4; DIAG within 2.0 sigma everywhere.
  F unchanged from the bin2 equal-weight result (DIAG 11.6 -> 6.2 between 2e-3 and 1e-4, not flat).
- Note: b1/b2 now taken from the 168/92 done prod seeds (previous tables: 200/100); NLO row shifts by 0.01.

## 9 Oct, 14:00-14:30: UPDATE 9 Oct (f29517e, psmc quad precision) -- bin4 built, checks, submitted
- Merged origin/2026-10-nnlo21 (f29517e) into 2026-10-nnlo21-cluster (fast-forward). Built on a compute node
  (job 48804374, `nnlo21/build_sliced21.sh`, same as bin3) in `$P/builds/sliced21-f29517e`; frozen
  `$P/bin4/{sliced21,nlo31}` (md5 in `$P/bin4/MD5SUMS`); bin, bin2, bin3 untouched.
- Checks (`$P/runs/chk4`, 30k points, 3 iterations, seed 7, VEGAS_EQUAL=2): sliced21 b0, b1, b2 outputs bit-identical
  bin3 vs bin4 (cmp). nlo31 lo (300k points x 6, seed 7): bin3 1644.22 +- 4.33 pb, bin4 1642.37 +- 4.31 pb
  (difference 1.9 pb = 0.3 sigma; the files differ point by point as expected: rotated transverse basis).
- Job lists (bin4, VEGAS_EQUAL=2, same seeds): `$P/runs/prod/{B/r-4 3500, B/rcorr-4 400, Bx/conv-4 200,
  Bx/edge-r-4 300, Bx/plain-r-4 300, Bx/plain-rcorr-4 100}`, made by copying the bin3 *-g cmd.sh with
  bin3 -> bin4; diff of every cmd.sh against its bin3 predecessor: exactly one changed line (the binary path), 0 others.
  No F runs. Groups file `$P/runs/prod/groups4.txt`.
- CPU estimate from bin3 per-job CPU: 3500x1.42 + 400x0.50 + 200x0.60 + 300x1.24 + 300x1.13 + 100x0.35
  = 4,970+200+120+372+339+35 = 6,036 core-h (target 6,050; within 20%, pre-approved).
- Submitted 14:28 (feeder 48804389 on groups4.txt): B r-4 48804390, B rcorr-4 48804488, Bx conv-4 48804489,
  edge r-4 48804490, plain r-4 48804491, plain rcorr-4 48804492 (4,800 jobs).

## 9-10 Oct: results of the bin4 (f29517e) reruns -- `REPORT-2026-10-09.md`
- Done: B r 3500, rcorr 400, conv 200, edge r 300, plain r 300, plain rcorr 100. 104 B r tasks of the first array never
  started (feeder resubmitted them 17:28); task 48811454_1664 hung RUNNING 3 h without starting the script (kt08), cancelled
  by id, s3624 rerun (48812664). An accidental duplicate submission 48812536 was cancelled by id. Feeder 48804389 chain
  stopped (stop file, feeder 48812262 cancelled).
- **Correction of the 8 Oct (8b) reading: the 1e5-1e6 tail was not incomplete cancellation in deep corners but accepted
  dipoles with negative-energy mapped-Born partons (psmc precision), see notebook 9 Oct.** bin4 confirms this: flagged B r
  seeds 52 -> 0, max |iteration| 1e6 -> 2.2e4, sample sigma at 1e-5 2.5e6 -> 173 pb (robust 95); rcorr 96 -> 4 flagged,
  conv 66 -> 2; plain-weights counts per seed 4724 -> 615 (B r); B NNLO plain mean now finite and within 1.7 sigma of
  NNLOJET at all tau_cut down to 1e-5 (error 5.9 pb there).
- Replays s5368 (r), s2077 (rcorr) with P2BDEBUG=1: no DBG BIG events at these iteration sizes (8e3-2e4), `results/dbg/h*.log`.
- Tools: garb_report.py (plain-weights counts added), combine_b4.sh.

## 10 Oct: UPDATE 10 Oct (Z, W, E NNLO validations against NNLOJET) -- phase 1 (build, tables, pilots, submission)
- Code: HEAD bf8d504 (origin/2026-10-nnlo21-ncc merged). bin5 = `$P/bin5/{sliced21,nlo31}` (md5 in `MD5SUMS`), built on a
  compute node (job 48819098) in `$P/builds/sliced21-bf8d504`; bin..bin4 untouched.
- Checks bin4 vs bin5 (`$P/runs/chk5`, VEGAS_EQUAL=2, seed 7): sliced21 b0 (zeus tables old 3-class for bin4, zeus5 for bin5,
  30k x 3) bit-identical. nlo31 lo mode 2 (100k x 3) is NOT bit-identical: iteration 1 is identical, iterations 2, 3 differ by
  <= 2e-3 relative (0.2 sigma), also without P2BSLICE; both binaries are deterministic on their own. Code differences in me31/
  virt31/nlo31 flavour sums (Z, W refactoring) change the integrand at rounding level and the VEGAS grid follows.
- Tables: `$P/tables/zeus5_{0..28}.tab` (photon/Z) and `zeus5w_*` (EW31="2 0"), 29 jobs each (`$P/tables/mk5`). `tabchkg`:
  zeus5 max deviation / max coefficient (up, ubar, down, dbar, gluon) 9.3e-5 4.9e-3 5.6e-4 3.6e-3 2.1e-3; zeus5w 4.8e-4 3.0e-3
  3.1e-4 3.8e-3 2.0e-3 (all < 1e-2).
- NNLOJET: Z, W with the existing OpenMP build (`$P/nnlojet`); E with `$P/src/omp-tauzq/nnlojet-v1.0.2` (fresh v1.0.2 source +
  `nnlojet-v1.0.2-dis_tauzq.patch`, OPENMP=ON, own install; job 48819103); the production install is not touched.
  Set directories `$P/runs/nnj-{Z,W,E}` made by `mk_nnj_set.sh` (template runcard, `nnlojet` link to the install used;
  `nnlojet_card.py` takes the process from the template, `nnlojet_prod.sh` the install from `$N/nnlojet`). Warmups: LO 1M[5],
  R V 2M[5], VV 1M[5], RV RRb 10M[6], RRa 8M[4] (Z RRa_1) / 5M[4] (others: cheaper, A had frozen RRa warmups at 4 iterations).
  E warmups LO, R, V, VV converge (chi2/it ~1, LO 970.3 pb, R 73.0, V 175.4 pb in the last iterations).
- Per-job CPU at the same ncall relative to the photon ZEUS production (A): E 0.8, Z 1.15 (LO, R, V, VV).
- Bug fixed on the way: `nnlojet_prod.sh`, `size_nnlojet.py`, `combine_nnlojet.py` globbed `DIS.*`; W's files are `DISWm.*`.
- Our side (`mk_ours_set.py`, lists `$P/runs/ours5/{E,Z,W}[-pilot]/<part>`, B commands with bin5): Z EW31="1 0", W EW31="2 0",
  P2BSLICE=1; E ZFIX=3 plain (+ P2BSLICE=1 P2BEXTRA=1 lo/vi/kp/r, a quarter of the seeds); b0 16 seeds (7.5M x 6), b1 200, b2 100,
  lo 200, vi 250, kp 200, r 3500 (seeds 1001-), rcorr 400 (TECHDIFF TECHMIN="1d-10 1d-8"), VEGAS_EQUAL=2 for r, rcorr.
  Pilots (4 seeds each, seeds 9001-9004): CPU per job [core-h] b1 0.40/0.44/0.42 (E/Z/W), b2 0.64/0.74/0.72, lo 0.12/0.14/0.14,
  kp 0.085/0.11/0.10, vi 1.01/2.92/3.16, r 1.10/2.89/1.91, rcorr 0.37/0.81/0.55; E extras lo-x 0.11, kp-x 0.075, vi-x 0.22,
  r-x 1.26. All outputs finite.
- CPU estimate (10 Oct, before the NNLOJET pilots exist; core-h). Ours (measured pilots x seeds): E 5.6k (incl. extras 1.1k),
  Z 11.4k, W 7.9k = 24.9k. NNLOJET: production sized as the photon ZEUS one (A, 12.8k) times the per-job CPU ratio (E 0.8,
  Z 1.15, W assumed 1.15) = 10.2k + 14.7k + 14.7k, warmups ~2.5k per set, pilots ~0.5k per set: ~48k. Total ~73k against
  the 70k budget (91k limit). Ours submitted 10 Oct (feeder 48820641, groups `$P/runs/ours5/groups5.txt`, order E, Z, W).
  NNLOJET production NOT yet submitted: waits for the warmups (E, W big ones still pending/running) and the production pilots
  (`$P/runs/nnj-*/prod/<part>/pilot.list`, 5 seeds per part, fed by `nnj_pilot_wait.sh`); sizing then with
  `OBS=dis_tauzq size_nnlojet.py` (E) or the default observables (Z, W).

## 11 Oct: d2 (second technical-cut decade), NNLOJET pilots and production sizing (phase 2)
- d2 (item 7): `mk_ours_set.py SET d2` (E, Z, W): 200 seeds 2001-2200 each (as rcorr/conv), lists `$P/runs/ours5/{E,Z,W}/d2`,
  `TECHDIFF=1 TECHMIN="1d-11 1d-9" TECHREF="1d-10 1d-8" VEGAS_EQUAL=2`; cmd.sh diffed against rcorr's (s2001): only TECHMIN changed
  and TECHREF added. Groups `$P/runs/ours5/groups-d2.txt`, submitted with the feeder (GROUPFILE=that file): feeder 48847583,
  arrays E 48847584, Z 48847585, W 48847586 (200 tasks each). CPU estimate from rcorr per job (0.37/0.81/0.55 core-h):
  E 74, Z 162, W 110 = 0.35k core-h (the lower cut may cost a little more). (First run of the script overwrote the per-set
  `groups.txt` of E/Z/W; restored from groups5.txt. d2 writes `groups-d2.txt` now.)
- NNLOJET pilots (5 seeds per part, 9001-9005): all done except E RV_7 s9003-9005 and W RRb_2 s9001-9005, which never reached
  the queue (nnj_pilot_wait.sh submitted only 2 of 5 tasks for E RV_7 and nothing for W RRb_2); resubmitting them was refused
  by the permission classifier and left to AK; the production covers those parts (W RRb_2 sized from W RRb_1 sigma and CPU,
  E RV_7 from 2 pilots). nnj_pilot_wait.sh never exits because of these (stopped my wait).
- Measured per-job CPU (core-h, production ncall of A) E/Z/W: LO 0.63/0.85/0.73, R 0.64/0.86/0.32, V 0.61/0.87/0.32,
  RRb_1 0.74/0.91/0.62, RRa_1 1.26/1.81/0.73, VV_1 0.16/0.18/0.15. Ratio to A (RRb_1 ~ 1): E 0.8, Z 1.15 as predicted from the
  warmups; **W measured: R/V/RR/RV/VV are 0.3-0.7 of the Z cost (W/Z about 0.4-0.75, mean ~0.55), not 1.15 relative to A**
  (so W is ~0.55 x Z ~ 0.63 x A per job).
- Sizing (`size_nnlojet.py`, now with env NNLO_TOT_FB (absolute NNLO total target, fb), FLOORFRAC (per-bin absolute floor =
  FLOORFRAC x total target; A: 1/6), NLO_FLOOR (NLO bin floor as fraction of the NLO total; A: none), FALLBACK=PART=DONOR,
  MINJOBS) with A's relative LO/NLO targets (LO 1e-4; NLO 0.3 % total / 1 % per bin plus floor 0.3 % of the NLO total per bin).
  NNLO totals: Z 0.8 pb (as A), W 30 fb (2 % of an assumed 1.5 pb = 40 pb x LO ratio), E 5.5 pb. Floors: Z, W 0.4 x T, E 0.5 x T.
  The A rule (floor T/6 for Z, T = 0.8 pb) would cost 21k for Z alone, bin-driven, so floors were relaxed; E cannot follow A's
  relative target (RR sigma/job 30-180 pb at LO 970 pb: 400k core-h at T = 0.8 pb). Sizing outputs `$P/runs/nnj-{E,Z,W}/sizing.txt`.
  Tool: `cluster/mk_nnj_prod.py` (production seeds 1001-, finished pilots count; lists `.../prod/<part>/prod.list`,
  groups `$P/runs/nnj-groups.txt`, 51 groups, 6 h limit, 12 h for RRa_3 of E/Z, 2000MB, exclude et[07-40]). NOT submitted.
- CPU table [core-h] (warmups ~2.5k per set incl. done, pilots ~0.5k per set): ours E 5.6k + 0.07k, Z 11.4k + 0.16k, W 7.9k + 0.11k
  (25.3k with d2); NNLOJET production E 18.6k, Z 7.7k, W 4.3k (30.6k); warmups 7.5k, pilots 1.5k: total ~65k against 70k (limit 91k).
- NNLOJET feeder (to be run by the coordinator, AK approval):
  `GROUPFILE=/ptmp/mpp/akarlber/disorder-nnlo21/runs/nnj-groups.txt DEADLINE=$(date -d 2026-10-16 +%s) /u/akarlber/work/disorder/nnlo21/cluster/feeder.sh`

## 11 Oct evening: /ptmp file quota hit, recovery
- 17:30: the /ptmp file-count hard limit (1,572,864 files) was reached. 3,556 running NNLOJET jobs
  (mostly E RRb_2) died writing their histograms ("Disk quota exceeded", exit 2; ~4,200 core-h lost),
  ~6,000 failed at launch (exit 0:53), and both feeder chains died (17:39, 18:07).
- Recovery (AK OK): finished runs archived to /ptmp/mpp/akarlber/archive/*.tar (cs-production prod-nnlo and
  prod-lonlo, runs/nnlojet/prod, runs/prod/{B,D,F,Bx}); entry counts checked before removing originals.
- run_task.sh now calls nnj_clean.sh on success in NNLOJET production dirs (path contains /prod/): keeps
  done, time.log, *.dat; packs run.log, attempts, cmd.sh, pattern, job.run, NNLOJET log into logs.tgz;
  deletes the copied warmup grid. 11 -> 5 files per job. Already-finished prod dirs cleaned the same way.
- 20:52: both chains restarted (dis-feeder-nnj with nnj-groups.txt, dis-feeder with ours5/groups-d2.txt,
  deadline 16 Oct). Failed jobs had no done marker and 1 attempt, so they were resubmitted.
- combine_nnlojet.py is the bin4 (ZEUS q2/ptavg/m12) combiner: it does not read the E tauzq histograms.

## 12 Oct: bin6, Ex6 and G6 (UPDATE 12 Oct steps 1-4)
- bin6 = `$P/bin6/{sliced21,nlo31,tau2_nlo}` (md5 in MD5SUMS), built from HEAD 7ddc208 (16-core job 48907174, `$P/builds/sliced21-7ddc208`;
  tau2_nlo from `$P/builds/disent`, cmake -DDISORDER_SLICING=ON). Checks `$P/runs/chk6/chk.sh` (VEGAS_EQUAL=2, seed 7): sliced21 b0 30k x 3 and nlo31 lo mode 2
  100k x 3 without VTARGET: outputs bit-identical to bin5 (cmp).
- Generator `mk_g6.py G6W|G6P|Ex6 pilot|prod` -> `$P/runs/ours6/<SET>[-pilot]/<part>`; groups `$P/runs/ours6/groups6.txt` (order G6 W, Ex6, G6 photon),
  feeder chain dis-feeder-g6 (deadline 19 Oct). DISENT photon NLO: `$P/runs/ours6/DISENT-ph` (40 x 3M, seeds 1-40, no -CC, 30 min each, done).
- Pilot CPU per job [core-h] (4 seeds; Ex6 r-x 19-20 seeds): G6 W / photon: b1 0.457/0.453, b2 0.774/0.748, lo 0.136/0.119, vi 3.08/1.37, kp 0.092/0.090,
  r 1.67/1.14, rcorr 0.42/0.34, d2 0.43/0.33; Ex6 r-x 1.43 (others from E pilots 0.11, 0.22, 0.075).
- Estimate (core-h): G6 W 4.45k, G6 photon 2.91k, Ex6 2.96k (r-x 2.87k): total 10.3k vs 9k budget (+14%, within +30%): submitted.
- Ex6 20 r-x pilot seeds (VTARGET=0): flagged (20x median) 1/19 finished at the time (s9013) = 5%, against 86/875 = 9.8% in bin5.
- File count: ~7 files per job dir x 9.0k jobs (incl. 2650 Ex6) = ~65k on top of 800k = ~865k (< 950k).

## 13 Oct: G6 / Ex6 / DISENT-ph combined -- `REPORT-2026-10-13.md`
- All G6 (W, photon), Ex6 and DISENT-ph seeds done and finite, no reruns. Tools: `combine6.py SET` (G6W|G6P|Ex6), `gen_tables6.py SET`; results `results/{G6W,G6P,Ex6}_6_*`, `outliers6_*`, `rx_tail_compare6.txt`.
- W NLO references of the instructions are in combine_tcut units already (used without dividing by x); photon DISENT = raw coefficient / x, plain mean over 40 seeds.
- Ex6 r-x with VTARGET=0: the printed `iteration` values are not on the scale of the result (median max|iteration| 170 vs 14 in bin5), so the 20x-median flag rate (69/2000 vs 86/875) is not like for like; per-seed sd of the total row is 0.73x bin5 at tau_cut 1e-5 and 3-7x at larger tau_cut.

## 13 Oct (evening): bin7 and Ex7 pilot (UPDATE 13 Oct); Ex7 production NOT fed (lowest priority, coordinator's change of plan)
- bin7 = `$P/bin7/{sliced21,nlo31}` (MD5SUMS there), built from HEAD e6d834c (contains 03f6981, e82e1f8; checked with merge-base --is-ancestor) with
  `nnlo21/build_sliced21.sh $P/src/MCFM-10.3 $P/src/nnlojet-v1.0.2/install/lib64 $P/builds/sliced21-e6d834c` as a 16-core job (48921377, ~1 min). No tau2_nlo copied.
- Checks `$P/runs/chk7/chk.sh`: sliced21 b0 (30k x 3) and nlo31 lo mode 2 (100k x 3), VEGAS_EQUAL=2, seed 7, without VTARGET: bin6 and bin7 outputs bit-identical (cmp).
  r-x short run (100k x 3) with `ZFIX=3 P2BSLICE=1 P2BEXTRA=1 VTARGET=-1`: finite.
- `mk_ex7.py pilot|prod` -> `$P/runs/ours7/Ex7[-pilot]/<part>`; Ex7: lo-x 200, vi-x 250, kp-x 200 (seeds 5001-), r-x 2000 (5001-7000), env as Ex6 but VTARGET=-1, bin7.
  Seeds disjoint from Ex6 (lo/vi/kp 2001-2250, r 1001-3000), bin5 E extras (r-x 1001-1875, lo/kp/vi-x 2001-2062) and the pilots (Ex6 9001-9020, Ex7 9101-9120).
  Groups file `$P/runs/ours7/groups7.txt` (copy of Ex7/groups.txt); job lists created, NO feeder started.
- Pilot: 20 r-x seeds 9101-9120, array 48921477 (dis-x7p-r-x), all done and finite, 0.65 core-h per job (Ex6 pilot 1.43). `scatter7.py DIR` prints sd and 1.48 MAD of the total row.
  Total row, tau_cut 1e-3 / 1e-4 (sd, 1.48 MAD, pb): Ex7 pilot 481/267 and 442/62.5; Ex6 (2000 seeds, VTARGET=0) 3942/30.4 and 3550/12.1; Ex6 pilot (20 seeds) 73.5/38.4 and 27.1/11.3.
  Ex7 pilot median -119 / -38 (Ex6 -6.3 / -1.5). The robust scatter is 9x (1e-3) and 5x (1e-4) LARGER than Ex6's, so the criterion (below Ex6) is not met; only the 2000-seed Ex6 sd (tail-driven) is above.
  Printed final `iteration` values are ~2.3e4 (VTARGET=-1 scale, not comparable).
- File projection: Ex6 seed dir = 6 files (attempts, cmd.sh, done, pattern, run.log, time.log) = 7 entries; Ex6 total 18,559 for 2,650 seeds; Ex7 would add ~18.6k.
  /ptmp quota at 925,130 (851k at the start of the day; rising with no job of mine running: not from this session beyond ~5k), so ~944k with Ex7 (< 1.0M).

## 13 Oct: J7 (ZEUS jets, window (x,Q2)=(0.01,400), bin6, ZFIX=2) and N7 (NNLOJET W top-up)
- Generators `mk_j7.py PH|W` (-> `$P/runs/ours7j/{PH,W}/<part>/{jobs.list,pilot.list}`, groups.txt, groups-pilot.txt) and `mk_n7.py`
  (new seeds 11001- in `$P/runs/nnj-W/prod/<part>/`, list `prod7.list`, groups `$P/runs/nnj-groups7-{RRa1,rest}.txt`).
- Window check (`ours7j/chk`, nlo31 lo, `2 0.01 400`, ZFIX=2): log prints `ZFIX window Q2, y: 3.998000E+02 4.002000E+02 3.950593E-01 3.954545E-01`.
- J7 commands as in the instructions: PH EW31="0 0", tables zeus5, seeds 1001- (b2 200, vi 240, kp 200, r 1600, rcorr 200, d2 100); W EW31="2 0",
  zeus5w, seeds 5001- (150, 180, 150, 1200, 200, 100). The pilot = the first 5 seeds of each part (they count as production).
- Pilot CPU [core-h/job]: PH b2 0.933, vi 1.858, kp 0.111, r 1.382, rcorr 0.454, d2 0.489 -> PH 3.0k core-h (instructions ~6k desktop-h).
  W b2 0.872, kp 0.105, r 1.852, rcorr 0.687, d2 0.662, vi pilot not finished when W was started (G6W vi 3.1 -> ~0.6k) -> W ~3.1k (instr. 5k).
  Gates passed (-50%, -38%).
- N7: seeds disjoint (existing 1001-3530, 9001-9005). Targets RRa_1 5500, RRb_1 5100, RRa_2 1850, RRa_3 420, RRa_5 550, RRb_5 940 -> new jobs 5228, 2570,
  929, 314, 412, 706 = 10159. Existing mean CPU/job (all seeds, size_nnlojet.py) 3295, 2974, 2903, 3859, 1111, 461 s -> 4785+2123+749+337+127+90 =
  8.2k core-h (instr. 8k, +3%). The errors follow sigma/job: RRa_1 887 fb -> 12 fb at 5500.
- nnj_clean.sh (called by run_task.sh at the end of a job) covers N7 item 2: 8 entries per finished job dir (dir, done, time.log, logs.tgz, 4 .dat), 9 while running.
  Measured: finished NNLOJET job 8 entries, ours job dir 8 (mkjobs precreates dir, cmd.sh, pattern = 3).
- File projection: usage 851k (start), 897k before RRa_1, 914k at 2.4k running RRa_1 jobs; after the first 1462 finished RRa_1 jobs 923.6k (other
  agents' files also grow); 947.5k at 14:00 with RRa_1 4.8k done. Projection of the rest: N7 rest 4931 x 5 + J7 W 1980 x 5 + PH tail ~ 0.99M
  (< 1.0M, soft 1,048,576). Z (if submitted) would add ~37k.
- Submitted 13 Oct: pilots PH 48921395-400 (names dis-jpp-*), N7 RRa_1 feeder dis-feeder-n7a (48921516, array 48921519, nice 500 to let PH go first),
  PH feeder dis-feeder-j7 (48926765; arrays dis-jp-*, vi list jobs-rest.list = seeds 1006-1240), W pilots dis-jwp-* 48926776-81, W feeder
  dis-feeder-j7w (48930792; arrays dis-jw-*, vi jobs-rest.list 5006-5180), N7 rest feeder dis-feeder-n7b (48932761; arrays dis-n7-RRb_1 48932806,
  RRa_2 07, RRa_3 08, RRa_5 09, RRb_5 10). Deadline 17 Oct. Note: alma outside the et nodes had only ~2.5k slots for us; et[07-40] idle but excluded.
- N7 Z estimate (not submitted; existing done seeds, mean CPU): RRa_3 721 -> x6 = +3605 jobs x 2.30 = 8.3k, RRa_2 131 -> x2 = +131 x 1.56 = 0.2k,
  RRb_3 882 -> x2 = +882 x 1.06 = 0.9k: 9.4k core-h (instructions 13k), ~4.6k jobs, ~37k files.
- 13 Oct 15:00, N7 Z approved by AK (item 5): `mk_n7.py Z` (same script, argument Z; seeds 11001- disjoint from existing Z 1001-, 9001-9005) in `$P/runs/nnj-Z/prod`,
  lists `prod/<part>/prod7.list`, groups `$P/runs/nnj-groups7-Z.txt`. New jobs RRa_3 3605 (to 4326, 12 h limit), RRa_2 131 (to 262), RRb_3 882 (to 1764) = 4618,
  est. 9.4k core-h. Feeder dis-feeder-n7z 48933872 (deadline 17 Oct), arrays dis-n7z-RRa_3 48933875, RRa_2 48933876, RRb_3 48933877, nice 500 (below J7 W and N7 W rest).
  Files: usage 667k before, 681k after (dirs precreated); final +23k for Z, ~0.77M with the rest of J7/N7 and Ex7 (~18.6k).

## 13 Oct (night): Ex7 pilot extended to 200 r-x seeds, decision of UPDATE 13 Oct (d): NOT fed
- `mk_ex7.py pilot200` extends `Ex7-pilot/r-x/jobs.list` to seeds 9101-9300; 180 new seeds as array 48932746 (nice 0, bin7, same command/env as the pilot); 9 tasks lost to NODE_FAIL on gt34 resubmitted (48934847). All 200 done, finite. ~0.65 core-h per job.
- `pilot200.py` -> `results/ex7_pilot200.txt`. Total row (CELL 0.05 0.50), pb, plain means:
  | set | tau_cut | N | mean +- sd/sqrtN | sd | 1.48 MAD | median | frac > 20 MAD |
  |---|---|---|---|---|---|---|---|
  | Ex7 | 1e-3 | 200 | -110.6 +- 25.0 | 353 | 145.5 | -138.6 | 0 |
  | Ex6 (2000) | 1e-3 | 2000 | 20.8 +- 88.1 | 3942 | 30.4 | -6.3 | 0.039 |
  | Ex6 first 200 | 1e-3 | 200 | 666 +- 741 | 10470 | 32.1 | -5.2 | 0.045 |
  | Ex7 | 1e-4 | 200 | -35.7 +- 14.3 | 202.7 | 66.5 | -26.0 | 0.005 |
  | Ex6 (2000) | 1e-4 | 2000 | -77.4 +- 79.4 | 3550 | 12.1 | -1.45 | 0.038 |
  | Ex6 first 200 | 1e-4 | 200 | -832 +- 786 | 11110 | 12.5 | -2.6 | 0.055 |
- (a) sd7^2*0.65 = 8.1e4 (1e-3), 2.7e4 (1e-4) against sd6^2*1.43 = 2.2e7, 1.8e7: holds by a factor ~270 (needs only a few % of Ex6's variance). (b) 0 and 0.5% against 3.9% and 3.8%: holds. (c) Ex7-Ex6 = -131 +- 92 (pull -1.43) at 1e-3, +42 +- 81 (+0.52) at 1e-4: consistent but 1e-3 is outside 1 sigma, and Ex6's error is a tail-driven sd/sqrtN (weak test).
- Caution: Ex7 medians (-139, -26) differ from Ex6's (-6, -1.5) by far more than the Ex6 MAD (the Ex6 distribution is mostly near zero with rare large seeds carrying the mean); Ex7 1.48 MAD is 4.8x (1e-3) and 5.5x (1e-4) larger than Ex6's. Ex7 mean is -4.4 sigma/-2.5 sigma from zero, Ex6's is not significant.
- Verdict: (a), (b) clear; (c) not clear at 1e-3 -> feeder NOT started (rule: marginal = do not feed). Production lists and groups7.txt remain prepared.

## 13 Oct (evening): UPDATE 13 Oct (e): Ex7 production fed
- `submit.py` takes the default nice from the environment variable `NICE` (default 100 as before). Feeder chain `dis-feeder-ex7` (first job 48939329, then hourly; self-resubmitting job 48939332 pending), `GROUPFILE=$P/runs/ours7/groups7.txt`, `NICE=500`, DEADLINE 17 Oct (00:00), excluding et[07-40].
- Arrays (all nice 500): lo-x 48939333 (200), vi-x 48939334 (250), kp-x 48939336 (200), r-x 48939344 (2000); seeds 5001-, bin7, `ZFIX=3 P2BSLICE=1 P2BEXTRA=1 VTARGET=-1`.
- Files before: 715k (limit 1.0M); Ex7 adds ~18.6k.
- CPU: r-x 2000 x 0.65, vi-x 250 and lo/kp-x from the Ex6/E pilots: ~1.6-1.8k core-h. Finish depends on how soon the nice-500 tasks get slots behind J7/N7.

## 14 Oct: completeness, J7/N7/Ex7 combination (REPORT-2026-10-14)
- All lists of J7 PH/W, N7 W/Z and Ex7 100% done with finite outputs (`check7.py`; the "active" counts of the last feeder passes were stale); no resubmission needed. RRa_1 s14135, s14146 complete.
- J7 NLO (b1, lo; not in the J7 lists): `mk_j7nlo.py PH|W` -> `ours7j/<SET>/{b1,lo}` 100 seeds each, groups `groups-nlo.txt`, direct `submit.py` (arrays 48950862-65), 0.53 / 0.15 core-h per job.
- `combine7n.py`, `combine7j.py`, `combine7e.py`, `gen_tables7.py`; results in `results/{N7,J7,Ex7_7}*`, raw J7 ZCELL in `results/j7/*.zcell.gz` (3.5 MB).
- N7 Z RRa_2 seed 11083 has 183.5 pb (cross value); heavy-tailed NNLOJET parts are normal, see `results/N7_{W,Z}_outliers.txt`.
- `/ptmp` file count 756k at the end of the round (limit 1,048,576).

## 14 Oct: J8 (UPDATE 14 Oct (b) + (d)), bin8 and bin9
- bin8 = `$P/bin8/{sliced21,nlo31}` from HEAD 971f224 (contains 37301af, ZWIN). Window checks passed (`$P/runs/chk8`): sliced21 b1 `ZFIX=2 ZWIN="0.1 5000"` prints `x 0.100000 Q2 5000.00`; nlo31 lo `2 0.1 5000` prints `ZFIX window Q2, y:  4.997500E+03  5.002500E+03  4.938241E-01  4.943182E-01`.
- bin8 pilots (5 seeds/part, `$P/runs/ours8j/PH`, ZEUS rows only, cost estimate only): core-h/job kp 0.113, rcorr 0.478, d2 0.479, b1 0.562, lo 0.141 (J7 PH 0.111, 0.454, 0.489, 0.53, 0.15: within +-6%); b2, vi, r still running at the start of production.
- bin9 = `$P/bin9/{sliced21,nlo31}` (MD5SUMS) from HEAD a3d5eee (contains c046b57, ZCUTS; build `$P/builds/sliced21-a3d5eee`, job 48952867). Check (`$P/runs/chk9/chk.sh`, sliced21 b0 30k x 3, seed 7, ZFIX=2): with ZCUTS=1 21 ZCELL rows, et4 rows equal to the ZEUS total and ptavg 8-15 rows, the 15 ZEUS rows and RESULT identical to ZCUTS unset; without ZCUTS identical to bin8. `$P/bin9/READY` marks the check passed.
- J8 production on bin9 with ZCUTS=1: `mk_j8.py PH`, `mk_j8nlo.py PH` -> `$P/runs/ours8j/PH9` (b2 200, vi 240, kp 200, r 1600, rcorr 200, d2 100, b1 100, lo 100; seeds 1001-), groups `groups-prod.txt`. Feeder `dis-feeder-j9` (first pass by hand 14 Oct 09:38, nice 100, DEADLINE 17 Oct 12:00): arrays b2 48952877, vi 48952978, kp 48952979, r 48952980, rcorr 48952981, d2 48952982, b1 48953082, lo 48953083.
