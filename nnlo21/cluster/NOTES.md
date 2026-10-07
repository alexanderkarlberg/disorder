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

### CPU (user time, approximate)
A ≈ 16k core-h (incl. warmups ≈ 2k); B ≈ 7.5k; F ≈ 2.2k; Bx ≈ 1.4k;
D ≈ 1.9k (incl. the duplicate lodef 0.54k and invalid locorr 0.1k); pilots and
diagnostics ≈ 0.5k; lost to the 10:18 outage and node failures ≈ 3k.
