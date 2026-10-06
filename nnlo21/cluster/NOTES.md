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
