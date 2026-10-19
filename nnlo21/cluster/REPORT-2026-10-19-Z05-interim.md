# nnlo21 cluster: Z05 interim report (19 Oct 2026, r at 35 %)

Branch `2026-10-nnlo21-cluster`. UPDATE 18 Oct (Z05: N3LO 1+1 P2B in the zeus05 set-up of 1803.09973, 7 scales). Facts only.
Summary file: `results/Z05-summary-interim.json.gz` (`summary_z05.py` from the ncc branch, run 19 Oct 10:52 on every finished job).

## 1. Build and checks (item 1)

bin12 = `$P/bin12/{sliced21,nlo31}` from 9670206 (HOPPET 2.3.0, NNPDF31_nnlo_as_0118). MD5 sliced21 a8b80120a15fd6c4fc1f607eec3c1b7f, nlo31 5bdc422cf65caad022045adcf837dfa5.

- **`nlo31 chk 1 1`:** kp_conv q -> q (b, g, lp) / closed form - 1 at xi 0.01 / 0.1 / 0.5: max 1.55e-8 / 2.1e-10 / 1.2e-11 (as bin11).
- **zeus05 r test** (`LABSET=zeus05 SCALES=3 VTARGET=0 VEGAS_MIX=0.2 VEGAS_EQUAL=2 nlo31 r 20000 2 61 M3`): 31 LCELL and 217 SCELL lines. The SCELL lines with scales (1.00, 1.00) equal the LCELL lines exactly.
- **J8 r seed, bin12 against bin11: not identical.**
  - Command: the PH9 command without ZCUTS (removed in 42a6c7a), `ZFIX=2 EW31="0 0" P2BSLICE=1 VEGAS_EQUAL=2 nlo31 r 100000 3 1001 1d-9 2 0.1 5000 psmc 0 0 1d-12`.
  - Iteration 1 is identical. From iteration 2 the numbers differ:

    | | bin11 | bin12 |
    |---|---|---|
    | iteration 2 | −1.06302372e-4 | −1.06330832e-4 |
    | result | −9.54785e-5 ± 2.41e-5 | −9.53925e-5 ± 2.39e-5 |
    | events dropped (garbage dipole) | 35 | 37 |

  - bin12 is reproducible: two runs in separate directories are cmp-identical.
  - A build of 42a6c7a (the clean-up) gives bin11's numbers exactly; only the removed `P2BDROP 1` header line differs. So the change lies in 42a6c7a..9670206. d5c2021 (born31 |z|^2) records last-bit changes (<= 2e-10) that VEGAS adaptation amplifies.
  - User time for this run: bin12 170 s, bin11 239 s.

## 2. Tables (item 2)

- **Beam grids:** `$P/tables/{z05,z05_x2,z05_xh,z05_p2,z05_ph}_{0..43}.tab`, 44 files each. All 220 jobs completed.
- **Soft table:** `$P/tables/softG.tab` is a copy of the earlier `soft.tab`. bin12 reads it without rebuilding.

## 3. Runs (items 3, 4)

### Job layout

The /ptmp file count stood at 982k of the 1.05M soft limit.

- **One file per job:** each job writes `$P/Z05/<part>/<part>_<seed>.out`, holding the run output and the `/usr/bin/time -v` report.
- **Success test:** the file is renamed from `.out.tmp` only on exit 0 with 217 SCELL lines.
- **Scripts:** `z05.py` (submit, status), `run_z05.sh`, `feeder_z05.sh` (hourly).
- **Excluded nodes:** et[07-40] (see 3.4), ct30, gt34, gt37, kt02, kt06, and gt25, gt27-29, gt39-40, kt01, kt03, kt10 as before.

### Wave 1 and the item-4 decision

Wave 1 was 5 % of every part, 1383 jobs. It was evaluated on 18 Oct 22:55:

- **Projected total:** 30.7k core-h, so no r cut.
- **Means at 1e-4:** within 1.4 sigma of the MPP pilot.
- **Decision:** full ranges submitted.

### State on 19 Oct 10:52

**Done:**

| part | done of total |
|---|---|
| vi | 2729 / 2729 |
| kp | 2479 / 2479 |
| rc1 | 1000 / 1000 |
| b0 | 200 / 200 |
| b1 | 400 / 400 |
| lo | 400 / 400 |
| b2 | 728 / 729 |
| rc2 | 994 / 1000 |
| r | 6545 / 18700 (35 %) |

About 2400 r jobs are running.

**Means at tau_cut 1e-4 against the MPP pilot** (central scale; sigma per bin in pb; rows 1 = >= 1 jet, 8 = eta 2-2.5, 10 = E_T 6-8):

| part | n | row 1 (pull) | row 8 (pull) | row 10 (pull) |
|---|---|---|---|---|
| b2 | 728 | 73268.0 +- 6.3 (-0.4) | 11130.3 +- 0.9 (-0.4) | 37001.5 +- 42.8 (-0.0) |
| vi | 2729 | -13472.4 +- 12.0 (-1.3) | -2836.5 +- 1.9 (+0.9) | -5557.5 +- 27.9 (-0.5) |
| kp | 2479 | 6881.0 +- 3.4 (-1.0) | 1235.7 +- 0.5 (-0.9) | 3591.4 +- 3.9 (-1.1) |
| r | 6580 | -66519 +- 267 (-1.3) | -9435.6 +- 10.6 (+0.8) | -34822 +- 303 (-1.2) |

**CPU per job** (user+sys, core-h, by node type; n in brackets):

| part | mean | ct | gt | kt | et | core-h so far | projected |
|---|---|---|---|---|---|---|---|
| r | 1.384 (max 4.88) | 1.325 (4464) | 1.551 (1554) | 1.408 (545) | 1.138 (26) | 9121 | 25887 |
| vi | 1.219 | 1.179 | 1.425 | 1.226 | 1.010 | 3326 | 3326 |
| b2 | 2.061 | 2.068 | 2.225 | 2.175 | 1.591 | 1500 | 1502 |
| b1 | 1.681 | 1.642 | 1.936 | 1.809 | 1.316 | 672 | 672 |
| lo | 0.928 | 0.883 | 1.063 | 0.947 | 0.739 | 371 | 371 |
| kp | 0.120 | 0.120 | 0.145 | 0.126 | 0.102 | 298 | 298 |
| rc1 | 0.278 | 0.251 | 0.295 | 0.267 | 0.212 | 278 | 278 |
| rc2 | 0.244 | 0.231 | 0.281 | 0.253 | 0.203 | 243 | 244 |
| b0 | 0.102 | 0.100 | 0.113 | – | 0.092 | 20 | 20 |

- **Totals:** projected 32.6k core-h; used so far 15.8k.
- **Against the instructions' estimates:** rc1 and rc2 are about 4x below their 1.0 h, vi 0.7x, b2 0.8x, r 0.7x.
- **et jobs:** the et jobs (all wave-1 non-r parts, 26 r test jobs) had the lowest CPU per job.

### Failed and hung jobs

No run has failed: no non-zero exit and no missing SCELL lines.

Hung tasks never started their job script: they ran for 1-2 h with no attempts line and no output. All were cancelled by id and resubmitted:

| time | tasks | nodes / part |
|---|---|---|
| 18 Oct 22:50 | 164 | et13, et15, et26 (wave 1: kp, rc1, rc2, b2, vi) |
| 19 Oct 05:05 | 1 | gt37 (b2) |
| 19 Oct 07:06 | 4 of a 30-task et test | et17 (r) |
| 19 Oct 10:51 | 59 | gt34, gt37, kt02, kt06 (r; started about 09:13) |

68 r attempts were killed by our own cancellation of a duplicate r array (18 Oct 22:55); they were resubmitted.

**et nodes:** 30 r seeds were tested on et17 on 19 Oct. 26 finished and 4 hung, and the idle et nodes then showed "not responding". et stays excluded.

### Seeds beyond 20 x 1.48 MAD

Total row (>= 1 jet), tau_cut 1e-4:

| r | rc1 | rc2 | b2, vi, kp, b0, b1, lo |
|---|---|---|---|
| 0.0012 | 0.0160 | 0.0181 | 0 |

## 4. N3LO P2B part per row (`summary_z05.py`, central scale, sigma per bin in pb; interim)

b2 + vi + kp + r, plus rc1 + rc2 below 1e-4; r with 6545 of 18700 seeds. Row names as printed by the script.

| row | 1e-3 | 5e-4 | 2e-4 | 1e-4 | 3e-5 | 1e-5 |
|---|---|---|---|---|---|---|
| >=1 jet | 682 +- 98 | 541 +- 145 | 215 +- 239 | 160 +- 269 | 80 +- 393 | 58 +- 654 |
| eta 1 | 49 +- 69 | 127 +- 90 | 64 +- 141 | -17 +- 186 | -1589 +- 1706 | 261 +- 828 |
| eta 2 | 56 +- 61 | -77 +- 79 | -293 +- 174 | -229 +- 170 | 1212 +- 1697 | -584 +- 786 |
| eta 3 | 157 +- 74 | 186 +- 103 | 253 +- 140 | 237 +- 167 | 100 +- 241 | 43 +- 350 |
| eta 4 | 123 +- 77 | -2 +- 95 | -140 +- 135 | -34 +- 159 | 151 +- 220 | -36 +- 362 |
| eta 5 | 182 +- 54 | 172 +- 60 | 156 +- 74 | 67 +- 87 | 65 +- 120 | 203 +- 150 |
| eta 6 | 105.8 +- 11.7 | 128.3 +- 15.8 | 145.1 +- 22.5 | 93.2 +- 25.5 | 154.7 +- 39.4 | 214.5 +- 63.0 |
| eta 7 | 89.3 +- 5.3 | 90.3 +- 6.7 | 89.9 +- 8.8 | 94.3 +- 10.9 | 100.8 +- 15.9 | 80.8 +- 26.4 |
| eta 8 | 57.6 +- 3.0 | 56.9 +- 3.6 | 65.8 +- 4.6 | 61.3 +- 5.6 | 67.9 +- 9.1 | 50.8 +- 15.2 |
| E_T 1 | 378 +- 120 | 201 +- 184 | -25 +- 271 | 215 +- 309 | 398 +- 486 | 455 +- 702 |
| E_T 2 | 268 +- 74 | 340 +- 123 | 234 +- 164 | 14 +- 196 | 134 +- 460 | -342 +- 608 |
| E_T 3 | 141.7 +- 27.1 | 128.9 +- 46.8 | 68.9 +- 95.8 | 22.8 +- 114.6 | -376.3 +- 342.8 | 167.9 +- 356.3 |
| E_T 4 | 39.0 +- 10.3 | 18.1 +- 32.1 | 67.4 +- 25.5 | 40.7 +- 33.2 | 89.9 +- 45.1 | -14.6 +- 160.8 |
| E_T 5 | -4.0 +- 3.9 | -6.5 +- 4.6 | -3.4 +- 5.6 | -13.4 +- 9.3 | 18.5 +- 13.5 | -34.0 +- 23.1 |
| E_T 6 | -2.2 +- 1.1 | -0.9 +- 1.3 | -0.9 +- 1.7 | -5.8 +- 2.2 | 0.8 +- 4.2 | -2.5 +- 6.4 |
| E_T 7 | 0.1 +- 0.1 | -0.0 +- 0.2 | -0.5 +- 0.3 | -0.4 +- 0.4 | -2.3 +- 2.0 | 1.2 +- 1.5 |
| E_T 8 | 0.0 +- 0.0 | 0.0 +- 0.0 | -0.0 +- 0.0 | -0.0 +- 0.0 | -0.1 +- 0.1 | 0.0 +- 0.1 |

At 3e-5, rows eta 1 and eta 2 (-1589 +- 1706, +1212 +- 1697) have errors 2-10x those at 1e-4 and 1e-5 in the same rows.

## 5. Next

The feeder continues: about 12k r seeds remain, at about 2400 running at a time. The final report follows at the end of the run.
