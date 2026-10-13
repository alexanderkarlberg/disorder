# NNLO DIS 2+1 by τ₂ slicing: cluster production — instructions for Claude

**Written 6 Oct 2026 on thA371a (draft of the night, updated in the evening
with projected slicing and the new technical cut).** For the Claude session that AK starts on the cluster (Slurm, partition
`alma`, 24 h wall time, at most 10k running and 30k queued jobs). Read this
file, `CLAUDE.md`, `docs/nnlo21-plan.md` and the nnlo21 entries of
`docs/notebook.md` (3–6 Oct) first. The page with the current status:
https://claude.ai/artifact/ChKctwxDdEau7H9DGMxwyW

## UPDATE 13 Oct (d): Ex7 pilot to 200 r-x seeds, new decision rule (after the bin7 pilot, 16766e3)

The 20-seed pilot (9101–9120) failed the criterion of UPDATE 13 Oct (robust scatter below Ex6's), and the
production was rightly not fed. That criterion was the wrong one: the error of the mean is set by the full
sd, tails included, and Ex6's sd comes from rare tail seeds (total row at 1e-4: sd 3550 over 2000 seeds,
27 over its 20-seed pilot), which 20 seeds cannot see. If Ex7 had no such tails, its sd of ≈ 450 would give
≈ 10 pb on 2000 seeds against Ex6's 79, at 0.65 instead of 1.43 core-h per job.
1. Extend the pilot to 200 r-x seeds: 9121–9300 (with 9101–9120), same command and bin7, ≈ 120 core-h.
   Priority: after J7 PH, N7 RRa_1 and J7 W (it is small; it may fill gaps).
2. Decision, total row at τ_cut 1e-3 and 1e-4: **feed Ex7** (the r-x seeds 5001–7000 and the lo/vi/kp-x
   parts) if (a) sd_Ex7² × 0.65 < sd_Ex6² × 1.43 with sd_Ex6 from the 2000 Ex6 seeds, (b) the fraction of
   seeds with |value − median| > 20 × (1.48 MAD) is not larger than Ex6's, and (c) the Ex7 and Ex6 means
   agree within their errors (sd/√N; same integrand, so they must). Otherwise **drop Ex7**: Ex6 + bin5 are
   the final E extras, and only the report of point 3 below (for Ex6 + bin5) is needed.
3. Report the 200-seed numbers (mean ± sd/√N, sd, 1.48 MAD, median, outlier fraction; Ex6 and its first
   200 seeds alongside) before production; feed it without waiting if (a)–(c) hold clearly.

## UPDATE 13 Oct (c): J7 — ZEUS jets in a narrow (x, Q²) window, photon and W, FIRST (AK: as soon as possible)

Why: at the fixed point (G6) the τ_zQ total is flat, but W with the ZEUS jet selection in a narrow window
around (x, Q²) = (0.01, 400) drifts like the full ZEUS set (MPP `runs/zw2`: +15–21% of the NNLO per unit
ln τ_cut, 2–3σ, soft jet bins). J7 gives the same for the photon at 4× the MPP statistics, plus more W, so
the two bosons can be compared jet bin by jet bin at fixed kinematics.
1. bin6 as it is (7ddc208; the window mode `ZFIX=2` of sliced21 and nlo31 is in it; `VTARGET` is not
   used). nlo31 takes the window centre from its arguments: `2 0.01 400`, **not** `2 0 0` (empty window).
   Check in one short run that the log prints `ZFIX window Q2, y:  3.998000E+02  4.002000E+02 ...`.
2. Commands (`$P` paths; `PH` = photon: `EW31="0 0"`, tables `$P/tables/zeus5`; `W`: `EW31="2 0"`,
   tables `$P/tables/zeus5w`; every job with `ZFIX=2`):
   ```
   b2:  sliced21 b2 5000000 6 SEED zeus <tables> $P/tables/soft.tab
   vi:  P2BSLICE=1 nlo31 vi 1000000 6 SEED 1d-9 2 0.01 400 psmc 0 0 1d-12
   kp:  P2BSLICE=1 nlo31 kp 2000000 6 SEED 1d-9 2 0.01 400 psmc 0 0 1d-12
   r:   P2BSLICE=1 VEGAS_EQUAL=2 nlo31 r 1000000 6 SEED 1d-9 2 0.01 400 psmc 0 0 1d-12
   rcorr, d2: as r, with the TECHDIFF settings of G6 (rcorr: TECHMIN "1d-10 1d-8"; d2: TECHMIN
        "1d-11 1d-9", TECHREF "1d-10 1d-8")
   ```
   PH: b2 200, vi 240, kp 200, r 1600, rcorr 200, d2 100 (seeds 1001–). W: the same counts minus MPP's
   (b2 150, vi 180, kp 150, r 1200, rcorr 200, d2 100), on seeds 5001– (MPP used 9001–9700). MPP's local
   jobs took 0.2 h (kp), about 1.5–2 h (vi), about 2.5–4 h (r) and about 5.5 h (b2) per job on the
   desktops: about 6k core-h for PH and 5k for W. Pilot of 5 jobs per part for the timing; submit if within
   +30%.
3. Order: J7 PH, then N7 RRa_1, then J7 W, then Ex7 and the rest of N7.
4. Report (as soon as PH is done; W and the others can follow): per ZEUS observable bin (cross, q2,
   ptavg_12, m12) the NNLO at every τ_cut (with d2 below 1e-4), the NLO, and the slopes per unit ln τ_cut
   in % of the bin NLO, windows 5e-3→5e-4, 2e-3→2e-4, 1e-3→1e-4, 1e-4→1e-5, pulls as in REPORT-2026-10-13.
   PH and W side by side, and the ratio of their slopes per bin. Put the raw `ZCELL` outputs in
   `results/j7/` so that MPP can combine them with `runs/zw2`.

## UPDATE 13 Oct (b): N7 — NNLOJET W (ZEUS dijets) to 10% at NNLO, the arbiter for the jet drift

Why: our ZEUS-jet NNLO rises in ln τ_cut at a constant rate over 2.5 decades (W total −1.2 to −1.3% of
LO per unit ln, 7–11σ; REPORT-2026-10-12 §5). A constant slope is not what power corrections do, which
points to a wrong ln τ_cut coefficient below the cut. NNLOJET decides which end is right, but its W
NNLO is 0.176 ± 0.056 pb (32%). Ours is 0.159 ± 0.009 at 1e-3, 0.180 ± 0.011 at 5e-4 and 0.357 ± 0.045 at
1e-5, so the two ends differ by 0.18 pb, about 3σ of NNLOJET's error. Brought to ≈ 0.017 pb, NNLOJET
separates the two ends at > 5σ.
1. Same runcard, build and workflow as the W production of UPDATE 10 Oct (`nnj5/W`); new seeds only, and
   combine them with the existing ones. Top up the parts that carry the error (`nnj5/W/parts.txt`):
   RRa_1 272 → 5500 seeds (error 54 → ≈ 12 fb; ≈ 4.7k core-h), RRb_1 2530 → 5100 (9.8 → 6.9 fb;
   2.1k), RRa_2 921 → 1850 (0.7k), RRa_3 106 → 420, RRa_5 138 → 550, RRb_5 234 → 940 (together
   ≈ 0.6k). Total ≈ 8k core-h, target NNLO error ≤ 0.02 pb. Size it with `size_nnlojet.py` from the
   existing CPU/job; submit if within +30%, otherwise report first.
2. **/ptmp file count** (11 Oct incident): keep only what `combine_nnlojet.py` needs per job and pack the
   rest per part (tar) as each job ends; check the count after the first 500 jobs.
3. Order: see J7 (UPDATE 13 Oct (c)) above.
4. Report: NNLOJET W per observable bin (cross, q2, ptavg_12, m12) with the new errors, against ours
   (bin5 `W`) at every τ_cut, pulls, and per bin the τ_cut at which ours crosses NNLOJET. In particular,
   for the total and ptavg 30–60, say which τ_cut range is compatible.
5. If the cluster is still free after N7: the same for Z (`nnj5/Z`; error 1.17 pb from RRa_3 0.93,
   RRa_2 0.40, RRb_3 0.38). RRa_3 ×6, RRa_2 and RRb_3 ×2 give ≈ 0.6 pb for ≈ 13k core-h. Report the
   estimate first, do not submit.

## UPDATE 13 Oct (after REPORT-2026-10-13): Ex7 = Ex6 again with `VTARGET=-1`

REPORT-2026-10-13 received (MPP). G6 is complete, no further G runs. Ex6 ran with `VTARGET=0` (it had
started before the correction of 13 Oct, 02:30), and the report confirms the MPP finding: the extras'
scatter is 2–7× larger than bin5's at τ_cut ≥ 3e-5 and smaller only at 1e-5, the column VEGAS adapts to.
The Ex6 results stay valid (the target changes the sampling, not the integrand) and are kept.
1. Build **bin7** from `origin/2026-10-nnlo21-ncc` (e82e1f8 or later; it contains 03f6981, which bin6,
   built from 7ddc208, does not). Check: short b0/lo runs bit-identical to bin6 without `VTARGET`, and one
   short r-x run with `VTARGET=-1` finite.
2. **Ex7**: the Ex6 commands and seed counts (lo-x 200, vi-x 250, kp-x 200, r-x 2000; ≈ 2.6k core-h) with
   `ZFIX=3 P2BSLICE=1 P2BEXTRA=1 VTARGET=-1` and bin7, on new seeds disjoint from Ex6 (lo/vi/kp-x 2001–,
   r-x 1001–3000) and from the bin5 E extras: e.g. 5001– for every part. Pilot of 20 r-x seeds first;
   continue if the per-seed scatter of the total row at 1e-4 and 1e-3 is below Ex6's.
3. Report (new file, e.g. REPORT-2026-10-14) for Ex7 alone and for the per-cell inverse-variance combination
   Ex6 + Ex7: NNLO_P2B and NLO_P2B against NNLOJET, slopes, item-6 window. Per-seed scatter (sd and
   1.48 MAD) per part and τ_cut for bin5, Ex6 and Ex7 side by side.
   **Flag outliers on the final per-seed values, not on the `iteration` lines**: with `VTARGET` set,
   VEGAS integrates the target (sum of |cells|), so the printed iterations are on that scale and are not
   comparable with runs without `VTARGET` (the caveat of REPORT-2026-10-13 §3.1).

## UPDATE 12 Oct (after REPORT-2026-10-12): two diagnostic sets, Ex6 and G6 (AK: the cluster is empty, go)

Context (notebook 11–12 Oct, REPORT-2026-10-12): the ZEUS-type jet NNLO (B, Z, W) rises in ln τ_cut
(W: 7–11σ, not levelling off at 1e-5); the P2B − plain term is not the cause (MPP, 11 Oct); E plain falls
towards NNLOJET, E P2B agrees but its extras were too noisy. In mode 3 we found (12 Oct) that heavy tails
of r come from the VEGAS target seeing only one row: new option `VTARGET=0` (VEGAS adapts to the sum of
|cells| over all rows at the smallest τ_cut; sampling only, default bit-identical): r errors 3–14× smaller
in jet bins, unbiased. The E extras (target = one τ_zQ row) have the same problem.

1. Build **bin6** from `origin/2026-10-nnlo21-ncc` (92195a5 or later; `nnlo21/build_sliced21.sh`).
   Check: b0 and lo short runs (mode 2, as for bin5) bit-identical to bin5 without `VTARGET`.
2. **Ex6 — E P2B extras, more statistics** (the plain parts of E from bin5 stay; P2B = bin5 plain +
   Ex6 extras): `ZFIX=3 P2BSLICE=1 P2BEXTRA=1`, **default VEGAS target — do NOT set `VTARGET=0`**
   (correction 13 Oct, 02:30: for the extras, θ(T ≤ τ_cut), `VTARGET=0` adapts to the smallest region
   and is worse than the default in many cells, MPP test `runs/exv`; if Ex6 already runs with
   `VTARGET=0`, stop it). **Update 13 Oct, 08:30 (AK OK): use `VTARGET=-1`** (all cells; per-seed spread
   of r-x 2–20×, lo-x 10–40× smaller than the default, MPP `runs/exw`) for Ex6 seeds not yet started.
   It needs nlo31 from 03f6981 or later (rebuild bin6 if older; default bit-identical, check as in
   step 1). Ex6 seeds already running or done with the default target are valid: let them finish, and
   keep them in separate output directories (`ex6_default/`, `ex6_vtm1/`) so they can be weighted
   separately.
   B/E commands of bin5 (`... 2 0 0 psmc 0 0 1d-12`, r with `VEGAS_EQUAL=2`): lo-x 200, vi-x 250,
   kp-x 200, r-x 2000 seeds (≈ 2.6k core-h).
3. **G6 — fixed-point NNLO at x = 0.1, Q² = 5000 (hard kinematics), photon and W⁻**, plain slicing,
   mode 1 (τ_zQ bins), `EW31="0 0"` and `EW31="2 0"`:
   - `sliced21 b1 3000000 6 SEED 0.1 5000 <softtable>`, `sliced21 b2 5000000 6 SEED 0.1 5000 <softtable>`;
   - `nlo31 lo|vi|kp|r NCALL 6 SEED 1d-9 1 0.1 5000 psmc 0 0 1d-12` (NCALL as B: lo 3M, vi 1M, kp 2M,
     r 1M; r with `VEGAS_EQUAL=2`); rcorr and d2 as for B (TECHDIFF), with `1 0.1 5000`;
   - seeds per boson: b1 100, lo 100, b2 100, vi 250, kp 200, r 2000, rcorr 300, d2 150;
   - pilots (4 seeds per part) for the CPU estimate; rough guess ≈ 3k core-h per boson.
   - NLO references (DISENT, `tau2_nlo -x 0.1 -Q2 5000 -nev 3000000`, 40 × 3M): W⁻ (with `-CC`) from
     MPP `runs/ccw1` per bin [NLO coefficient, x dσ/dx dQ² units as combine_tcut]: 0.05–0.1 7.606e-4
     ± 3.1e-6, 0.1–0.2 5.873e-4 ± 1.4e-6, 0.2–0.3 2.374e-4 ± 6.4e-7, 0.3–0.4 1.212e-4 ± 3.9e-7,
     0.4–0.5 7.177e-5 ± 2.7e-7, sum 1.7783e-3 ± 2.8e-6; photon: MPP runs it tonight (`runs/fph`), or
     run 40 seeds of `tau2_nlo` without `-CC` yourself.
4. Budget ≈ 9k core-h for Ex6 + G6; submit if the pilot estimate is within +30%, otherwise report first.
   Order: G6 W, Ex6, G6 photon. Watch the /ptmp file count (11 Oct incident): no NNLOJET this time.
5. Report as REPORT-2026-10-12 (items 6–8 format):
   - G6: per bin and τ_cut every part, NLO against DISENT, NNLO and its change per unit ln τ_cut in % of
     the bin's NLO (W and photon separately) — the question is whether the W rise appears at a fixed
     point in τ_zQ bins (then it is in the 2+1 pieces, not in the jet selection or the Q²/y integration);
   - Ex6: E NNLO P2B (bin5 plain + Ex6) per bin and total against NNLOJET (59.49 ± 12 pb), with the
     slopes, and the flagged-seed statistics of the extras with VTARGET=0.

## UPDATE 10 Oct (after the 9 Oct B round): three NNLO validations against NNLOJET

AK gave free rein over the weekend (9 Oct). When the bin4 B round of
UPDATE 9 Oct is complete and reported, prepare and run these three sets. All
the code is on `origin/2026-10-nnlo21-ncc` (it contains this branch, the
psmc fix, photon + Z, W exchange, the τ_zQ event-shape option `ZFIX=3` and
P2B for mode 1). Notebook 8–10 Oct for the local checks.

| set | NNLOJET | ours | local checks |
|---|---|---|---|
| **Z**: ZEUS dijets with photon + Z (e⁻) | `nnlojet_epLJJ_zeus2jZ.run` (V_NC = Z+GAMMA), existing build | `EW31="1 0"`, B commands, 5-class tables | LO per Q² bin within NNLOJET's errors; NLO 2+1 with Z = DISENT (fixed point) |
| **W**: ZEUS dijets with W⁻ (e⁻ p → ν) | `nnlojet_epNJJ_zeus2jW.run` (DISWm, identity CKM), existing build | `EW31="2 0"`, B commands, W tables | LO total 3.7554 ± 0.0011 vs 3.7508 ± 0.0028 pb; NLO 2+1 with W vs DISENT (x = 0.1, Q² = 5000, τ_zQ bins): total bin e⁻ −0.40 ± 0.24 % (24 b1 + 40 lo seeds), ν −0.46 ± 0.27 % at τ_cut 2e-4, every bin within 1.6σ; one-loop four-quark W functions = NNLOJET's pointwise |
| **E**: event shape τ_zQ in the inclusive cuts | `nnlojet_epLJJ_tauzq.run`, **patched build** (local observable `dis_tauzq`) | photon, `ZFIX=3`, plain slicing (+ P2B-extra runs) | LO 969.46 ± 0.93 vs 968.87 ± 0.55 pb; NLO vs NNLOJET (V + R, 8 + 24 seeds): total 248.9 ± 2.0 (ours, plain, τ_cut 1e-4) vs 249.7 ± 1.1 pb, every bin within 1.1σ |

1. `git fetch`; build `bin5` from `origin/2026-10-nnlo21-ncc` (new
   directory; `nnlo21/build_sliced21.sh`): sliced21, nlo31. Check: photon
   b0 and lo (mode 2, short runs) bit-identical to bin4.
2. Beam tables (5 classes; the old 3-class ones are refused):
   - photon/Z: `sliced21 mktab i i $P/tables/zeus5` for i = 0 … 28;
   - W: `EW31="2 0" sliced21 mktab i i $P/tables/zeus5w` (b is not in the
     down-type classes for W; these tables carry a flag word, and a W run
     refuses photon/Z tables and vice versa);
   - `sliced21 tabchkg <prefix>` (with `EW31="2 0"` for zeus5w): largest
     deviation ≲ 1e-2 of the largest coefficient (thA371a: ≤ 7.7e-3).
3. NNLOJET: Z and W with the existing build; E with a patched copy: in a
   copy of the NNLOJET v1.0.2 source `patch -p1 <
   nnlo21/validation/nnlojet-v1.0.2-dis_tauzq.patch`, its own build and
   install directory (do not touch the production build; hard21 links
   `libnnlojet_core.so`). Workflow and sizing as the photon ZEUS production
   (warmups per channel, LO R V RR RV VV). E has no jet requirement (only
   τ_zQ ≥ 0.05): check the warmups converge before sizing.
4. Our side with bin5, the B commands and seeds (ncall, itmx,
   `VEGAS_EQUAL=2` for r, the TECHDIFF rcorr set):
   - Z: `EW31="1 0"`, `zeus` with `$P/tables/zeus5`, `P2BSLICE=1` for
     lo/vi/kp/r;
   - W: `EW31="2 0"`, `zeus` with `$P/tables/zeus5w`, `P2BSLICE=1`;
   - E: `ZFIX=3`, `zeus` with `$P/tables/zeus5` (the Q², y cuts only; τ_zQ
     bins replace the ZEUS selection), photon, **plain** slicing for b0 … r;
     in addition lo/vi/kp/r with `P2BSLICE=1 P2BEXTRA=1` (the P2B − plain
     term alone) on about a quarter of the seeds.
   The matrix elements with Z and W are ≈ 1.6× slower than the photon's
   (more channels per point); E costs about the B set.
5. Pilots first, then the CPU estimate per set in NOTES. Budget for all
   three (NNLOJET and ours): ≈ 70k core-h. Submit if the estimate is
   within +30% of that; otherwise report first. Order: E, Z, W.
6. Report, as for B: per observable bin and τ_cut against NNLOJET (Z, W:
   total, Q², p̄_T, m12; E: the five τ_zQ bins and their sum, plain and
   P2B), plain means, flagged seeds, tail statistics. For E add, per bin,
   the value at the largest τ_cut with τ_cut ≤ 5e-3 × (lower bin edge) and
   τ_cut ≥ 1e-4, with the pulls against the smaller τ_cut using
   √|σ_a² − σ_b²| (the τ_cut points share events), and the sum of these
   per-bin values (notebook 10 Oct: at the fixed point each τ_zQ bin is flat
   in that window; the total only once the lowest bin is).
7. (11 Oct) Technical cut below τ_cut 1e-4: the residual is a series of
   decades in s_min (notebook 11 Oct); the second decade is not negligible
   (ZEUS: +2.3 ± 1.3 pb at 3e-5, +9 ± 2 pb at 1e-5). For each set run, besides
   rcorr, the second decade `TECHDIFF=1 TECHMIN="1d-11 1d-9" TECHREF="1d-10
   1d-8" VEGAS_EQUAL=2` (about 200 seeds, as conv) and add it to r for
   τ_cut < 1e-4; the third decade is ≈ 0 (median), do not go below 1e-12 W².
8. (11 Oct) Report additionally, per bin: each part (b2, vi, kp, r, rcorr,
   d2) per τ_cut, and the NNLO change per unit ln τ_cut in the windows
   5e-3→5e-4, 2e-3→2e-4, 1e-3→1e-4 in % of the bin's LO (pull errors
   √|σ_a² − σ_b²|). Bin 4 shows an upward drift of ≈ 1% of LO per unit ln in
   the ZEUS total (3–4.6% in the hard bins, ≈ 3σ); E (plain, τ_zQ) against Z
   (P2B, jets) and the E P2B-extra runs tell where it comes from.

## UPDATE 9 Oct: psmc fix, rerun the B set (AK: go)

Thanks for the 8 Oct (8b) report. Its replays (seeds 2042, 2049) found the
cause of the 1e5–1e6 tail, but it is not "incomplete cancellation in deep
corners": in all 81 DBG BIG events the real is cut away (τ₂ ≈ 4e-9…3e-8)
and an accepted dipole has a mapped Born with a **negative-energy parton**
(e.g. E = 183, 173.7, −31.8 GeV in the DBG `dip` lines). psmc built the
second emission with m²/E² up to 1e-6, so tiny p_i·p_j came out negative,
y > 1 in the dipole mapping, and the P2B fallback to plain weights accepted
the garbage dipole above τ_cut. Fixed in `dis31/psmc.f90` (emission in quad
precision; commit "psmc: emission in quad precision …"). On the thservs,
3 seeds × 3M points of B r: such accepted dipoles 453/568/809 per job before,
0 after. Details: notebook, 9 Oct.

The fixed point (F) had none of these (0 before and after); the fix does
not explain its non-flatness, which we study at MPP first. **Do not rerun F.**

1. `git pull`, rebuild nlo31/sliced21 as `bin4` (new directory; leave bin3).
   Checks before submitting:
   - sliced21 (b0, b1, b2) does not use the changed routines: a short run
     must be bit-identical to bin3;
   - nlo31 (lo, vi, kp, r) changes point by point (the transverse basis is
     rotated): a short lo run must agree with bin3 statistically. The
     existing lo, vi and kp results stay valid (no dipoles, no rerun).
2. Rerun with bin4, `VEGAS_EQUAL=2`, the same seeds and otherwise the same
   settings as the bin3 round:
   - B r (3,500), B rcorr (400);
   - Bx: conv (200), edge r (300), plain r (300), plain rcorr (100).
   From your bin3 timings ≈ 6,050 core-h; submit without asking again if
   the estimate is within 20% of that, otherwise report first.
3. Report, as for bin3 (same tables, so that bin3 and bin4 sit side by side):
   - flagged seeds (20× rule), the garbage-dipole counts, the "plain
     weights" counts per seed;
   - mean, median, trimmed mean, robust against sample σ: is the tail gone?
   - B NNLO per bin and τ_cut against NNLOJET (plain mean = result), down to
     1e-5;
   - if any seed still blows up: replay it with `P2BDEBUG=1`, keep the `DBG`
     lines, and check the energies in the `dip` lines.
4. Commit to `2026-10-nnlo21-cluster`.

## UPDATE 8 Oct: garbage-dipole guard, then rerun r with it

Thanks for the 8 Oct report. Your replays (seeds 3197, 4637) show what the
blow-ups are. Some dipoles have a mapped Born with an almost-zero-energy
parton (E ≈ 6e-4 GeV, finite and positive) and values of 1e61…1e179, next
to a real of ≈ 4e13 (the legitimate large terms). New in `nlo31`: an event
with a dipole |D| > 1e30 and |D| > 1e12 |R|, or a non-finite D, is
dropped with all its dipoles, in every mode, and counted ("r: events
dropped (garbage dipole) N of M" at the end of the output; thresholds via
`DIPGARB="1d30 1d12"`). On thA371a the output is bit-identical on normal
runs: the ≈ 3e-4 of events it counts there were already dropped as
non-finite.

1. `git pull` (at least the commit that added this section), rebuild
   nlo31/sliced21 (`bin3`). Check that a short b0/lo/r run (mode 2) is
   bit-identical to bin2 apart from the new count line.
2. Rerun with bin3, `VEGAS_EQUAL=2`, same seeds and otherwise the same
   settings:
   - B r (3,500), F r (1,200);
   - all TECHDIFF sets (B rcorr 400, F rcorr 300, plain rcorr, cut
     convergence 192);
   - the r cross-checks (edge, plain).
   CPU ≈ the last round. Report the estimate before submitting (AK).
3. Report, per set:
   - the number of flagged seeds (20× rule) before (bin2) and after (bin3),
     and the per-seed garbage counts;
   - mean, median, trimmed mean, and robust against sample σ (your 54 vs
     555 pb at 2e-4): does the heavy tail go away?
   - B: the NNLO per bin and τ_cut against NNLOJET (plain mean = result);
   - F: the NNLO against τ_cut, whether it becomes flat below 2e-4. It
     falls from 11.6 to 6.1 between 2e-3 and 1e-4 with equal weights: power
     corrections or something else?
   - if any seed still blows up: replay it with `P2BDEBUG=1` and keep the
     `DBG` lines.
4. Commit to `2026-10-nnlo21-cluster`.

## UPDATE 7 Oct evening (read this first): rerun r with equal iteration weights

**Cause of the ZEUS NNLO excess.** By default `nlo31` combines the histogram
cells over iterations with weights 1/σ²_it of the VEGAS target. For the
heavy-tailed part r (real − dipoles) this is biased upwards. In the target
cell (total, τ_cut 1e-5), the MPP P2B production had an equal-weight
iteration average of −8620.4 ± 8.5 against the reported −8570.6: a bias of
**+49.8 ± 7.1 pb**. The other parts (b0, b1, b2, lo, vi, kp) have
|bias| ≤ 0.5 pb.

With r rerun on thA371a with `VEGAS_EQUAL=2` (291 seeds, otherwise the MPP
settings), plus the existing b2, vi, kp and the TECHDIFF correction, the
ZEUS NNLO total is 10.6 ± 1.5, 12.6, 11.2, 11.2, 12.3, 11.8 ± 4.4,
13.7 ± 4.2, 14.6 ± 4.7 pb from τ_cut 2e-2 to 1e-4. NNLOJET: 11.12 ± 1.19.
At 3e-5 and 1e-5: −1.7 ± 7.5 and −8.5 ± 9.2, low by 1.7–2σ. There the
mean of r is ≈ 2σ below its median: heavy-tailed outliers, which now count
fully. A few high-p̄_T/high-m₁₂ bins still show excesses of 2–4σ at
2e-4…1e-4 (p̄_T 30–60: 0.8 ± 0.2 against −0.18 ± 0.04). The mode-2
window against the fixed point now agrees within errors.

What to do:
1. `git pull` (branch `2026-10-nnlo21`, at least the commit that added
   this section) and rebuild sliced21/nlo31. The new commits fix the
   `kt_jets` NaN crash (non-finite momenta are rejected; seed 4711
   completes) and the `hcacc` bounds. They also add `ZFIX=2`, the
   diagnostic window with the ZEUS selection.
2. First, with existing outputs: for B r, F r and plain r, compute per
   seed the equal-weight average of iterations 2–6 (the ` iteration` lines
   are the target cell: total at 1e-5 in mode 2, all τ_zQ bins at 1e-5 in
   mode 1). Compare it with the reported cell. Expect ≈ +50 pb for B r.
   Report it.
3. Rerun every r with `VEGAS_EQUAL=2` and otherwise unchanged settings
   (seeds may be reused):
   - B r, P2BSLICE = 1;
   - F r;
   - the B cross-checks with r: psmc edge 1e-13 (the measured "edge
     shift" may itself be a weighting artefact), plain slicing, cut
     convergence.
   The TECHDIFF correction runs already used `VEGAS_EQUAL=2`; keep them.
   b0, b1, b2, lo, vi, kp and NNLOJET need no rerun.
4. Outliers. With equal weights a garbage iteration counts fully. Per part
   and column, report plain mean, median and 1% trimmed mean, and the seeds
   with any iteration beyond 20× the median |iteration| of that part. The
   plain mean stays the result (AK); list the flagged seeds. Replay one or
   two flagged seeds with `P2BDEBUG=1` and keep the `DBG BIG` lines for us.
5. Report: the NNLO coefficient per bin and τ_cut against NNLOJET (B), the
   fixed-point NNLO (F), the cross-checks. Commit to `2026-10-nnlo21-cluster`.
   Push only with AK's OK.

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
