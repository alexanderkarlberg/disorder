# NNLO DIS 2+1 by τ₂ slicing: cluster production — instructions for Claude

**Written 6 Oct 2026 on thA371a (draft of the night, updated in the evening
with projected slicing and the new technical cut).** For the Claude session that AK starts on the cluster (Slurm, partition
`alma`, 24 h wall time, at most 10k running and 30k queued jobs). Read this
file, `CLAUDE.md`, `docs/nnlo21-plan.md` and the nnlo21 entries of
`docs/notebook.md` (3–6 Oct) first. The page with the current status:
https://claude.ai/artifact/ChKctwxDdEau7H9DGMxwyW

## UPDATE 16 Oct (b): KREDO — rerun every kp and kp-x with the fixed nlo31 (AK: go)

A bug in nlo31's K + P was found at MPP (notebook 16 Oct, commit 40c013a): in `kp_conv` a local variable `Cf`
shadowed the colour factor `CF` (Fortran is case-insensitive), so the q → q integrated plus-distribution term
was wrong since the first nlo31 (2 Oct). **Only the parts `kp` and `kp-x` are affected** (every mode); b0, b1,
b2, lo, vi, r, rcorr, d2, lo-x, vi-x, r-x and all NNLOJET runs stay as they are. At MPP (same seeds, old and
fixed binary) the fix removes about half of the jet drift (J8 42–55%, J7 PH 47–54%, W5 56–67%); a smaller drift
remains in the quark channel, so **Q8 continues unchanged** (it runs only r, rcorr, d2).
1. Build **bin11** from `origin/2026-10-nnlo21-ncc` at this commit or later. Against bin10 only nlo31's kp
   changes (all other parts, and sliced21, bit-identical: checked at MPP on a J8 kp run, new build = the fixed
   binary of the impact runs). Checks:
   - `nlo31 chk 1 1` ends with three lines `kp_conv q -> q (b, g, lp) / closed form - 1, xi ...`, all
     |values| ≤ 1e-7 (bin10 gives 2e-3 … 1 in the b and lp columns);
   - one J8 kp seed with bin11 and bin10 (same command): Δ = bin11 − bin10 of the total row at τ_cut 1e-3,
     1e-4, 1e-5 should be close to MPP's mean −7.41e-8, −1.276e-7, −1.650e-7 pb/GeV² (per-seed scatter ≈ 1%).
2. Rerun **every kp and kp-x job that enters a reported result** with bin11: the same commands, environment,
   ncall, itmx and seeds as the original jobs; only the binary changes. New directories next to the old ones
   (e.g. `kp11`, `kp-x11`); keep the old outputs. Order:
   - (a) CH8 q and g, J8 (Q8's report needs the CH8 q kp);
   - (b) J7 PH, J7 W, J7v;
   - (c) the bin5 validations Z, W, E (and E's kp-x), and B (bin4, with its extras);
   - (d) G6 photon and W; Ex6, Ex7, X8 (kp-x).
   Not needed: G8 (flat and psmc kp carry the same bug term, so the comparison stands), D (NNLO 1+1, no kp),
   N7 (NNLOJET). The local N3LO 1+1 P2B set (MPP `runs/n3p`) is rerun at MPP.
   Cost: ≈ 3.5k jobs of ≈ 0.1 core-h (X8 kp-x 0.074): a few hundred core-h. Submit if ≤ 2k core-h, otherwise
   report the estimate first. Q8 has priority if the two compete.
3. Recombine each set with the bin11 kp / kp-x as in its last report (NNLO per bin at every τ_cut, slopes,
   pulls against NNLOJET for Z, W, E and the N7 references). Per set also Δkp = bin11 − bin10 per cell with its
   slope per unit ln τ_cut; compare with MPP's same-seed values (total row, τ_cut 1e-3 / 1e-4 / 1e-5):
   J8 −7.41e-8 / −1.276e-7 / −1.650e-7, J7 PH −9.21e-7 / −1.532e-6 / −2.190e-6 (pb/GeV²), W5 = bin5 W
   −0.1043 / −0.1683 / −0.2029 pb. Agreement within errors shows the reruns used the right binary and commands.
4. Report (REPORT-…-KREDO, as soon as (a) and (b) are done; (c), (d) may follow in a second report): per set the
   corrected tables next to the bin10 ones, and the Δkp checks of item 3. In the Q8 report use the bin11 CH8 q kp.
5. Order: Q8 and KREDO (a) now, then KREDO (b)–(d).

## UPDATE 16 Oct: Q8 — does the quark-channel drift stop below τ_cut 1e-4? (AK: go)

REPORT-2026-10-15-CH8 received (MPP), thanks: the drift is in the quark channel (gluon share 0.03 ± 0.03).
In every set the last window (1e-4 → 1e-5) is consistent with zero (CH8 q −11.3 ± 11.7, J8 +1.9 ± 13.1, against
≈ −44 and −55 from 1e-3 → 1e-4). That decides between two readings: a wrong ln τ_cut term (the slope continues)
or large quark-channel power corrections that die out by 1e-4 (the slope stops, and the discrepancy with NNLOJET
is a constant offset at τ_cut → 0). Q8 measures the quark-channel slope 1e-4 → 1e-5 (and 3e-5) precisely.
1. bin10, CH8 q set-up exactly (mask 1: `zeus5q`, sliced21 arguments `0.02 1`, nlo31 argument 10 = `1`; J8
   window and commands). **Only the parts that carry the error below 1e-4: r, rcorr and d2** (and b2, vi, kp only
   if your budget says they matter; at 1e-5 in J8 the errors were r 2.7e-8, rcorr 1.6e-8, d2 8.5e-9, b2 0.9e-9,
   vi 2.0e-9, kp 0.4e-9 pb/GeV²).
2. Sizing: from the CH8 q seeds, per part, the per-seed sd of N(1e-4) − N(1e-5) (and of N(1e-4) − N(3e-5)); choose
   the extra seeds so that the error of the total-row NNLO slope 1e-4 → 1e-5, combined with CH8 q, is **≤ 6e-9
   pb/GeV² per unit ln** (now 11.7e-9; the two readings differ by ≈ 40e-9). CPU per job as in CH8 (r 2.03, rcorr
   0.61, d2 0.60 core-h). Submit if ≤ 15k core-h, otherwise report the estimate first.
3. Seeds: disjoint from everything of CH8 (CH8 q used 7001–8600 for r, g 8001–9600: the two channels share
   8001–8600, harmless since the integrands differ, but avoid it now): 21001– for every part (five digits, the
   four-digit range is used up).
4. Report: combined with CH8 q (equal-weight seed means over old + new seeds), the quark-channel NNLO total,
   p̄_T 8–15 and both Q² bins at 2e-4, 1e-4, 3e-5, 1e-5, the slopes 2e-4 → 2e-5-like windows we have (2e-4 → 1e-4,
   1e-4 → 3e-5, 1e-4 → 1e-5, 3e-5 → 1e-5), and the per-part decomposition of the 1e-4 → 1e-5 change. With the
   rcorr/d2 contributions shown separately (the technical-cut decades grow fastest at small τ_cut).
5. Order: now.

## UPDATE 15 Oct (b): CH8 — the J8 window split into quark- and gluon-initiated channels (AK: go)

REPORT-2026-10-15-G8 received (MPP), thanks; the correction of the sizing is noted. Reading at MPP: where flat
is reliable (5e-3 → 5e-4) it drifts like psmc; no sign of a psmc gap, and no further generator runs.
CH8 splits the drift by the incoming parton: quarks + antiquarks (mask 1) and gluon (mask 2). Every part uses the
same mask (Born luminosity, beam functions, above-cut PDFs), so the two channels add up exactly to J8.
1. Build **bin10** from `origin/2026-10-nnlo21-ncc` at this commit or later. New: `sliced21 mktab` and `tabchkg`
   take `PDFMASK=1|2` from the environment (beam tables of one channel; default unchanged). Checked at MPP
   (node 10): unmasked bit-identical to the production table, mask 1 + mask 2 = unmasked to 9e-11.
2. Tables: build `$P/tables/zeus5q` (`PDFMASK=1`) and `$P/tables/zeus5g` (`PDFMASK=2`) like `zeus5` (29 nodes,
   photon); run `tabchkg` on each with the same PDFMASK.
3. Commands: those of J8 (`ZFIX=2 ZCUTS=1 EW31="0 0"`, window (0.1, 5000), `ZWIN="0.1 5000"` for sliced21,
   psmc), with the channel mask in every part:
   - sliced21 (b2, b1): tables `zeus5q` / `zeus5g` and arguments 8, 9 = `0.02 <mask>` (argument 8 is not used
     in mode 2 but must be there): `sliced21 b2 5000000 6 SEED zeus $P/tables/zeus5q $P/tables/soft.tab 0.02 1`;
   - nlo31 (vi, kp, r, rcorr, d2, lo): argument 10 = mask, e.g. `nlo31 r 1000000 6 SEED 1d-9 2 0.1 5000 psmc 1 0 1d-12`.
   Per channel the J8 seed counts: b2 200, vi 240, kp 200, r 1600, rcorr 200, d2 100, b1 100, lo 100 (seeds
   7001– for mask 1, 8001– for mask 2). About 2 × 3k core-h; pilot 5 per part, submit if within +30% of J8.
4. Report: per channel and their sum next to J8 (total, both Q² bins, p̄_T and m12 bins): NLO, NNLO at every τ_cut,
   slopes as in J8 (absolute and in % of the channel's |NLO(1e-5)|), and the pulls (q + g) − J8 per cell (must be
   consistent). The question: does the drift sit in one channel?
5. Order: CH8 now (G8 is finished).

## UPDATE 15 Oct: G8 — the J8 window with nlo31's other phase-space generator (AK: go)

REPORTs J7v and J8 received (MPP, notebook 15 Oct), thanks. J7v excludes a statistical VEGAS bias (parts
agree cell by cell, near-Gaussian per-seed values); the cut-design rows show the drift sits entirely in the
Born-level jet kinematics. Open: a coverage gap of the psmc phase-space generator near the 2+1 limit (a region
psmc generates badly or not at all gives a deficit that grows as τ_cut falls, whatever VEGAS does). G8 reruns
the above-cut parts in the J8 window, where the drift is 12–16σ, with the generator of nlo31's argument 9
set to `flat` instead of `psmc` (Breit-frame lepton + `breit` map, the generator before 5 Oct).
1. bin9; the J8 commands for **vi, kp, r** (`ZFIX=2 ZCUTS=1 EW31="0 0" P2BSLICE=1`, r with
   `VEGAS_EQUAL=2`, window `2 0.1 5000`) with argument 9 `flat`:
   `nlo31 <part> <ncall> 6 SEED 1d-9 2 0.1 5000 flat 0 0 1d-12`. Seeds 6001–.
   b2, rcorr, d2, b1, lo: J8's as they are (b2 does not use the generator; rcorr and d2 are ≈ 1/10 of the
   drift there, keep J8's and say so).
2. Pilot 5 per part first. At MPP (short runs, J8 window, 200k × 4 for r): `flat` r is ≈ 3× slower per call
   than psmc, with ≈ 3× the per-iteration error before adaptation. Size from the pilot: N per part such that
   the error of the total NNLO slope 1e-3 → 1e-4 is ≤ 18e-9 pb/GeV² per unit ln (a third of the J8 drift,
   −54.9e-9 ± 4.6e-9). Submit if ≤ 6k core-h, otherwise report the estimate first.
3. Report: per part (vi, kp, r) and cell (total, p̄_T 8–15, at every τ_cut) the means for `flat` and J8
   (psmc) with pulls; per-seed sd, 1.48 MAD, median, fraction beyond 20 MAD; the NNLO total and p̄_T 8–15 and
   their slopes with the `flat` parts next to J8. The question: does the drift change with the generator?
4. After J7v is complete (its remaining seeds).
5. Caveat from the MPP short runs (2 flat seeds, vi 100k × 4, r 200k × 4, total row): flat agrees with psmc at
   τ_cut ≥ 2e-3 but is 15–20% low in vi at 1e-3 and below (vi at 3e-5: −8.4e-6, −7.0e-6 against −10.7e-6 ±
   0.1e-6 with psmc), with RESULT errors that do not cover the difference. Most likely flat undersamples small τ_2
   (no log mapping), which would mimic a coverage gap in psmc. So, in addition: 5 flat seeds per part at
   **4× ncall** in the pilot; if their means move towards psmc, flat is undersampling there, and only the
   τ_cut where flat is stable under the 4× test count for the comparison (say which). Report per-seed medians
   next to the means.

## UPDATE 14 Oct (d): cut-design rows `ZCUTS=1` attached to J8 and J7v (AK: go)

Why: the literature (docs/literature-fiducial-power-corrections.md, 14 Oct) names √τ·L^k power corrections
from the symmetric per-jet E_T cut as a candidate for the drift; their slope in ln τ_cut is almost constant
over our range. `ZCUTS=1` (nlo31 and sliced21, mode 2, also with ZFIX=2) adds six rows after the 15 ZEUS
rows: `ZCELL et4`, `noeta`, `et4ne`, each as total (`0 0`) and p̄_T 8–15. et4: jets E_T > 4, p̄_T > 8 (η and
m12 cuts as ZEUS); noeta: ZEUS E_T cut, no η_lab cuts; et4ne: both. Without ZCUTS the output is
bit-identical; with it the ZEUS rows and RESULT are unchanged (checked at MPP, notebook 14 Oct).
1. Build **bin9** from `origin/2026-10-nnlo21-ncc` at this commit or later (it has `VTARGET`, `ZWIN` and
   `ZCUTS`). Check: a short sliced21 b0 with `ZFIX=2 ZCUTS=1` prints 21 ZCELL rows, with `et4` equal to
   the ZEUS rows (balanced Born jets).
2. **J8**: use bin9 instead of bin8 and add `ZCUTS=1` to every job (b2, vi, kp, r, rcorr, d2, b1, lo). If J8
   jobs are already running on bin8, stop feeding them, keep what is done for the ZEUS rows, and run the
   remaining seeds on bin9 with ZCUTS=1; report the extra rows from the bin9 seeds only.
3. **J7v**: bin9 instead of bin7, with `ZCUTS=1` (vi, kp, r as in UPDATE (c)). The extra rows also need
   the other parts in the window: rerun **J7 PH b2 200, rcorr 200, d2 100, b1 100, lo 100** with bin9,
   `ZCUTS=1` and the J7 PH commands, seeds 3001– (≈ 400 core-h). The ZEUS rows of J7v can still be combined
   with J7 PH's old parts as in (c).
4. Report per set (J8, J7v): for the three variants, the NNLO total and p̄_T 8–15 at every τ_cut, and the
   slopes as in J7, next to the ZEUS rows. The question: does the drift vanish or shrink in et4?
5. Order unchanged: J8, J7v (with the reruns of 3), X8.

## UPDATE 14 Oct (c): J7v — J7 photon above the cut with VEGAS on all cells (AK: go)

Why: a pointwise RG check at MPP (notebook 14 Oct) found no error below the cut. The drift is ≈ 0.1% of
the b2 ↔ (vi + kp + r) cancellation (J7 PH total at 1e-4: b2 1.1e-3, NNLO 6e-6), so a sampling bias of
0.1–0.2% in r (or vi, kp) would explain it. J7 aimed VEGAS at one cell (total row at 1e-5), and vi, kp
combined iterations with 1/σ² weights; Ex7 showed that cells VEGAS does not adapt to can be far off.
1. bin7 (it has `VTARGET`), window `2 0.01 400`, `ZFIX=2 EW31="0 0"`, the J7 PH commands for vi, kp, r
   only, each with `P2BSLICE=1 VTARGET=-1 VEGAS_EQUAL=2` (r already had VEGAS_EQUAL=2). Seeds 3001–
   (J7 PH used 1001–): vi 240, kp 200, r 1600. Pilot 5 per part; submit if the CPU per job is within
   +50% of J7 PH's (VTARGET=-1 changes the sampling, not the work per point).
2. Combine with J7 PH's b2, rcorr, d2 (and b1, lo for the NLO) as they are.
3. Report: per bin (as J7) the NNLO and the slopes for J7v next to J7 PH; per part and cell the pulls
   J7v − J7 PH (as `Ex7_7_consistency.txt`); per-seed sd, 1.48 MAD, median and the fraction beyond
   20 MAD for vi, kp, r in both sets (total row and ptavg 8–15, at 1e-3 and 1e-4).
4. Order: J8, then J7v, then X8.

## UPDATE 14 Oct (b): J8 — ZEUS jet cuts in a window at the G6 point (x, Q²) = (0.1, 5000), photon (AK: go)

Why: the drift is there for jets at (0.01, 400) (J7) and in the ZEUS set, while the G6 fixed point (0.1, 5000)
without jets is flat in the total. J8 puts the ZEUS jet cuts at the G6 point. If the jets drift there too, the
cause comes with the jet selection; if they are flat, it depends on x (gluon content at x = 0.01).
1. Build **bin8** from `origin/2026-10-nnlo21-ncc` (this commit or later): sliced21 now takes the window
   centre from `ZWIN="x Q2"` (default (0.01, 400), bit-identical; checked at MPP against the zw2 binary).
   Check: one short sliced21 b1 run with `ZFIX=2 ZWIN="0.1 5000"` prints `x 0.100000 Q2 5000.00`, and one
   short nlo31 lo run with `2 0.1 5000` prints `ZFIX window Q2, y:  4.997500E+03  5.002500E+03 ...`.
2. Commands: those of J7 PH (`EW31="0 0"`, tables `$P/tables/zeus5`, `ZFIX=2`) with bin8, plus
   `ZWIN="0.1 5000"` for every sliced21 job and the window arguments `2 0.1 5000` (instead of `2 0.01 400`)
   for every nlo31 job. Seed counts as J7 PH, including the NLO parts you added (b1 100, lo 100): b2 200,
   vi 240, kp 200, r 1600, rcorr 200, d2 100; seeds 1001–. Pilot 5 per part, submit if within +30% of the
   J7 PH cost per part (the cross section is smaller here; the CPU per job need not be).
3. The window straddles the Q² bin edge 5000, so both `q2 2000-5000` and `q2 5000-20000` are filled; use the
   `total` row as the main number.
4. Report (as soon as it is done): the J7 tables and slopes for J8, side by side with J7 PH, per bin; in % of
   |NLO(1e-5)| as before.
5. X8 of UPDATE 14 Oct goes after J8 (X8 is small).

## UPDATE 14 Oct (after REPORT-2026-10-14): hold production; one Ex7 check (X8)

REPORT-2026-10-14 received (MPP), thanks; the J7 NLO addition (b1, lo) was the right call. N7 settles the
jet drift: NNLOJET W 0.148 ± 0.015 pb, ours 4–5σ above at τ_cut ≤ 2e-4 and still rising, and J7 shows the
photon drifting like W. The cause is in our code and is hunted at MPP. Until then:
1. **No new production of ours** (sliced21/nlo31) and no further NNLOJET top-ups. Keep all outputs.
2. **X8, a check of Ex7 at large τ_cut** (REPORT §4.4: Ex7 vs Ex6/bin5 up to 29σ apart at τ_cut ≥ 2e-3).
   Our reading: bin5 (target cell 60 = total row at 1e-5) and Ex6 (`VTARGET=0`) hardly adapt to the large
   τ_cut cells, so their means there come from rare seeds and their errors are underestimated. X8 tests this
   with VEGAS aimed at the total row at τ_cut 5e-3: `VTARGET=53` (cell = column + 10 × (row − 1); total row
   0.05–0.5 = row 6, 5e-3 = column 3). Check in one short run that the log prints `VEGAS target cell   53`.
   - Ex7 commands and bin7 (`ZFIX=3 P2BSLICE=1 P2BEXTRA=1`, r-x with `VEGAS_EQUAL=2`), only `VTARGET=53`
     instead of −1; seeds 8001– (disjoint from Ex7 5001–7000, the pilots 9001–9300, Ex6 and bin5):
     lo-x 200, vi-x 250, kp-x 200, r-x 400. Estimate first (Ex7 CPU per job); submit if ≤ 1k core-h.
   - Report (short, e.g. REPORT-2026-10-15): per extra part, total row, τ_cut 2e-2 … 2e-4: X8, Ex7, Ex6, bin5
     means ± sd/√N, per-seed sd, 1.48 MAD, median, and the pulls X8 − Ex7, X8 − Ex6, X8 − bin5. At 5e-3
     (the target) X8 should be the most reliable of the four.
3. N7 Z: no action (the RRa_2 seed 11083 is noted; plain means stay).

## UPDATE 13 Oct (e): feed Ex7 production (after the 200-seed pilot, 0bb89c8)

The 200-seed pilot (`results/ex7_pilot200.txt`) meets all three conditions of (d). (c) asks for agreement
within the errors: pulls −1.43 (1e-3) and +0.52 (1e-4) are that. The Ex6 error is a tail-driven sd/√N, and
the medians differ only because the two distributions have different shapes (Ex6: most seeds near zero, rare
large seeds carry the mean), so neither is a reason to hold back. Ex7 has 3.5× (1e-3) and 5.5× (1e-4) smaller
errors than Ex6 with a tenth of the seeds.
1. Start the prepared Ex7 feeder (`groups7.txt`: lo-x 200, vi-x 250, kp-x 200, r-x 2000, seeds 5001–, bin7,
   `VTARGET=-1`), after J7 PH, N7 W and J7 W in priority (with or below N7 Z).
2. Report as in UPDATE 13 Oct, item 3: Ex7 alone (r-x 5001–7000 together with the pilot 9101–9300), and
   Ex6 + Ex7 combined per cell by inverse variance. Also give Ex7 alone with the bin5 E lo/vi/kp extras, as
   the inverse-variance weight of Ex6 is close to zero and its tails make its sd/√N unreliable.

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
