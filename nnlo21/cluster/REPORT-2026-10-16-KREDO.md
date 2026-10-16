# nnlo21 cluster: KREDO (a)+(b) -- kp rerun with the fixed nlo31 (bin11) -- report (16 Oct 2026)

Branch `2026-10-nnlo21-cluster`. UPDATE 16 Oct (b). Facts and significances only; no interpretation. (c) and (d) are also complete and follow in a second report.
Scripts (`nnlo21/cluster/`): `kredo.py` (reruns), `kpdelta.py` (dkp = bin11 - bin10, same seeds); the combine scripts of the earlier reports (`combine8j.py`, `combinech8.py`, `combine7j.py`, `combine7v.py`, and for (c)/(d) `combine5.py`, `combine6.py`, `combine7e.py`, `x8.py`) take `KP11=1`: they then read `kp11` / `kp-x11` instead of `kp` / `kp-x` and write to `results/kp11/` (same file names as the bin10 results in `results/`). Without `KP11` the outputs are unchanged (checked for `combine8j.py PH9`: every text output identical, raw dumps identical in content).
Numbers: `results/kp11/` (corrected tables of each set, `KPDELTA_<set>.txt`, `KPDELTA.json`; raw kp ZCELL rows of the bin11 seeds in `results/kp11/{j8,ch8,j7}/*_kp.zcell.gz`). Conventions of each set as in its last report (J8, J7: slopes in % of the set's bin |NLO(1e-5)|; CH8: also absolute in 1e-9 pb/GeV^2 per unit ln).

## 0. Summary

- **bin11** = 581e4fd: `nlo31 chk 1 1` kp_conv q -> q / closed form - 1: xi 0.01 -1.6e-8, -7.3e-9, -8.1e-9; xi 0.1 <= 2.1e-10; xi 0.5 <= 1.2e-11 (all <= 1e-7). Against bin10 (J8 window, `$P/runs/chk11/chk.sh`): sliced21 b0, nlo31 lo, nlo31 kp with mask 2 cmp-identical; kp with mask 1 differs. MD5 sliced21 e595fdd0997a29caecac51889b4fd998, nlo31 40b98b37e1934e570cee47ee67ce8739. (A first build from 40c013a, used for the CH8 q kp of the Q8 report, gives cmp-identical output.)
- **Reruns**: same cmd.sh as the original jobs with only the nlo31 path changed (checked by diff for every job), every done old seed. (a) J8 kp 200, CH8 q kp 200, CH8 g kp 200; (b) J7 PH kp 200, J7 W kp 150, J7v kp 200 (bin9 seeds) + 5 (bin7 seeds). All done, none failed. (c), (d): Z, W, E kp 200 each, E kp-x 50, the 4-seed pilots, B kp 200 + Bx plain-kp 16 (from the 11 Oct archive), G6 photon / W kp 200 each + pilots, Ex6, Ex7, X8 kp-x 200 each: all done. Not rerun: G8, G8x4 (item 2), J8's five bin8 kp seeds (`ours8j/PH/kp`, not in the J8 combination). Cost ~260 core-h.
- **dkp check (item 3)**, total row, tau_cut 1e-3 / 1e-4 / 1e-5, pb/GeV^2, ours (same seeds, N) vs MPP:

  | set | ours | MPP |
  |---|---|---|
  | J8 (200) | -7.403e-8 / -1.276e-7 / -1.654e-7 (+- 0.8e-10 ... 1.2e-10) | -7.41e-8 / -1.276e-7 / -1.650e-7 |
  | J7 PH (200) | -9.140e-7 / -1.533e-6 / -2.190e-6 (+- 2.5e-9 ... 3.6e-9) | -9.21e-7 / -1.532e-6 / -2.190e-6 |
  | CH8 q (200) | -7.400e-8 / -1.276e-7 / -1.653e-7 | (= J8) |
  | CH8 g (200) | 0 exactly at every tau_cut and cell | |
  | J7v bin9 seeds (200) | -9.152e-7 / -1.534e-6 / -2.190e-6 | |
  | J7 W (150) | -8.302e-9 / -1.386e-8 / -1.967e-8 | |

  J7 PH at 1e-3 differs from MPP's value by 0.8% (MPP's 24-seed per-seed scatter is ~1%).
  dkp slopes per unit ln (total row), windows 5e-3->5e-4 / 2e-3->2e-4 / 1e-3->1e-4 / 1e-4->1e-5: J8 +2.316e-8 / +2.371e-8 / +2.327e-8 / +1.639e-8 (MPP: +2.31, +2.36, +2.33e-8); J7 PH +2.403e-7 / +2.588e-7 / +2.686e-7 / +2.854e-7 (MPP: +2.44, +2.54, +2.65e-7).
- **NLO unchanged** in every set (b1, lo not rerun; checked: identical).

## 1. Corrected NNLO slopes, total row and main bins (bin10 -> bin11 kp)

### J8 (photon, ZEUS jets at (0.1, 5000)); % of the bin |NLO(1e-5)| per unit ln, (pull)

| bin | kp | 5e-3->5e-4 | 2e-3->2e-4 | 1e-3->1e-4 | 1e-4->1e-5 |
|---|---|---|---|---|---|
| total | bin10 | -44.11 +- 2.76 (-16.0) | -51.99 +- 3.72 (-14.0) | -57.67 +- 4.73 (-12.2) | +1.97 +- 13.64 (+0.1) |
| | **bin11** | **-19.95 +- 2.76 (-7.2)** | **-27.26 +- 3.72 (-7.3)** | **-33.38 +- 4.73 (-7.1)** | **+19.07 +- 13.64 (+1.4)** |
| Q^2 2000-5000 | bin10 | -45.17 +- 4.10 (-11.0) | -60.70 +- 5.50 (-11.0) | -68.32 +- 7.00 (-9.8) | +13.69 +- 19.15 (+0.7) |
| | bin11 | -20.65 +- 4.10 (-5.0) | -35.60 +- 5.50 (-6.5) | -43.58 +- 7.00 (-6.2) | +31.07 +- 19.15 (+1.6) |
| Q^2 5000-20000 | bin10 | -43.09 +- 3.75 (-11.5) | -43.53 +- 5.25 (-8.3) | -47.31 +- 6.77 (-7.0) | -9.43 +- 16.87 (-0.6) |
| | bin11 | -19.27 +- 3.75 (-5.1) | -19.15 +- 5.25 (-3.6) | -23.46 +- 6.77 (-3.5) | +7.40 +- 16.87 (+0.4) |
| ptavg 8-15 | bin10 | -28.67 +- 3.47 (-8.3) | -39.45 +- 4.88 (-8.1) | -47.60 +- 5.94 (-8.0) | +5.85 +- 15.94 (+0.4) |
| | bin11 | -12.05 +- 3.47 (-3.5) | -22.00 +- 4.88 (-4.5) | -30.14 +- 5.94 (-5.1) | +19.27 +- 15.94 (+1.2) |

Total row absolute (1e-9 pb/GeV^2 per unit ln): bin11 -19.12 +- 2.64, -26.12 +- 3.57, -31.99 +- 4.53, +18.28 +- 13.07 (bin10 -42.28, -49.83, -55.27, +1.89). All bins and the ZCUTS rows: `results/kp11/J8_PH9_slopes.txt`.

### CH8 (quark / gluon channel), total row: % of the channel's |NLO(1e-5)| | absolute 1e-9 per unit ln

| set | kp | 5e-3->5e-4 | 2e-3->2e-4 | 1e-3->1e-4 | 1e-4->1e-5 |
|---|---|---|---|---|---|
| q | bin10 | -65.9 +- 5.3 \| -33.73 +- 2.74 | -87.4 +- 7.4 \| -44.76 +- 3.77 | -85.2 +- 10.3 \| -43.65 +- 5.28 | -22.0 +- 22.9 \| -11.29 +- 11.73 |
| q | **bin11** | **-20.8 +- 5.3 \| -10.63 +- 2.74 (-3.9)** | **-41.2 +- 7.4 \| -21.09 +- 3.77 (-5.6)** | **-39.8 +- 10.3 \| -20.38 +- 5.28 (-3.9)** | **+9.9 +- 22.9 \| +5.09 +- 11.73 (+0.4)** |
| g | both | -1.7 +- 2.4 \| -0.70 +- 1.02 | -6.4 +- 3.3 \| -2.70 +- 1.41 | -3.6 +- 3.4 \| -1.53 +- 1.42 | +3.5 +- 6.2 \| +1.46 +- 2.61 |
| q+g | bin11 | -12.1 +- 3.1 \| -11.34 +- 2.92 (-3.9) | -25.4 +- 4.3 \| -23.79 +- 4.03 (-5.9) | -23.4 +- 5.8 \| -21.91 +- 5.47 (-4.0) | +7.0 +- 12.8 \| +6.55 +- 12.02 (+0.5) |
| J8 | bin11 | -19.9 +- 2.8 \| -19.12 +- 2.64 (-7.2) | -27.3 +- 3.7 \| -26.12 +- 3.57 (-7.3) | -33.4 +- 4.7 \| -31.99 +- 4.53 (-7.1) | +19.1 +- 13.6 \| +18.28 +- 13.07 (+1.4) |

The gluon channel is unchanged (dkp = 0). The (q+g) - J8 pulls are unchanged (340 pulls: mean -0.45, rms 0.89, max 2.6), since q and J8 get the same dkp. Gluon share of the q+g slope (1e-3 -> 1e-4) is now 0.07 +- 0.07. With the Q8 seeds (REPORT-2026-10-16-Q8): q 1e-3 -> 1e-4 -16.90 +- 2.98, 1e-4 -> 1e-5 +8.59 +- 7.00.

### J7 PH and J7 W (photon / W-, window (0.01, 400)); % of the bin |NLO(1e-5)| per unit ln, (pull)

| set, bin | kp | 5e-3->5e-4 | 2e-3->2e-4 | 1e-3->1e-4 | 1e-4->1e-5 |
|---|---|---|---|---|---|
| PH total | bin10 | -17.45 +- 6.14 (-2.8) | -17.65 +- 6.78 (-2.6) | -21.02 +- 7.07 (-3.0) | -26.25 +- 24.13 (-1.1) |
| | **bin11** | **-8.47 +- 6.14 (-1.4)** | **-7.98 +- 6.78 (-1.2)** | **-10.98 +- 7.07 (-1.6)** | **-15.59 +- 24.13 (-0.6)** |
| PH ptavg 8-15 | bin10 | -74.77 +- 27.29 (-2.7) | -87.41 +- 31.53 (-2.8) | -89.63 +- 34.82 (-2.6) | -77.33 +- 123.13 (-0.6) |
| | bin11 | -49.40 +- 27.29 (-1.8) | -59.77 +- 31.53 (-1.9) | -60.69 +- 34.82 (-1.7) | -46.00 +- 123.13 (-0.4) |
| W total | bin10 | -21.25 +- 7.94 (-2.7) | -20.35 +- 9.12 (-2.2) | -26.46 +- 9.68 (-2.7) | -1.89 +- 17.91 (-0.1) |
| | **bin11** | **-8.20 +- 7.94 (-1.0)** | **-6.40 +- 9.12 (-0.7)** | **-11.96 +- 9.68 (-1.2)** | **+13.26 +- 17.91 (+0.7)** |
| W ptavg 8-15 | bin10 | -92.52 +- 54.68 (-1.7) | -80.54 +- 65.69 (-1.2) | -121.54 +- 69.00 (-1.8) | -43.19 +- 109.66 (-0.4) |
| | bin11 | -45.93 +- 54.68 (-0.8) | -30.11 +- 65.69 (-0.5) | -68.41 +- 69.00 (-1.0) | +13.63 +- 109.66 (+0.1) |

### J7v (sets PH, Z, V9 as in REPORT-2026-10-15-J7v), total row, % per unit ln, (pull)

| set | kp | 5e-3->5e-4 | 2e-3->2e-4 | 1e-3->1e-4 | 1e-4->1e-5 |
|---|---|---|---|---|---|
| PH | bin10 | -17.45 +- 6.14 (-2.8) | -17.65 +- 6.78 (-2.6) | -21.02 +- 7.07 (-3.0) | -26.25 +- 24.13 (-1.1) |
| | bin11 | -8.47 +- 6.14 (-1.4) | -7.98 +- 6.78 (-1.2) | -10.98 +- 7.07 (-1.6) | -15.59 +- 24.13 (-0.6) |
| Z | bin10 | -13.26 +- 3.78 (-3.5) | -15.29 +- 5.15 (-3.0) | -11.31 +- 6.08 (-1.9) | -26.48 +- 22.90 (-1.2) |
| | bin11 | -4.32 +- 3.78 (-1.1) | -5.69 +- 5.15 (-1.1) | -1.27 +- 6.08 (-0.2) | -15.83 +- 22.90 (-0.7) |
| V9 | bin10 | -12.60 +- 4.08 (-3.1) | -16.00 +- 5.16 (-3.1) | -11.97 +- 6.30 (-1.9) | -2.79 +- 18.13 (-0.2) |
| | bin11 | -3.63 +- 4.08 (-0.9) | -6.36 +- 5.16 (-1.2) | -1.90 +- 6.30 (-0.3) | +7.89 +- 18.13 (+0.4) |

J7v's VEGAS-consistency results (parts PHv vs PH cell by cell) concern vi, kp, r with the same binary on both sides; `results/kp11/J7v_consistency.txt` is the bin11 version.

## 2. Full tables

`results/kp11/`: `J8_PH9_{main,parts,slopes}.txt`, `CH8_{main,slopes,pulls}.txt`, `J7_{PH,W}_{main,parts,slopes}.txt`, `J7v_{main,main_zeus,parts,slopes,consistency,scatter}.txt` and the JSON files, next to the bin10 versions of the same names in `results/`. dkp per cell and tau_cut with its slopes: `results/kp11/KPDELTA_{J8,CH8q,CH8g,J7PH,J7W,J7v9,J7v7}.txt`.
