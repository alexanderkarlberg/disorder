# nnlo21 cluster: KREDO (c)+(d) -- the bin5 validations, B, G6, Ex6, Ex7, X8 with the bin11 kp / kp-x -- report (16 Oct 2026)

Branch `2026-10-nnlo21-cluster`. UPDATE 16 Oct (b), second report (first: `REPORT-2026-10-16-KREDO.md`, (a)+(b), with bin11, the reruns and the method). Facts and significances only; no interpretation.
Each set is recombined with the script of its last report, run once without and once with `KP11=1`. Without it every committed result is reproduced (text and JSON identical; for B from the archived run logs, see below). With it, outputs go to `results/kp11/` under the same names. Same-seed dkp: `results/kp11/KPDELTA_{Z5,W5,E5,E5x,B,G6P,G6W,Ex6x,Ex7x,X8x}.txt` (`x` = kp-x).
Script changes for this round: `combine5.py`, `combine6.py`, `combine7e.py`, `x8.py` and `combine7n.py` read their inputs (NNLOJET references, bin5 JSON) from `results/` when `KP11=1` is set. `combine_b4.sh` takes `BPROD` (the parent of the archived `prod/B`, `prod/Bx`) and `KP11`. `kpdelta.py` also reads `CELL` rows.

## 0. Summary

- **Reruns**: Z, W, E kp 200 each; E kp-x 50; B kp 200 and Bx plain-kp 16; G6P, G6W kp 200 each; Ex6, Ex7, X8 kp-x 200 each; the 4-seed pilots. All done, every old done seed has a bin11 twin.
  The pilots and Bx plain-kp are not in any combination; they are rerun only for completeness.
- **B**: the B runs were archived in the 11 Oct file-quota recovery (`/ptmp/mpp/akarlber/archive/disorder-nnlo21_runs_prod_{B,Bx}.tar`). Their run.log/done files (b1, lo, b2, vi, kp, r-4, rcorr-4; Bx conv-4, edge-r-4, plain-r-4, plain-rcorr-4) were extracted to local scratch, not /ptmp (the file count is at 981k of the 1.05M soft limit). `combine_b4.sh plain` on them reproduces `results/B_4_plain.{txt,json}` exactly.
- **dkp checks**, total row / bin, tau_cut 1e-3 / 1e-4 / 1e-5, same seeds:

  | set | ours | MPP |
  |---|---|---|
  | W5 (ZEUS total, pb) | -0.1034 / -0.1676 / -0.2024 (+- 0.0003) | -0.1043 / -0.1683 / -0.2029 |
  | Z5 (ZEUS total, pb) | -2.996 / -4.962 / -6.963 | |
  | B (ZEUS total, pb) | -2.959 / -4.870 / -6.836 | |
  | E (tau_zQ 0.05-0.5) | -7.895 / -15.41 / -23.72 | |
  | G6P / G6W (0.05-0.5, pb/GeV^2) | -1.857e-4 / -3.209e-4 / -4.170e-4; -1.205e-4 / -2.067e-4 / -2.665e-4 | |

  The kp-x changes (E's P2B extras; total bin 0.05-0.5) are small next to their scatter. At 1e-3 / 1e-4 / 1e-5:
  - Ex7: -0.20 +- 0.13 / -0.03 +- 0.08 / +0.10 +- 0.05;
  - X8: +0.07 +- 0.12 / +0.01 +- 0.08 / +0.04 +- 0.05;
  - bin5: -1.6 +- 1.2 / +0.09 +- 0.14 / +0.004 +- 0.004;
  - Ex6: +3.7 +- 2.9 / -0.28 +- 0.58 / -0.013 +- 0.016 (Ex6 kp-x has per-seed sd ~90 at 2e-2).
- **The corrected results against NNLOJET** (total row; NNLO in pb; pulls at tau_cut 2e-2 1e-2 5e-3 2e-3 1e-3 5e-4 2e-4 1e-4 3e-5 1e-5):

  | set | kp | NNLO at 1e-3 / 1e-4 / 1e-5 | pulls |
  |---|---|---|---|
  | W vs N7 (0.1480 +- 0.0151) | bin10 | 0.1586 / 0.2724 / 0.3572 | -2.8 -4.5 -4.4 -1.0 +0.6 +1.7 +4.1 +5.3 +0.5 +4.4 |
  | | **bin11** | **0.0552 / 0.1048 / 0.1548** | **-4.3 -6.9 -7.9 -5.9 -5.3 -5.0 -3.1 -1.9 -1.6 +0.1** |
  | Z vs N7 (12.92 +- 1.03) | bin10 | 15.00 / 16.08 / 20.41 | -0.7 -0.6 -1.0 -0.2 +1.5 +1.8 +2.1 +1.7 +2.7 +2.0 |
  | | **bin11** | **12.00 / 11.12 / 13.45** | **-1.5 -1.8 -2.5 -2.1 -0.6 -0.6 -0.5 -1.0 +0.6 +0.1** |
  | Z vs bin5 NNLOJET (10.87 +- 1.2) | bin10 | | +1.1 +1.0 +0.7 +1.3 +2.7 +3.0 +3.1 +2.6 +3.4 +2.5 |
  | | bin11 | | +0.3 -0.0 -0.7 -0.5 +0.7 +0.7 +0.7 +0.1 +1.3 +0.7 |
  | W vs bin5 NNLOJET (0.1757 +- 0.056) | bin10 | | -1.3 -1.8 -1.7 -0.8 -0.3 +0.1 +1.0 +1.6 +0.1 +2.5 |
  | | bin11 | | -1.7 -2.4 -2.7 -2.3 -2.1 -2.1 -1.6 -1.2 -1.6 -0.3 |
  | B vs NNLOJET (11.12 +- 1.2) | bin10 | 12.36 / 14.77 / 7.74 | +0.8 +0.6 -0.1 +1.3 +0.8 +1.7 +1.6 +1.5 +1.1 -0.6 |
  | | **bin11** | **9.40 / 9.90 / 0.90** | **+0.1 -0.4 -1.4 -0.4 -1.0 -0.3 -0.5 -0.5 -0.5 -1.7** |
  | E plain, tau_zQ 0.05-0.5 (59.49 +- 12) | bin10 | 159.6 / 107.96 / 73.62 | -13.4 -4.1 +6.0 +8.2 +6.4 +3.8 +2.8 +2.5 +1.5 +0.4 |
  | | bin11 | 151.7 / 92.56 / 49.90 | -13.5 -4.3 +5.7 +7.7 +5.9 +3.2 +2.1 +1.7 +0.7 -0.3 |
  | E P2B (bin5 extras) | bin10 | 86.98 / 83.26 / 90.77 | -3.9 -1.9 +1.2 +1.8 +0.8 +1.0 +0.6 +0.6 +1.4 +0.8 |
  | | bin11 | 77.50 / 67.94 / 67.05 | -3.9 -2.0 +1.1 +1.7 +0.5 +0.7 +0.3 +0.2 +0.8 +0.2 |
  | E P2B [Ex7] | bin10 | 1.09 / 65.29 / 67.33 | -7.2 -7.2 -5.5 -4.8 -3.0 -3.2 -0.9 +0.3 +0.7 +0.2 |
  | | bin11 | -7.00 / 49.85 / 43.71 | -7.3 -7.3 -5.7 -5.1 -3.4 -3.6 -1.5 -0.5 -0.1 -0.4 |

  The W, Z and B rows above are the ZEUS total. Every bin, with the crossings and the |pull| < 2 ranges, is in `results/kp11/N7_{W,Z}_{main,summary}.txt`, `{Z,W,E}5_main.txt`, `B_4_plain.txt`, `Ex6_6_main.txt` and `Ex7_7_main.txt`. N7 per-tau_cut summary over the 15 W bins with bin11 (worst-bin pull / rms):
  - 2e-2 to 5e-4: worst -4.4 ... -7.9, rms 2.8-4.4;
  - 2e-4: worst -3.3, rms 1.9;
  - 1e-4: worst -2.6, rms 1.5;
  - 3e-5: worst -2.3, rms 1.3;
  - 1e-5: worst +2.3, rms 0.86.

  The same for Z: worst -2.0 ... -3.2 at tau_cut >= 2e-3, -2.5 ... +2.3 below; rms 0.8-1.7.
- **G6 (fixed point x = 0.1, Q^2 = 5000, inclusive tau_zQ cells; % of the bin NLO (DISENT / x) per unit ln, (pull))**: with the bin11 kp the total-bin slopes change sign and are positive at 3-10 sigma above 1e-4:

  | set, bin | kp | 5e-3->5e-4 | 2e-3->2e-4 | 1e-3->1e-4 | 1e-4->1e-5 |
  |---|---|---|---|---|---|
  | G6P 0.05-0.5 | bin10 | -1.25 +- 0.25 (-5.1) | +0.75 +- 0.70 (+1.1) | +0.35 +- 0.40 (+0.9) | -2.22 +- 1.25 (-1.8) |
  | | **bin11** | **+0.77 +- 0.25 (+3.1)** | **+2.80 +- 0.70 (+4.0)** | **+2.36 +- 0.40 (+5.8)** | **-0.78 +- 1.25 (-0.6)** |
  | G6W 0.05-0.5 | bin10 | -1.35 +- 0.26 (-5.3) | +1.32 +- 0.33 (+4.0) | -0.35 +- 0.44 (-0.8) | -1.18 +- 0.88 (-1.3) |
  | | **bin11** | **+0.78 +- 0.26 (+3.0)** | **+3.48 +- 0.33 (+10.5)** | **+1.76 +- 0.44 (+4.0)** | **+0.29 +- 0.88 (+0.3)** |

  By bin (bin11, 1e-3 -> 1e-4): G6P 0.05-0.1 +3.95 +- 0.76 (+5.2), 0.1-0.2 +2.30 +- 0.62 (+3.7), 0.2-0.3 -0.91 +- 0.80, 0.3-0.4 -1.34 +- 1.11, 0.4-0.5 +3.56 +- 1.27 (+2.8). G6W 0.05-0.1 +3.59 +- 0.85 (+4.2), the other bins within 1.6 sigma. The largest single pulls (2e-3 -> 2e-4, the 0.05-0.1 bin: G6P +9.9, G6W +10.5) were already +6.8 / +7.5 with bin10. The dkp slope of the total bin is +5.9e-5 (G6P) / +3.7e-5 (G6W) pb/GeV^2 per unit ln, about 2.0% of the bin NLO.
- **Ex6, Ex7, X8**: only kp-x changes, and its change is within about 2 sigma of zero (above), so the P2B extras tables are the same within errors. The E P2B shifts in the table above come from the plain kp (bin5 E kp), not from kp-x. `results/kp11/X8_total.txt`, `X8_tables.md`, `Ex7_7_{main,slopes,window,consistency}.txt`, `Ex6_6_*`.

## 1. Files

`results/kp11/`: `{Z5,W5,E5}_{main,parts,slopes}.txt` + JSON; `N7_{W,Z}_{main,summary}.txt` + JSON; `B_4_plain.{txt,json}`; `{G6P,G6W,Ex6}_6_*`; `Ex7_7_*`; `X8_total.{txt,json}`, `X8_tables.md`; `KPDELTA_*.txt`, `KPDELTA.json`. The bin10 versions of the same names are in `results/`.
