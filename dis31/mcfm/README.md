# MCFM 10.3 routines for DIS 3+1 / 4+1 (GPL-3.0-or-later)

Copied from MCFM 10.3 (`src/Z2jet`, `src/Zbb`, `src/Wbb`, `src/Need`,
`src/Inc`), crossed to DIS by the callers (`../me31.f90`, `../me41.f90`;
incoming parton and incoming lepton as negative momenta in MCFM's
all-incoming convention). Changes are marked with "disorder (date)":
- `spinoru.f`: Minkowski product inline (MCFM's `dot` clashes with DISENT's
  `DOT`).
- `subqcdn.f` (`src/W2jet`), `spinork.f`, `checkndotp.f` (`src/Need`),
  `Inc/nwz.f`: q qbar g g with one gluon contracted with a vector n (spin
  correlations; `nwz` must be 0 for photon/Z).
- `msq_gqqQQg.f`, `makemb_photon.f`: photon-exchange versions of
  `msq_ZqqQQg` and `makemb` (with `makem` inlined); the charges of the two
  quark lines are arguments instead of the electroweak common blocks.

- `loop/`: the one-loop amplitudes of Bern, Dixon, Kosower for q qbar g g + V
  (hep-ph/9708239) and q qbar Q Qbar + V (hep-ph/9610370) as in MCFM 10.3
  (`src/BDK`, `src/W2jet`, `src/Wbb`, `src/Z2jet`, `src/Need/lnrat.f`,
  `lib/SpecialFunctions/ddilog.f`, `dclaus.f`; 51 files, the dependency
  closure of the routines used, found from MCFM's object files), plus
  `loop/xzqqgg_wrappers.f` (helicity wrappers from `src/Zbb/xzqqgg*.f`).
  8 Oct 2026: `loop/fvs.f` (`src/BDK`), `loop/fvf.f` (`src/W2jet`) and
  `a64v` (in the wrappers, from `src/Zbb/xzqqgg_v.f`), unchanged: the
  boson on a closed quark loop with vector coupling in q qbar g g (BDK's
  A6^v); checked against MCFM's own `xzqqgg_v` (`mqqb_vec0`) through
  virt31 at DIS points (same normalisation as its main term, photon, e+,
  nu, all three channels).
  Changed: `a6routine.f` and `a61g.f` stop instead of calling the exact
  top-loop routines (not ported; we use toploops = none, n_f = 5).
  Includes added: `epinv.f`, `epinv2.f` (MCFM's poles: epinv = epinv2 =
  1/eps, double pole epinv*epinv2), `heldefs.f`, `scale.f`, `toploops.f`,
  `masses.f`.

Routines:
- 3+1: `z2jetsq.f` (+ `storecsz.f`, `subqcd.f`): q qbar g g; `ampqqb_qqb.f`,
  `aqqb_zbb.f` (`src/Zbb`): q qbar Q Qbar. The four-quark routine takes the
  spinor products from `/zprods/` and has the lepton pair fixed in slots 3, 4.
- 4+1: `xzqqggg.f` (+ `amp_qqggg.f`): q qbar g g g (couplings g_s, e from
  `/qcdcouple/`, `/ewcouple/`, all colour structures with `/ColC/` = 0);
  `msq_gqqQQg.f` (+ `makemb_photon.f`, `nagyqqqqg.f` from `src/Wbb`):
  q qbar Q Qbar g, non-identical (MN) and identical (MI) quarks.

## Crossing rules (photon exchange)

Slots: 1 = -incoming parton, 2 = outgoing quark of the incoming line,
3 = outgoing lepton, 4 = -incoming lepton, then the other partons.

3+1 (`me31`):
- q g g: `z2jetsq(2, 1, 3, 4, 5, 6)` (the outgoing quark in MCFM's quark
  slot, the incoming quark in its antiquark slot); g -> q qbar g:
  `z2jetsq(2, 5, 3, 4, 1, 6)` (q at 2, qbar at 5).
- q Q Qbar: `ampqqb_qqb(2, 1, 5, 6)` (Q at 5, Qbar at 6). The incoming line
  is read as (2,1), opposite to MCFM's orientation (`qqb_z2jet` uses (1,2)
  for the same channel), so the amplitude is e_q A - e_Q B. (Reversing one
  quark line flips the sign of the charge-odd e_q e_Q term; checked:
  `ampqqb_qqb(1,2,5,6)` with + equals `(2,1,5,6)` with -.)
- identical quarks: direct D = A - B from `ampqqb_qqb(2,1,5,6)`, exchange
  E = Ae - Be from `ampqqb_qqb(5,1,2,6)`; |M|^2 = 4V e_q^2 [sum |D|^2 +
  sum |E|^2 + (2/N) sum_{j,j3} Re D(j,swap(j),j3) E*(j,swap(j),j3)].

4+1 (`me41`), MCFM's own sign conventions (as `qqb_z2jet_g`):
- q g g g: `xzqqggg(2, 5, 6, 7, 1, 3, 4)`; g -> q qbar g g:
  `xzqqggg(2, 1, 6, 7, 5, 3, 4)` (q at 2, qbar at 5).
- q Q Qbar g: `msq_gqqQQg(2, 1, 5, 6, 7, 4, 3, e_q, e_Q)` (Q at 5, Qbar at
  6, g at 7), MN; identical quarks the same with e_Q = e_q, MI.
- g -> q qbar Q Qbar: `msq_gqqQQg(2, 5, 6, 7, 1, 4, 3, e_q, e_Q)` (q, qbar,
  Q, Qbar at 2, 5, 6, 7).

Colour- and spin-correlated 3+1 Borns (`born31`): the colour-ordered
amplitudes of `subqcd` (A1: gluon order (A, B) of the call, A2: (B, A))
belong to c1 = (T^A T^B)_{q qbar}, c2 = (T^B T^A), with q = MCFM's first
quark slot (our outgoing quark) and qbar = its second (our incoming quark,
or the outgoing antiquark of g -> q qbar g). `subqcdn`'s qcdab/qcdba have
the opposite assignment (qcdab <-> the order (contracted, other)), and its
n-contraction is normalised to half the polarisation sum. For an incoming
antiquark the q qbar g g roles stay as for a quark (charge conjugation
swaps q <-> qbar and the colour order). Four quarks: the direct colour
structure D = T^a_{q1 qb2} T^a_{q3 qb4} with q1 = the outgoing quark of the
incoming line, qb2 = the incoming quark, q3, qb4 = the pair (for an
incoming antiquark q1 = the incoming antiquark, qb2 = the outgoing one, q3 =
the pair quark), and the identical-quark exchange E = -(MCFM's exchange
amplitude) in this basis.

One loop (`virt31`): q qbar g g from `a6treeg1`, `a61g1lc`, `a61g1slc`,
`a61g1nf`, `a63g1` with xzqqgg_v's colour weights (colourchoice 0), MCFM
labels 1 = outgoing quark, 2, 3 = gluons, 4 = -incoming quark (or the
outgoing antiquark), 5 = outgoing lepton, 6 = -incoming lepton; four quarks
from `atreez`, `a61z`, `a62z` as qqb_z2jet_v's q q channel, MCFM slots
1 = -incoming quark, 5 = outgoing quark of that line, 2 = the pair's
antiquark, 6 = its quark (identical quarks: 5 <-> 6). The vector-loop
(a63z, a64v) and axial pieces vanish for the photon and are left out.

## Checks (`../tests`, 2 Oct 2026)

- me31 against DISENT's MATFOR (`harness_me31`), summed over the
  labellings of the outgoing partons: ratio 1 to 1e-13 for all incoming
  flavours. This sum is blind to the charge-odd e_q e_Q terms of the four-
  quark channels (they are odd under Q <-> Qbar and cancel in it).
- me31 four-quark channels pointwise against Feynman diagrams with explicit
  Dirac matrices (`fd31.py`): ratio 1 to 1e-12 (d -> d u ubar, u -> u u
  ubar, dbar -> dbar ubar u). This fixed the sign of the e_q e_Q term and
  the identical-quark interference.
- born31 (`harness_born31`): msq = me31 (4e-16); colour conservation
  (1e-15); polarisation sums of the spin-correlated Borns (3e-14); soft-
  gluon limit of me41 against the eikonal sum with born31's T_i.T_k, nine
  channels including antiquarks, gluon-initiated and identical quarks
  (2e-4 at lambda = 1e-5); gluon-splitting collinear limits at fixed
  azimuth with the spin-correlated Borns (final g -> g g, g -> q qbar;
  initial q -> g, g -> g; 5e-4 at y = 1e-9). Colour matrices from
  `colour.py` (explicit SU(3)).
- me41 against me31 in single-collinear limits (`harness_lim41`): final-state
  q||g, g||g, g -> q qbar, initial-state q -> q g, g -> q qbar, all channels;
  ratio 1 to < 5e-4 at y = 1e-10, approaching 1 linearly in y.
- virt31 (`harness_virt31`), nine channels, two scales: tree = me31
  (2e-14); double pole -sum_i C_i |M0|^2 (3e-14); single pole of the
  renormalised virtual = the I-operator prediction with born31's T_i.T_k
  (7e-13). Finite part against NNLOJET v1.0.2 (outside the repository,
  `~/cernbox/disorder-comparisons/dis31_nnlojet`): equal up to
  (pi^2/12) sum_i C_i |M0|^2 (NNLOJET's normalisation) in q -> q g g,
  g -> q qbar g, four quarks (all charge structures) and identical quarks,
  at two scales.
