# MCFM 10.3 routines for DIS 3+1 / 4+1 (GPL-3.0-or-later)

Copied from MCFM 10.3 (`src/Z2jet`, `src/W2jet`, `src/Need`, `src/Inc`),
crossed to DIS by the callers (incoming parton and incoming lepton as
negative momenta in MCFM's all-incoming convention). Changes are marked
with "disorder (date)":
- `spinoru.f`: Minkowski product inline (MCFM's `dot` clashes with DISENT's
  `DOT`).

Checks (scratch harness, 2 Oct 2026): `z2jetsq` crossed to gamma* g -> q
qbar g agrees with DISENT's MATFOR (photon exchange, incoming gluon,
summed over the labellings of the outgoing partons) up to a constant
normalisation, 3e-14 over 20 random points.
- `ampqqb_qqb.f`, `aqqb_zbb.f` (`src/Zbb`): four-quark amplitudes. They
  take the spinor products from the common block `/zprods/` and have the
  lepton pair fixed in slots 3, 4.

Crossing rules found with the harnesses (`dis31/tests`, photon exchange,
against DISENT's MATFOR, all to 1e-13 up to one constant):
- slots: 1 = -incoming parton, 2 = outgoing quark of the incoming line,
  3 = outgoing lepton, 4 = -incoming lepton, 5, 6 = the other two partons;
- q g g: `z2jetsq(2, 1, 3, 4, 5, 6)` (the outgoing quark in MCFM's quark
  slot, the incoming quark in its antiquark slot);
- q Q Qbar: `ampqqb_qqb(2, 1, 5, 6)` (Q at 5, Qbar at 6); the boson on the
  incoming line (A) and on the pair (B), amplitude e_q A + e_Q B;
- identical quarks: the squares from `ampqqb_qqb(2,1,5,6)` and the
  exchanged `ampqqb_qqb(5,1,2,6)`; the interference only from MCFM's own
  construction, crossed (MCFM slots 1,2,5,6 = ours 1,6,2,5): direct
  `ampqqb_qqb(1,2,6,5)` with j2 swapped and B negated, exchange
  `ampqqb_qqb(1,5,2,6)`, interference (2/N) Re[(Ad - Bd)(j,swap(j)) (Ae + Be)*(j,j)]
  (other argument orders give the right squares but inconsistent phases);
- weights relative to each other: q g g (N/4)(1/2) sum(msq), four-quark
  4 |.|^2, identical (1/2) 4 [squares + interference].
