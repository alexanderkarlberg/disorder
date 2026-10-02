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
