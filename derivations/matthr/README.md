# DISENT's three-parton matrix element for gamma/Z and W exchange

FORM derivations behind the generalisation of `MATTHR` (`src/libdisent.f`)
from photon exchange to NC (γ, Z, γ/Z) and CC exchange, for charged
leptons and neutrinos of either charge.

For each process the tree-level |M|² of lepton + parton → lepton + 2
partons is computed with chiral projectors on the lepton line (helicity l)
and on the quark line (field chirality h). Couplings, colour factors and
1/q⁴ are stripped off, the gluon polarisation sum is −g, and momenta are
labelled as in DISENT: 1 incoming parton, 2 and 3 outgoing partons,
k = P(,6) incoming lepton, k' = P(,7) outgoing lepton, q = k − k'.

| file | process | result (A(l,h) = L_{μν}H^{μν}) |
|---|---|---|
| `quark.frm` | l q(p1) → l q(p2) g(p3) | 32 (−q²) PAIR / (s13 s23) |
| `antiquark.frm` | l q̄(p1) → l q̄(p2) g(p3) | the same with the two pairs swapped |
| `gluon.frm` | l g(p1) → l q(p2) q̄(p3) | 32 (−q²) PAIR / (s12 s13) |

with s_ij = 2 p_i·p_j and

* quark: PAIR = (k·p1)² + (k'·p2)² if l = h (LL, RR),
  PAIR = (k'·p1)² + (k·p2)² if l ≠ h (LR, RL);
* antiquark: the other way round;
* gluon (quark at p2, antiquark at p3): PAIR = (k·p3)² + (k'·p2)² if
  l = h, (k'·p3)² + (k·p2)² otherwise.

`check_quark.py`, `check_antiquark.py` and `check_gluon.py` run FORM,
eliminate p3 and p1·p2 through momentum conservation and p3² = 0, and
verify these identities exactly with sympy (each prints the ratio, 32).

## Consequences for MATTHR

The sum over the four helicity combinations with photon couplings gives
back DISENT's `QQ` (quark) and `GQ` (gluon); including the spin and colour
averages and g_s² = 4π α_s, `QQ` equals |M|²/(α_s/2π) exactly. For
general couplings write, for incoming parton i,

    M(i) = C2(i) QQ + C3(i) QQ3,
    QQ3  = QQ with (k·p1)² + (k'·p2)² − (k'·p1)² − (k·p2)² in the numerator,

where C2 = (S + O)/2 and C3 = (S − O)/2, S (O) being the coupling of the
same (opposite) helicity lepton-parton configuration, in units of the
photon coupling e_q². At LO the same combination gives
Y₊ C2 + Y₋ C3 (Y± = 1 ± (1−y)²), so C2 and C3 are the per-parton
couplings of F2/x and F3 (`parton_couplings` in
`src/mod_matrix_element.f90`, written to mirror `eval_matrix_element_new`
and checked against the quark-parton model in `tests/test_matrix_element`).
Antiquarks have C3 → −C3.

In the gluon channel the parity-violating part is odd under the exchange
of the quark and the antiquark (p2 ↔ p3). DISENT's three-parton phase space
(`GENTHR`/`GENDEC`) is symmetric under this exchange, and no observable
distinguishes the flavour of partons 2 and 3, so it integrates to zero and
is dropped. The gluon coupling per final-state flavour pair is
(C2(i) + C2(−i))/2 (e_q² for photons).

For photon exchange C2 = EQ², C3 = 0 and the gluon coupling is EQ², all
exact in floating point, so MATTHR reproduces the original DISENT
result bit for bit (`tests/test_matthr`, and the p2b validation outputs).

Run: `python3 check_quark.py` etc. (needs `form` in PATH and sympy).
