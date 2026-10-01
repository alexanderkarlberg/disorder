# N-jettiness slicing for DIS (branch 2026-10-tau-slicing)

Goal: an NNLO DIS 2+1 calculation for fully differential N3LO via P2B
(disorder, and proVBFH line by line). First step: the SCET ingredients at
NLO, against exact references.

Build: `cmake -S . -B build -DDISORDER_SLICING=ON` (off by default), then
`cmake --build build --target test_scet tau2_nlo`.

## 1+1 at NLO: `tau1b_nlo.py`

tau_1^b of Kang, Lee, Stewart (arXiv:1303.6952) at fixed (x, Q^2, y),
KLS cumulant (173)-(174) below the cut, O(alpha_s) real emission above,
against the exact NLO (MS-bar coefficient functions, checked against
hoppet). Converges as tau_cut ln tau_cut (notebook, 2026-10-01).
`--cumulant-only` prints the tau_1^a and tau_1^b cumulants (used as the
reference of `test_scet`).

## 2+1 at NLO: `mod_slicing_scet.f90`, `tau2_nlo.f90`

Conventions of Gaunt, Stahlhofen, Tackmann, Walsh (GSTW, arXiv:1505.04794):
geometric measure with Q_i = 2 E_i in the Breit frame (proton along +z),

    T_2 = sum_k min_i n_i.p_k ,   n_i = (1, nhat_i),
    T_2(3 final partons) = min( min_j (E_j - p_zj), min_{j<k} (E_j + E_k - |p_j + p_k|) ),

tau_2 = T_2/Q. With xi = T_cut the slicing cumulant is GSTW's C_{-1}^(1)
(3.32), here in units of alpha_s/(2 pi) (GSTW use alpha_s/(4 pi)) at
mu = mu_F = Q:

- jets (GSTW (A.9)-(A.12)): lambda_J = 2 E_J T_cut / mu^2;
  quark (1/2)[CF(7 - pi^2) - 3 CF L + 2 CF L^2],
  gluon (1/2)[CA(4/3 - pi^2) + 5 beta0/3 - beta0 L + 2 CA L^2].
- beam (GSTW (A.16), (A.19)): lambda_B = 2 E_a T_cut / mu^2, E_a the
  energy of the incoming Born parton;
  B = sum_j int dz [ I^(1)_aj(z) + (1/2) I_aj,0(z) L + (1/2) I_aj,1(z) L^2/2 ] (x f_j)(eta/z)
  with I^(1)_qq, I^(1)_qg (arXiv:1401.5478, = KLS (171)) and I^(1)_gg,
  I^(1)_gq (arXiv:1405.1044 (A.4)); (1/2) I_qq,0 = CF (1+z^2)[1/(1-z)]_+,
  (1/2) I_gg,0 = CA 2(1-z+z^2)^2/z [1/(1-z)]_+, (1/2) I_ii,1 L^2/2 = C_i L^2 delta(1-z).
- soft (GSTW (A.24), (A.31), (A.32)): lambda_S = T_cut/mu, s_ij = n_i.n_j/2;
  (1/2)[ sum_{i/=j} T_i.T_j (ln^2 s_ij - zeta2 + 4 sum_m I_ij,m) - 4 L sum_{i/=j} T_i.T_j ln s_ij - 4 L^2 sum_i C_i ],
  I_ij,m = I0(a,b) ln a + I1(a,b), a = s_jm/s_ij, b = s_im/s_ij, with the
  finite integrals of Jouttenus, Stewart, Tackmann, Waalewijn
  (arXiv:1102.4344, (56)-(60)), here as one-dimensional phi integrals
  with the y integral done analytically (complex dilogarithm). The
  region in (60) is y > sqrt(b/a) (from x = p_j/p_i > 1); the PDF text of
  the paper shows it garbled.
- hard: built from DISENT. VIRTHR's constant QQ (GG) is the finite part
  of [V + I_CS] (factorising virtual plus the Catani-Seymour I operator,
  CS (7.31): prefactor (4 pi mu^2)^eps / Gamma(1-eps), which equals the
  MS-bar one-loop prefactor up to O(eps^3) after mu^2 -> mu^2 e^gamma/(4 pi)).
  Writing V = Born (alpha_s/2pi) [mu^2/Q^2]^eps (1 - zeta2 eps^2/2)(A/eps^2 + B/eps + V0)
  with A = -sum_i C_i and minimal subtraction of the poles:
      H^(1)/Born = QQ + I0_CS + (pi^2/12) sum_i C_i        (mu = Q)
      I0_CS = sum_i [C_i pi^2/3 - gamma_i - K_i] + sum_{i/=k} T_i.T_k [l_ik^2/2 - (gamma_i/C_i) l_ik],
  l_ik = ln(2 p_i.p_k/Q^2), plus DISENT's non-factorising virtual (LEIV,
  ERTV), which is finite and not proportional to the Born. Checks: for the
  two-parton Born VIRTWO's QQ = CF(2 - pi^2) gives CF(-8 + pi^2/6), the
  DIS quark form factor; for three partons all single logs l_ik cancel,
  leaving const + sum T_i.T_k l_ik^2 as expected at mu = Q.

`test_scet` checks the complex dilogarithm (mpmath), I0 and I1 (mpmath,
30 digits; a plain scipy 2D quadrature was off by 3e-5), the two-parton
hard function, and the complete 1+1 cumulant (beam, jet, two-direction
soft, hard) against `tau1b_nlo.py --cumulant-only` (5e-7).

`tau2_nlo`: one DISENT run (DISENTFULL, ORDER 2, photon exchange, fixed x
and Q^2) gives in the same events
- the reference: all O(alpha_s^2) contributions;
- above the cut: the real four-parton events with tau_2 > tau_cut
  (counter-events and collinear terms have tau_2 = 0 exactly);
- below the cut: the three-parton Born reweighted with the cumulant.
Observable: tau_zQ bins above 0.05. `tau2_combine.py` combines runs.

## 1+1 at NNLO three ways: `nnlo11.f90`, `mcfm_lp/`, `nnlo11_ref.py`, `nnlo11_combine.py`

One DISENT run (fixed x, Q^2) gives, for lab-frame jet observables of 1+1
(anti-k_t R = 1, p_t > 5 GeV, -1 < eta < 2.5), the O(alpha_s^2) coefficient
- by DISENT + P2B (DISENT's O(alpha_s^2) 2+1 weights minus their Born
  projection, plus the inclusive part C2_incl O(Born) from hoppet),
- by tau_2 slicing + P2B (the same with the 2+1 part sliced with the invariant
  measure, absolute cuts or cuts relative to tau_1 of the Born),
- by pure tau_1 slicing: tau_1 = min over partitions [2 x P.p(beam) + m^2(jet)]/Q^2
  (recoil-free), DISENT's O(alpha_s^2) weights above the cut and the NNLO
  leading-power cumulant below it.
`nnlo11` needs DISENT's two-parton Born given to USER (`USER(2,0,0)`,
commented out in `src/libdisent.f`): build a copy of the fixed DISENT with
that call and link it instead of `libdisent.f.o`.
`mcfm_lp/build.sh /path/to/MCFM-10.3` builds `dis_tau1_lp`, the NNLO
leading-power cumulant (per Born, in (alpha_s/2pi)^k at mu = Q) from MCFM's
0-jettiness pieces: beam = quark beam function, the second "beam" = quark jet
function (MCFM's `jetq` is in alpha_s/4pi: J1/2, J2/4), soft = qqbar
0-jettiness soft function (equal for DIS at O(alpha_s^2), arXiv:1501.04110),
hard = `hardqq` at Qsq = -Q^2 (spacelike). Its NLO part equals the tau_1^a
cumulant of `tau1b_nlo.py` to 2e-6. `nnlo11_ref.py`: the exact inclusive
coefficients from hoppet. `nnlo11_combine.py` combines runs.
