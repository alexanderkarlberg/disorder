"""Generate FORM programs for squared tree amplitudes of
lepton + partons with one quark line (or two) coupled to the boson,
resolved into lepton helicity and quark-line chirality.

A diagram on a fermion line is a list of factors read from the
outgoing-fermion end (ubar or vbar side) to the incoming end:
  ("g", idx)             gamma^idx
  ("S", [(+1,"p2"),...], "D23")  (sum of signed momenta)slash * D23
plus a scalar prefactor string that may carry Lorentz indices
(non-abelian vertex). Conjugate amplitudes are built by reversing the
factor list and renaming indices with CONJ."""

CONJ = {"mu": "nu", "a3": "b3", "a4": "b4", "s1": "s2", "a2": "b2"}

def conj_idx(s):
    for a, b in CONJ.items():
        s = s.replace(f"({a})", f"({b})").replace(f"({a},", f"({b},").replace(f",{a})", f",{b})")
    return s

def slash(mom, line):
    return "(" + "+".join(f"{'-' if s < 0 else ''}g_({line},{p})" for s, p in mom).replace("+-", "-") + ")"

def string(factors, line, conj=False):
    out = []
    fs = list(reversed(factors)) if conj else factors
    for f in fs:
        if f[0] == "g":
            i = CONJ.get(f[1], f[1]) if conj else f[1]
            out.append(f"g_({line},{i})")
        else:
            out.append(slash(f[1], line) + "*" + f[2])
    return "*".join(out)

def amp(diagrams, line, conj=False):
    terms = []
    for coef, factors in diagrams:
        c = conj_idx(coef) if conj else coef
        terms.append(f"({c})*{string(factors, line, conj)}")
    return "\n + ".join(terms)
