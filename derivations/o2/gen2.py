"""All-outgoing amplitude builder for  V*(q) -> partons  with a lepton
current attached (see gen.py for the conventions of factor lists).
Momenta are generic labels (qa: outgoing quark, qb: outgoing antiquark,
pc, pd: gluons; qc, qd: second quark pair) that are mapped onto signed
DIS momenta per channel, e.g. an incoming quark p1 is qb = -p1.
Spin sums: u ubar = v vbar = pslash for the physical momentum p; in
terms of the all-outgoing label (-p1) this is -(label)slash, so every
crossed fermion gives a factor -1 (included in `sign`)."""

def lin(mom, sub):
    """signed momentum list in generic labels -> list of (sign, dis momentum)"""
    out = {}
    for s, lab in mom:
        for s2, p in sub[lab]:
            out[p] = out.get(p, 0) + s * s2
    return [(s, p) for p, s in out.items() if s != 0]

def slash(mom, line):
    t = []
    for s, p in mom:
        t.append(("+" if s > 0 else "-") + (f"{abs(s)}*" if abs(s) != 1 else "") + f"g_({line},{p})")
    return "(" + "".join(t).lstrip("+") + ")"

def vec(mom, idx):
    t = []
    for s, p in mom:
        t.append(("+" if s > 0 else "-") + (f"{abs(s)}*" if abs(s) != 1 else "") + f"{p}({idx})")
    return "(" + "".join(t).lstrip("+") + ")"

CONJ = {"mu": "nu", "ac": "bc", "ad": "bd", "s1": "s2", "t1": "t2"}

def cidx(i, conj):
    return CONJ.get(i, i) if conj else i

def string(factors, line, sub, conj=False):
    """factors: ("g", idx) or ("S", genericmom, densym) or ("P", h) projector"""
    fs = list(reversed(factors)) if conj else factors
    out = []
    for f in fs:
        if f[0] == "g":
            out.append(f"g_({line},{cidx(f[1], conj)})")
        elif f[0] == "S":
            out.append(slash(lin(f[1], sub), line) + "*" + f[2])
        elif f[0] == "P":
            out.append(f"(1/2)*g{f[1]}_({line})")
        elif f[0] == "p":   # external spinor sum, physical momentum
            out.append(slash(lin(f[1], sub), line))
    return "*".join(out)

def scalar(expr, sub, conj=False):
    """scalar prefactor: string with {vec:label:idx} placeholders"""
    import re
    def rep(m):
        lab, idx = m.group(1), m.group(2)
        return vec(lin([(1, lab)], sub), cidx(idx, conj))
    s = re.sub(r"\{(\w+):(\w+)\}", rep, expr)
    if conj:
        for a, b in CONJ.items():
            s = re.sub(rf"\b{a}\b", b, s)
    return s
