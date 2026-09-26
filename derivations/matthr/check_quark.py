"""Check the FORM result (quark.frm) against the DISENT MATTHR form:
X(l,h) = A(l,h) s13^2 s23^2 must equal c * PAIR(l,h) * s13*s23 * (-q^2)
with PAIR(same) = (k.p1)^2 + (kp.p2)^2, PAIR(opp) = (kp.p1)^2 + (k.p2)^2."""
import re, subprocess, sympy as sp
out = subprocess.run(["form", "-q", "quark.frm"], capture_output=True, text=True).stdout
kkp, kp1, kp2, kpp1, kpp2 = sp.symbols("kkp kp1 kp2 kpp1 kpp2")
names = {"k.kp": "kkp", "k.p1": "kp1", "k.p2": "kp2", "kp.p1": "kpp1", "kp.p2": "kpp2"}
exprs = {}
for m in re.finditer(r"(X\d\d)=\s*(.*?);", out, re.S):
    s = re.sub(r"\s+", "", m.group(2))
    for a, b in sorted(names.items(), key=lambda t: -len(t[0])):
        s = s.replace(a, b)
    exprs[m.group(1)] = sp.sympify(s.replace("^", "**"))
p1p2 = kp1 - kkp - kp2 - kpp1 + kpp2          # from p3^2 = 0
p1p3 = kp1 - kpp1 - p1p2                       # p1.(k+p1-kp-p2)
p2p3 = kp2 + p1p2 - kpp2                       # p2.(k+p1-kp-p2)
s13, s23, mq2 = 2 * p1p3, 2 * p2p3, 2 * kkp    # -q^2 = 2 k.kp
same = kp1**2 + kpp2**2
opp = kpp1**2 + kp2**2
for key, pair in (("X66", same), ("X77", same), ("X67", opp), ("X76", opp)):
    r = sp.factor(sp.cancel(exprs[key] / (pair * s13 * s23 * mq2)))
    print(key, "ratio to PAIR*s13*s23*(-q^2):", r)
