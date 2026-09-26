"""Check gluon.frm against the DISENT GQ form: X(l,h) = c * PAIR(l,h) *
s12*s13*(-q^2), PAIR(same) = (k.p3)^2 + (kp.p2)^2 (q at p2, qbar at p3),
PAIR(opp) = (kp.p3)^2 + (k.p2)^2."""
import re, subprocess, sympy as sp
out = subprocess.run(["form", "-q", "gluon.frm"], capture_output=True, text=True).stdout
kkp, kp1, kp2, kpp1, kpp2 = sp.symbols("kkp kp1 kp2 kpp1 kpp2")
names = {"k.kp": "kkp", "k.p1": "kp1", "k.p2": "kp2", "kp.p1": "kpp1", "kp.p2": "kpp2"}
ex = {}
for m in re.finditer(r"(X\d\d)=\s*(.*?);", out, re.S):
    s = re.sub(r"\s+", "", m.group(2))
    for a, b in sorted(names.items(), key=lambda t: -len(t[0])):
        s = s.replace(a, b)
    ex[m.group(1)] = sp.sympify(s.replace("^", "**"))
p1p2 = kp1 - kkp - kp2 - kpp1 + kpp2
p1p3 = kp1 - kpp1 - p1p2
kp3 = kp1 - kp2 - kkp          # k.(k+p1-kp-p2)
kpp3 = kkp + kpp1 - kpp2       # kp.(k+p1-kp-p2)
s12, s13, mq2 = 2 * p1p2, 2 * p1p3, 2 * kkp
same = kp3**2 + kpp2**2
opp = kpp3**2 + kp2**2
for key, pair in (("X66", same), ("X77", same), ("X67", opp), ("X76", opp)):
    print(key, sp.factor(sp.cancel(ex[key] / (pair * s12 * s13 * mq2))))
