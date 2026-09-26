"""Numerical e q(p1) -> e q(p2) q(p3) qbar(p4) with photon couplings:
direct pairs (1->2),(3,4); exchanged (1->3),(2,4). Returns spin sums of
|A_dir(boson on line 1)|^2 (D1-type), |A_dir(boson on pair)|^2 (D2-type)
and 2 Re A_dir A_exch^* (E-type), all without couplings/colour/1/q^4."""
import numpy as np
from spinors import G, slash, spinors, bar, up
def mdot(a, b):
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3]
def msq(p): pu = up(p); return mdot(pu, pu)
def add(*ps):
    return [sum(p[i] for p in ps) for i in range(4)]
def neg(p): return [-x for x in p]
def S(p): return slash(p) / msq(p)
LOW = [1, -1, -1, -1]

def line_boson(ubar, v, o, i, g1g2, mu, sg):
    """ubar(o) [g^sg S(o+g) g^mu + g^mu S(-i-g) g^sg] v(i), g = g1g2 sum; i given as its label momentum (signed)"""
    return ubar @ (G[sg] @ S(add(o, g1g2)) @ G[mu] + G[mu] @ S(neg(add(i, g1g2))) @ G[sg]) @ v

def amps(P):
    k, kp, p1, p2, p3, p4 = P[6], P[7], P[1], P[2], P[3], P[4]
    q = [k[i] - kp[i] for i in range(4)]
    res = {"D1": 0.0, "D2": 0.0, "E": 0.0}
    Lc = [[bar(ukp) @ G[m] @ uk for m in range(4)] for uk in spinors(k) for ukp in spinors(kp)]
    m1 = neg(p1)   # all-outgoing label of the incoming quark
    for L in Lc:
        for u1 in spinors(p1):
            for u2 in spinors(p2):
                for u3 in spinors(p3):
                    for v4 in spinors(p4):
                        b2, b3 = bar(u2), bar(u3)
                        # contract boson index with lepton current (lower index)
                        def dirl1():  # boson on line 1->2, gluon -> (3,4)
                            return sum(LOW[m]*LOW[s] * L[m] * line_boson(b2, u1, p2, m1, add(p3, p4), m, s) * (b3 @ G[s] @ v4)
                                       for m in range(4) for s in range(4)) / msq(add(p3, p4))
                        def dirl2():  # boson on pair (3,4), gluon from line 1->2
                            return sum(LOW[m]*LOW[s] * L[m] * (b2 @ G[s] @ u1) * line_boson(b3, v4, p3, p4, add(p2, m1), m, s)
                                       for m in range(4) for s in range(4)) / msq(add(p2, m1))
                        def exl1():   # exchanged: boson on line 1->3, gluon -> (2,4)
                            return sum(LOW[m]*LOW[s] * L[m] * line_boson(b3, u1, p3, m1, add(p2, p4), m, s) * (b2 @ G[s] @ v4)
                                       for m in range(4) for s in range(4)) / msq(add(p2, p4))
                        def exl2():   # exchanged: boson on pair (2,4), gluon from line 1->3
                            return sum(LOW[m]*LOW[s] * L[m] * (b3 @ G[s] @ u1) * line_boson(b2, v4, p2, p4, add(p3, m1), m, s)
                                       for m in range(4) for s in range(4)) / msq(add(p3, m1))
                        a1, a2, e1, e2 = dirl1(), dirl2(), exl1(), exl2()
                        res["D1"] += abs(a1)**2; res["D2"] += abs(a2)**2
                        res["E"] += 2 * ((a1 + a2) * np.conj(e1 + e2)).real
    return res
