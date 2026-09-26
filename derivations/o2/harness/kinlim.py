"""Four-parton DIS points approaching single-unresolved limits, in
DISENT's layout/Breit frame (see kin.py). Built in the rest frame of
W = p1 + q (p1 along +z there) and boosted along z.
kind: 'fs:a,b' (partons a,b in {2,3,4} collinear), 'soft:a' (parton a
soft), 'is:a' (parton a collinear to the incoming parton)."""
import math, random, sys
from kin import boostz, write

def rot(n):
    """random unit vector"""
    c = 2*n.random()-1; ph = 2*math.pi*n.random(); s = math.sqrt(1-c*c)
    return [s*math.cos(ph), s*math.sin(ph), c]

def twobody(Ptot, rng, m1=0.0, dirn=None):
    """massless two-body decay of Ptot (4-vector) with random or given direction"""
    M2 = Ptot[3]**2 - sum(x*x for x in Ptot[:3]); M = math.sqrt(M2)
    d = dirn or rot(rng); E = M/2
    a = [E*d[0], E*d[1], E*d[2], E]; b = [-E*d[0], -E*d[1], -E*d[2], E]
    return boost(a, Ptot), boost(b, Ptot)

def boost(p, P):
    """boost p from the rest frame of P to the frame where P has momentum P"""
    M = math.sqrt(P[3]**2 - sum(x*x for x in P[:3]))
    bv = [P[i]/P[3] for i in range(3)]; b2 = sum(x*x for x in bv)
    if b2 < 1e-300: return p[:]
    g = P[3]/M; bp = sum(bv[i]*p[i] for i in range(3))
    q = [p[i] + ((g-1)*bp/b2 + g*p[3])*bv[i] for i in range(3)]
    return q + [g*(p[3] + bp)]

def point(kind, lam, rng, Q=40.0, y=0.4, xi=0.3):
    E = Q/2; E1 = E/xi; W = math.sqrt(4*E*E1 - 4*E**2); beta = (E1 - 2*E)/E1
    P = [[0.0]*4 for _ in range(8)]
    P[1] = [0, 0, E1, E1]; P[5] = [0, 0, -2*E, 0]
    P[6] = [E/y*2*math.sqrt(1-y), 0, -E, E/y*(2-y)]; P[7] = [E/y*2*math.sqrt(1-y), 0, E, E/y*(2-y)]
    Wv = [0, 0, 0, W]
    t, arg = kind.split(":")
    if t == "fs":
        a, b = map(int, arg.split(",")); c = ({2, 3, 4} - {a, b}).pop()
        M2 = lam*W*W; pc = (W*W - M2)/(2*W); d = rot(rng)
        fc = [pc*d[0], pc*d[1], pc*d[2], pc]
        Pab = [Wv[i]-fc[i] for i in range(4)]
        fa, fb = twobody(Pab, rng)
        f = {a: fa, b: fb, c: fc}
    elif t == "soft":
        a = int(arg); d = rot(rng); e = lam*W
        fa = [e*d[0], e*d[1], e*d[2], e]
        rest = [Wv[i]-fa[i] for i in range(4)]
        o = sorted({2, 3, 4} - {a}); fb, fc = twobody(rest, rng)
        f = {a: fa, o[0]: fb, o[1]: fc}
    elif t == "is":
        a = int(arg); th = math.sqrt(lam); e = 0.3*W; ph = 2*math.pi*rng.random()
        fa = [e*math.sin(th)*math.cos(ph), e*math.sin(th)*math.sin(ph), e*math.cos(th), e]
        rest = [Wv[i]-fa[i] for i in range(4)]
        o = sorted({2, 3, 4} - {a}); fb, fc = twobody(rest, rng)
        f = {a: fa, o[0]: fb, o[1]: fc}
    for i in (2, 3, 4):
        P[i] = boostz(f[i], beta)
    return P

if __name__ == "__main__":
    kind, seed, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    pts = []
    for lam in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        rng = random.Random(seed)
        pts.append(point(kind, lam, rng))
    write(out, pts)
