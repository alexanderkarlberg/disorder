"""e+ e- -> q(P1) qbar(P2) g(P3) points in DISENT layout: P5 = q (timelike),
P6 = electron, P7 = positron."""
import math, random, sys
from kin import rambo, write
rng = random.Random(int(sys.argv[2]))
pts = []
for _ in range(int(sys.argv[1])):
    rs = 10 ** rng.uniform(1.2, 2.0)
    fs = rambo(3, rs, rng)
    c = 2*rng.random()-1; ph = 2*math.pi*rng.random(); s_ = math.sqrt(1-c*c)
    e = [rs/2*s_*math.cos(ph), rs/2*s_*math.sin(ph), rs/2*c, rs/2]
    eb = [-e[0], -e[1], -e[2], rs/2]
    P = [None, fs[0], fs[1], fs[2], [0,0,0,0], [0,0,0,rs], e, eb]
    pts.append([None]+[P[i] for i in range(1,8)])
write(sys.argv[3], pts)
