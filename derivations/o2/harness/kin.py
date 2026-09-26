"""Random DIS phase-space points in DISENT's layout and Breit frame:
P[i] = (px, py, pz, E) for i = 1..7 (1 incoming parton, 2..n outgoing
partons, 5 = q = k - k', 6/7 incoming/outgoing lepton), as in GENTWO."""
import math, random

def boostz(p, beta):
    g = 1 / math.sqrt(1 - beta**2)
    return [p[0], p[1], g * (p[2] + beta * p[3]), g * (p[3] + beta * p[2])]

def rambo(n, W, rng):
    """n massless momenta with total (0,0,0,W), flat (RAMBO)."""
    qs = []
    for _ in range(n):
        c = 2 * rng.random() - 1; ph = 2 * math.pi * rng.random()
        e = -math.log(rng.random() * rng.random())
        s = math.sqrt(1 - c * c)
        qs.append([e * s * math.cos(ph), e * s * math.sin(ph), e * c, e])
    Q = [sum(q[i] for q in qs) for i in range(4)]
    M = math.sqrt(Q[3]**2 - Q[0]**2 - Q[1]**2 - Q[2]**2)
    b = [-Q[i] / M for i in range(3)]; g = Q[3] / M; a = 1 / (1 + g); x = W / M
    out = []
    for q in qs:
        bq = sum(b[i] * q[i] for i in range(3))
        out.append([x * (q[i] + b[i] * q[3] + a * bq * b[i]) for i in range(3)] + [x * (g * q[3] + bq)])
    return out

def point(n, rng, Q=None, y=None, xi=None):
    """n = number of partons including the incoming one (3 or 4)."""
    Q = Q or 10 ** rng.uniform(1.1, 2.0)
    y = y or rng.uniform(0.1, 0.9)
    xi = xi or rng.uniform(0.05, 0.9)
    E = Q / 2; E1 = E / xi
    P = [[0.0] * 4 for _ in range(8)]
    P[1] = [0, 0, E1, E1]
    P[5] = [0, 0, -2 * E, 0]
    P[6] = [E / y * 2 * math.sqrt(1 - y), 0, -E, E / y * (2 - y)]
    P[7] = [E / y * 2 * math.sqrt(1 - y), 0, E, E / y * (2 - y)]
    W = math.sqrt(4 * E * E1 - 4 * E**2)
    beta = (E1 - 2 * E) / E1
    fs = rambo(n - 1, W, rng)
    # random rotation of the final state about z keeps the leptons in x-z
    for i, p in enumerate(fs):
        P[2 + i] = boostz(p, beta)
    return P

def write(fname, pts):
    with open(fname, "w") as f:
        for P in pts:
            f.write(" ".join(f"{P[i][j]:.17e}" for i in range(1, 8) for j in range(4)) + "\n")

if __name__ == "__main__":
    import sys
    n, N, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(seed)
    write(sys.argv[4], [point(n, rng) for _ in range(N)])
