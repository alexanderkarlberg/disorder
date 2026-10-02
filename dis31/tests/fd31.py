#!/usr/bin/env python3
"""Independent check of me31's four-quark channels (photon exchange):
Feynman diagrams with explicit Dirac matrices and helicity spinors, physical
momenta (no crossing), explicit SU(3) colour sums, for
    e(k) q(p1) -> e(k') q(p2) Q(p3) Qbar(p4)
(photon on the incoming line, two diagrams, and on the pair, two diagrams;
identical quarks: minus the same with p2 <-> p3). Reads fd31_points.txt
(written by dump_fd31) and requires me31/FD = 1 to 1e-10.
Normalisation as me31: average over lepton and quark spin and colour,
alpha = 1/137, divided by (alpha_s/2pi)^2."""
import itertools, sys
import numpy as np

I2, Z2 = np.eye(2), np.zeros((2, 2))
sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]),
       np.array([[1, 0], [0, -1]], complex)]
G = [np.block([[I2, Z2], [Z2, -I2]]).astype(complex)] + \
    [np.block([[Z2, s], [-s, Z2]]) for s in sig]          # gamma^mu, Dirac rep.

def dot(a, b):
    return a[0]*b[0] - a[1:] @ b[1:]

def slash(p):                                             # p = (E, px, py, pz)
    return p[0]*G[0] - p[1]*G[1] - p[2]*G[2] - p[3]*G[3]

def u(p, h):                     # massless spinor of helicity h
    th = np.arccos(np.clip(p[3]/p[0], -1, 1)); ph = np.arctan2(p[2], p[1])
    c = (np.array([np.cos(th/2), np.exp(1j*ph)*np.sin(th/2)]) if h > 0 else
         np.array([-np.exp(-1j*ph)*np.sin(th/2), np.cos(th/2)]))
    return np.sqrt(p[0])*np.concatenate([c, h*c])

def v(p, h):                     # antiquark spinor (up to a phase)
    return u(p, -h)

def bar(s):
    return s.conj() @ G[0]

def vslash(sb, s):               # gamma^mu (sb gamma_mu s)
    j = np.array([sb @ G[m] @ s for m in range(4)])
    return slash(j)

T = np.zeros((8, 3, 3), complex)  # Gell-Mann / 2
T[0][0, 1] = T[0][1, 0] = 1; T[1][0, 1] = -1j; T[1][1, 0] = 1j
T[2][0, 0] = 1; T[2][1, 1] = -1
T[3][0, 2] = T[3][2, 0] = 1; T[4][0, 2] = -1j; T[4][2, 0] = 1j
T[5][1, 2] = T[5][2, 1] = 1; T[6][1, 2] = -1j; T[6][2, 1] = 1j
T[7] = np.diag([1, 1, -2])/np.sqrt(3)
T /= 2
cD = np.einsum('aji,akl->ijkl', T, T)   # T^a_{c2 c1} T^a_{c3 c4}
cE = np.einsum('aki,ajl->ijkl', T, T)   # p2 <-> p3
CDD = np.sum(abs(cD)**2).real
CDE = np.sum(cD*cE.conj()).real

def amp(p1, p2, p3, p4, k, k2, h, eq, eQ):
    h1, h2, h3, h4, hk, hk2 = h
    q = k - k2
    L = vslash(bar(u(k2, hk2)), u(k, hk))/dot(q, q)
    # photon on the incoming line, gluon p3 + p4
    g = p3 + p4
    J = vslash(bar(u(p3, h3)), v(p4, h4))/dot(g, g)
    a = bar(u(p2, h2)) @ (J @ slash(p1 + q) @ L/dot(p1 + q, p1 + q)
                          + L @ slash(p1 - g) @ J/dot(p1 - g, p1 - g)) @ u(p1, h1)
    # photon on the pair, gluon p1 - p2
    g = p1 - p2
    J = vslash(bar(u(p2, h2)), u(p1, h1))/dot(g, g)
    b = bar(u(p3, h3)) @ (L @ slash(p3 - q) @ J/dot(p3 - q, p3 - q)
                          + J @ slash(q - p4) @ L/dot(q - p4, q - p4)) @ v(p4, h4)
    return eq*a + eQ*b

def fd(P, eq, eQ, ident):
    four = lambda x: np.array([x[3], x[0], x[1], x[2]])
    p1, p2, p3, p4 = (four(P[i]) for i in range(4))
    k, k2 = four(P[5]), four(P[6])
    s = 0
    for h in itertools.product((1, -1), repeat=6):
        D = amp(p1, p2, p3, p4, k, k2, h, eq, eQ)
        if not ident:
            s += CDD*abs(D)**2
            continue
        E = amp(p1, p3, p2, p4, k, k2, (h[0], h[2], h[1], h[3], h[4], h[5]), eq, eQ)
        s += CDD*(abs(D)**2 + abs(E)**2) - 2*CDE*(D*np.conj(E)).real
    return s*(4*np.pi/137)**2*(8*np.pi**2)**2/(2*2*3)

worst = 0
for row in np.loadtxt(sys.argv[1] if len(sys.argv) > 1 else 'fd31_points.txt'):
    P = row[:28].reshape(7, 4)      # P[i] = (px, py, pz, E) of DISENT slot i+1
    m = row[28:]
    r = [m[0]/fd(P, -1/3, 2/3, False), m[1]/fd(P, 2/3, 2/3, True),
         m[2]/fd(P, -1/3, 2/3, False)]   # dbar: C-conjugate of d (photon)
    print('me31/FD: d->d u ubar %.13f  u->u u ubar %.13f  dbar->dbar ubar u %.13f' % tuple(r))
    worst = max(worst, max(abs(x - 1) for x in r))
print('largest |me31/FD - 1|: %.2e' % worst)
sys.exit(0 if worst < 1e-10 else 1)
