"""Explicit Dirac-spinor evaluation of tree amplitudes (independent check).
Metric (+,-,-,-), Dirac representation; momenta given as (px,py,pz,E)."""
import numpy as np
I2 = np.eye(2); Z2 = np.zeros((2, 2))
sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
G = [np.block([[I2, Z2], [Z2, -I2]]).astype(complex)] + [np.block([[Z2, s], [-s, Z2]]) for s in sig]
G5 = 1j * G[0] @ G[1] @ G[2] @ G[3]
METRIC = np.diag([1, -1, -1, -1])

def up(p):   # (px,py,pz,E) -> contravariant (E,px,py,pz)
    return np.array([p[3], p[0], p[1], p[2]], float)

def slash(p):
    pu = up(p); return pu[0] * G[0] - pu[1] * G[1] - pu[2] * G[2] - pu[3] * G[3]

def spinors(p):
    """two independent massless spinors u(p) (columns) with u ubar summed = pslash;
    helicity eigenstates not needed since we sum over spins."""
    pu = up(p); E = pu[0]; ps = slash(p)
    # build from the projector pslash: take columns of pslash*gamma0 applied to basis
    M = ps @ G[0]
    # pick two independent columns, orthonormalise w.r.t. sum rule
    w, v = np.linalg.eigh((M + M.conj().T) / 2)
    # the non-zero eigenvalues are 2E (twice); eigenvectors give u with u u^dagger = 2E P
    idx = np.argsort(-w)[:2]
    return [v[:, i] * np.sqrt(w[i]) for i in idx]

def bar(u):
    return u.conj() @ G[0]

def vspinors(p):
    # for massless v(p) with sum v vbar = pslash as well: same construction
    return spinors(p)

def check():
    p = [0.3, -0.4, 1.2, np.sqrt(0.09 + 0.16 + 1.44)]
    s = sum(np.outer(u, bar(u)) for u in spinors(p))
    print("spin sum error", np.abs(s - slash(p)).max())
    for u in spinors(p): print("Dirac eq", np.abs(slash(p) @ u).max())
if __name__ == "__main__":
    check()
