"""Size of the non-tree part of (standard - real-log continuation) of the
MCFM photon one-loop interference in DIS, compared with DISENT's
non-factorising part and with the tree, per colour structure."""
import sys, math, numpy as np, copy
sys.path.insert(0, "oneloop")
import bdk
from bdk import spinoru, Amp
import subprocess
pts = "harness/pts3v.txt"
def harness(cf, ca):
    out = subprocess.run(["harness/harness", "virt", pts, str(cf), str(ca), "0.5"], capture_output=True, text=True).stdout
    return [list(map(float, l.split())) for l in out.strip().split("\n")]
h10, h11 = harness(1, 0), harness(1, 1)
std_lnrat = bdk.lnrat
def realln(x, y): return complex(math.log(abs(x/y)), 0.0)
rows = []
for line, a, b in zip(open(pts), h10, h11):
    v = list(map(float, line.split())); P = [None] + [v[4*i:4*i+4] for i in range(7)]
    emsq = a[4]
    NX = (2*a[1] - emsq/2*a[0])/emsq; NY = (2*(b[1]-a[1]) - emsq/2*(b[0]-a[0]))/emsq
    neg = lambda x: [-t for t in x]
    za, zb, s = spinoru([None, neg(P[1]), P[2], P[7], neg(P[6]), P[3]])
    r = {}
    for tag, f in (("std", std_lnrat), ("re", realln)):
        bdk.lnrat = f
        A = Amp(za, zb, s, musq=emsq)
        tLL, x1, y1 = A.virt5((1, 2, 3, 4, 5)); tLR, x2, y2 = A.virt5((1, 2, 4, 3, 5))
        r[tag] = (tLL+tLR, x1+x2, y1+y2)
    bdk.lnrat = std_lnrat
    d = lambda i, j: 2*(P[i][3]*P[j][3]-P[i][0]*P[j][0]-P[i][1]*P[j][1]-P[i][2]*P[j][2])
    L = [math.log(abs(d(1, 2))/emsq), math.log(abs(d(1, 3))/emsq), math.log(abs(d(2, 3))/emsq)]
    rows.append(dict(T=r["std"][0], D51=r["std"][1]-r["re"][1], D52=r["std"][2]-r["re"][2], L=L,
                     N51=0.00633257*(NX+2*NY), N52=-0.00633257*NX))
for k, n in (("D51", "N51"), ("D52", "N52")):
    X = np.array([[r["T"]] + [r["T"]*l for l in r["L"]] + [r["T"]*r["L"][i]*r["L"][j] for i in range(3) for j in range(i, 3)] for r in rows])
    y = np.array([r[k] for r in rows]); c, *_ = np.linalg.lstsq(X, y, rcond=None)
    R = y - X @ c
    print(k, "tree-part coeffs", np.round(c, 4))
    print("   non-tree residual / DISENT non-fact:", np.round(R / np.array([abs(r[n]) for r in rows]), 3)[:10])
    print("   non-tree residual / tree:", np.round(R / np.array([r["T"] for r in rows]), 3)[:10])
