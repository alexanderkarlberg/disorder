"""Lepton-direction average (exact for quadratic dependence) of the MCFM
photon one-loop interference in e+e- kinematics, which projects onto
-g_{mu nu} H^{mu nu}: compare with DISENT's ERTV alone."""
import sys, subprocess, math, numpy as np
sys.path.insert(0, "oneloop")
from bdk import spinoru, Amp
pts = "harness/ptsee.txt"
def harness(cf, ca):
    out = subprocess.run(["harness/harness", "virtee", pts, str(cf), str(ca), "0.5"], capture_output=True, text=True).stdout
    return [list(map(float, l.split())) for l in out.strip().split("\n")]
h10, h11 = harness(1, 0), harness(1, 1)
rows = []
for line, a, b in zip(open(pts), h10, h11):
    v = list(map(float, line.split())); P = [None] + [v[4*i:4*i+4] for i in range(7)]
    s = a[4]; E = math.sqrt(s)/2
    EX, EY = a[0], b[0]-a[0]
    T = I51 = I52 = 0.0
    R = np.linalg.qr(np.array([[0.3,0.5,0.8],[0.9,-0.2,0.1],[0.1,0.7,-0.4]]))[0]
    for d in [sg*R[:, i] for i in range(3) for sg in (1, -1)]:
        e = [E*d[0], E*d[1], E*d[2], E]; eb = [-E*d[0], -E*d[1], -E*d[2], E]
        neg = lambda x: [-t for t in x]
        za, zb, ss = spinoru([None, P[2], P[1], neg(eb), neg(e), P[3]])
        A = Amp(za, zb, ss, musq=s)
        for ip in ((1,2,3,4,5), (1,2,4,3,5)):
            t, x51, x52 = A.virt5(ip); T += t/6; I51 += x51/6; I52 += x52/6
    dd = lambda i, j: 2*(P[i][3]*P[j][3]-P[i][0]*P[j][0]-P[i][1]*P[j][1]-P[i][2]*P[j][2])
    L = [math.log(dd(1, 2)/s), math.log(dd(1, 3)/s), math.log(dd(2, 3)/s)]
    rows.append(dict(EX=EX, EY=EY, T=T, I51=I51, I52=I52, L=L))
for tgt in ("I51", "I52"):
    X = np.array([[r["EX"], r["EY"], r["T"]] + [r["T"]*l for l in r["L"]] +
                  [r["T"]*r["L"][i]*r["L"][j] for i in range(3) for j in range(i, 3)] for r in rows])
    y = np.array([r[tgt] for r in rows])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    print(tgt, "ERTV coeffs", np.round(c[:2], 8), "tree", np.round(c[2:], 5), "max rel resid", np.max(np.abs(X @ c - y)/np.abs(y)))
