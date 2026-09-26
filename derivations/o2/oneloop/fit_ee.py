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
    s = a[4]
    NX = (2*a[1] + s/2*a[0])/s; NY = (2*(b[1]-a[1]) + s/2*(b[0]-a[0]))/s
    neg = lambda x: [-t for t in x]
    za, zb, ss = spinoru([None, P[2], P[1], neg(P[7]), neg(P[6]), P[3]])
    A = Amp(za, zb, ss, musq=s)
    tLL, v51LL, v52LL = A.virt5((1, 2, 3, 4, 5)); tLR, v51LR, v52LR = A.virt5((1, 2, 4, 3, 5))
    d = lambda i, j: 2*(P[i][3]*P[j][3]-P[i][0]*P[j][0]-P[i][1]*P[j][1]-P[i][2]*P[j][2])
    L = [math.log(d(1, 2)/s), math.log(d(1, 3)/s), math.log(d(2, 3)/s)]
    rows.append(dict(NX=NX, NY=NY, T=tLL+tLR, I51=v51LL+v51LR, I52=v52LL+v52LR, L=L))
for tgt in ("I51", "I52"):
    X = np.array([[r["NX"], r["NY"], r["T"], r["T"]*r["L"][0], r["T"]*r["L"][1], r["T"]*r["L"][2],
                   r["T"]*r["L"][0]**2, r["T"]*r["L"][1]**2, r["T"]*r["L"][2]**2] for r in rows])
    y = np.array([r[tgt] for r in rows])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    print(tgt, np.round(c, 8), "max rel resid", np.max(np.abs(X @ c - y)/np.abs(y)))
print("with cross logs:")
for tgt in ("I51", "I52"):
    X = np.array([[r["NX"], r["NY"], r["T"]] + [r["T"]*l for l in r["L"]] +
                  [r["T"]*r["L"][i]*r["L"][j] for i in range(3) for j in range(i, 3)] for r in rows])
    y = np.array([r[tgt] for r in rows])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    print(tgt, np.round(c[:2], 8), "max rel resid", np.max(np.abs(X @ c - y)/np.abs(y)))
