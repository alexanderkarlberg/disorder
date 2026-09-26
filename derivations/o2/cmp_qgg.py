import sys, subprocess
sys.path.insert(0, "harness")
import formeval as fe
ex = {n: fe.compile_expr(s) for n, s in fe.parse("qgg.out").items()}
def form_vals(ptsfile):
    res = []
    for line in open(ptsfile):
        v = list(map(float, line.split())); P = [None] + [v[4*i:4*i+4] for i in range(7)]
        mom = {"k": P[6], "kp": P[7], "p1": P[1], "p2": P[2], "p3": P[3], "p4": P[4]}
        env = fe.Env(mom); d = lambda a, b: fe.dot(mom[a], mom[b])
        S = {"D23": 1/(2*d("p2","p3")), "D24": 1/(2*d("p2","p4")), "D234": 1/(2*(d("p2","p3")+d("p2","p4")+d("p3","p4"))),
             "D14": -1/(2*d("p1","p4")), "D13": -1/(2*d("p1","p3")), "D134": 1/(2*(-d("p1","p3")-d("p1","p4")+d("p3","p4"))),
             "D34": 1/(2*d("p3","p4")), "i34": 1/d("p3","p4")}
        res.append({n: fe.evaluate(c, env, S) for n, c in ex.items()})
    return res
F = form_vals(sys.argv[1])
for CF, CA in ((4/3, 3.0), (1.0, 0.0), (0.7, 2.2)):
    out = subprocess.run(["harness/harness", "matfor", sys.argv[1], str(CF), str(CA), "0.5"], capture_output=True, text=True).stdout.split("\n")
    print(f"CF={CF:.3f} CA={CA}: ratio DISENT/FORM (quark q->qgg):")
    for f, line in zip(F, out):
        QG, EMSQ = float(line.split()[0]), float(line.split()[5])
        mine = CF*(f["X1h77"]+f["X1h76"]) + (CF-CA/2)*(f["X2h77"]+f["X2h76"])
        print(f"   {QG/mine:.15f}   (EMSQ {EMSQ:.1f})")
