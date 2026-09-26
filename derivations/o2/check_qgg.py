import sys, subprocess, math
sys.path.insert(0, "harness")
import kin, formeval as fe
ex = fe.parse("qgg.out")
split = {}
for name, s in ex.items():
    # separate terms containing e_ (evaluate with E set to 0 vs full)
    split[name] = fe.compile_expr(s)
pts = [l.split() for l in open("harness/pts4.txt")]
for line in pts:
    v = list(map(float, line)); P = [None] + [v[4*i:4*i+4] for i in range(7)]
    mom = {"k": P[6], "kp": P[7], "p1": P[1], "p2": P[2], "p3": P[3], "p4": P[4]}
    env = fe.Env(mom)
    d = lambda a, b: fe.dot(mom[a], mom[b])
    S = {"D23": 1/(2*d("p2","p3")), "D24": 1/(2*d("p2","p4")), "D234": 1/(2*(d("p2","p3")+d("p2","p4")+d("p3","p4"))),
         "D14": -1/(2*d("p1","p4")), "D13": -1/(2*d("p1","p3")), "D134": 1/(2*(-d("p1","p3")-d("p1","p4")+d("p3","p4"))),
         "D34": 1/(2*d("p3","p4")), "i34": 1/d("p3","p4")}
    full = {n: fe.evaluate(c, env, S) for n, c in split.items()}
    env.E = {k: 0.0 for k in []} ; 
    class Z(dict):
        def __missing__(self, k): return 0.0
    env.E = Z()
    noeps = {n: fe.evaluate(c, env, S) for n, c in split.items()}
    print(" ".join(f"{n}={full[n]:.6e}(eps part {full[n]-noeps[n]:.1e})" for n in ("X1h77","X1h76","X2h77","X2h76")))
