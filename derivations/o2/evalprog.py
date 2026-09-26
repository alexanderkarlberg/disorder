"""Numerical values of FORM program outputs at DISENT-layout points."""
import sys
sys.path.insert(0, "harness")
import formeval as fe

MAPS = {
 "qgg":   {"qa": [(1,"p2")], "qb": [(-1,"p1")], "pc": [(1,"p3")], "pd": [(1,"p4")]},
 "qgg_g": {"qa": [(1,"p2")], "qb": [(1,"p3")], "pc": [(-1,"p1")], "pd": [(1,"p4")]},
 "d_a":   {"qa": [(1,"p2")], "qb": [(-1,"p1")], "qc": [(1,"p3")], "qd": [(1,"p4")]},
 "d_b":   {"qa": [(1,"p3")], "qb": [(1,"p4")], "qc": [(1,"p2")], "qd": [(-1,"p1")]},
 "e":     {"qa": [(1,"p2")], "qb": [(-1,"p1")], "qc": [(1,"p3")], "qd": [(1,"p4")]},
 "e2":    {"qa": [(1,"p2")], "qb": [(-1,"p1")], "qc": [(1,"p4")], "qd": [(1,"p3")]},
}
DENS = {  # symbol: labels summed
 "Dac": "qa pc", "Dad": "qa pd", "Dbc": "qb pc", "Dbd": "qb pd", "Dacd": "qa pc pd",
 "Dbcd": "qb pc pd", "Dcd": "pc pd",
}
DENS4 = {"Dcd": "qc qd", "Dab": "qa qb", "Dx1": "qa qc qd", "Dx2": "qb qc qd", "Dx3": "qc qa qb",
         "Dx4": "qd qa qb", "Dcab": "qc qa qd", "Dx1b": "qb qa qd", "Dad2": "qa qd",
         "Dabd": "qa qc qb", "Dx2b": "qd qc qb", "Dcb": "qc qb"}

def labmom(lab, mp, P):
    v = [0.0]*4
    for s, p in mp[lab]:
        i = int(p[1:])
        for j in range(4): v[j] += s * P[i][j]
    return v

def symbols(prog, P):
    mp = MAPS[prog]
    S = {}
    table = DENS if prog.startswith("qgg") else DENS4
    for name, labs in table.items():
        labs = labs.split()
        if not all(l in mp for l in labs): continue
        v = [sum(labmom(l, mp, P)[j] for l in labs) for j in range(4)]
        S[name] = 1.0 / fe.dot(v, v)
    if prog.startswith("qgg"):
        S["icd"] = 1.0 / fe.dot(labmom("pc", mp, P), labmom("pd", mp, P))
    return S

def values(prog, ptsfile, outfile=None):
    ex = {n: fe.compile_expr(s) for n, s in fe.parse(outfile or prog + ".out").items()}
    res = []
    for line in open(ptsfile):
        v = list(map(float, line.split())); P = [None] + [v[4*i:4*i+4] for i in range(7)]
        mom = {"k": P[6], "kp": P[7], "p1": P[1], "p2": P[2], "p3": P[3], "p4": P[4]}
        env = fe.Env(mom)
        S = symbols(prog, P)
        res.append({n: fe.evaluate(c, env, S) for n, c in ex.items()})
    return res
