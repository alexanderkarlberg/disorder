"""Evaluate FORM output (Format 255, nospaces) numerically at DISENT-layout
phase-space points. Expressions are compiled into Python functions of a
dict of dot products, inverse propagators and Levi-Civita values."""
import re, math

def parse(fname):
    out = open(fname).read()
    ex = {}
    for m in re.finditer(r"^\s*(\w+)=\s*(.*?);", out, re.S | re.M):
        ex[m.group(1)] = re.sub(r"\s+", "", m.group(2))
    return ex

def split_terms(s):
    """split at top-level + and - signs"""
    terms, depth, start = [], 0, 0
    for i, ch in enumerate(s):
        if ch == "(": depth += 1
        elif ch == ")": depth -= 1
        elif ch in "+-" and depth == 0 and i > start and s[i-1] not in "*/^(":
            terms.append(s[start:i]); start = i
    terms.append(s[start:])
    return terms

def compile_expr(s, chunk=150):
    t = split_terms(s)
    return [_compile1("".join(t[i:i+chunk])) for i in range(0, len(t), chunk)]

def _compile1(s):
    s = s.replace("^", "**")
    s = re.sub(r"e_\((\w+),(\w+),(\w+),(\w+)\)", r"E['\1','\2','\3','\4']", s)
    s = re.sub(r"\b(\w+)\.(\w+)\b", r"V['\1.\2']", s)
    s = re.sub(r"(?<![\w'\[])([A-Za-z]\w*)(?![\w'(\[\]])", r"S['\1']", s)
    return compile(s, "<form>", "eval")

def dot(a, b):
    return a[3]*b[3] - a[0]*b[0] - a[1]*b[1] - a[2]*b[2]

def eps(a, b, c, d):
    """Levi-Civita eps_{mu nu rho sigma} a b c d with eps_{0123} = +1,
    vectors given as (px,py,pz,E); lower indices via metric (+,-,-,-)."""
    import itertools
    v = [[x[3], x[0], x[1], x[2]] for x in (a, b, c, d)]   # upper components
    lo = [[x[0], -x[1], -x[2], -x[3]] for x in v]            # lower components
    tot = 0.0
    for perm in itertools.permutations(range(4)):
        sgn = 1
        p = list(perm)
        for i in range(4):
            for j in range(i + 1, 4):
                if p[i] > p[j]: sgn = -sgn
        tot += sgn * lo[0][p[0]] * lo[1][p[1]] * lo[2][p[2]] * lo[3][p[3]]
    return tot

class Env:
    def __init__(self, mom):
        self.mom = mom
        names = list(mom)
        self.V = {}
        for a in names:
            for b in names:
                self.V[f"{a}.{b}"] = dot(mom[a], mom[b])
        self.E = _Eps(mom)

class _Eps(dict):
    def __init__(self, mom): self.mom = mom
    def __missing__(self, key):
        v = eps(*(self.mom[k] for k in key)); self[key] = v; return v

def evaluate(codes, env, S):
    g = {"V": env.V, "E": env.E, "S": S, "math": math}
    return sum(eval(c, g) for c in codes)
