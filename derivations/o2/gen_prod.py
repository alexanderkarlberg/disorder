"""Production FORM programs: helicity-difference (parity-violating) and,
for CC, helicity-sum structures, written out as optimised Fortran.
Uses the amplitude builders of gen_o2.py (validated against DISENT)."""
import sys, re
import gen_o2 as g

DOTS = [("k","p1"),("k","p2"),("k","p3"),("k","p4"),("p1","p2"),("p1","p3"),("p1","p4"),
        ("p2","p3"),("p2","p4"),("p3","p4"),("k","vv"),("p1","vv"),("p2","vv"),("p3","vv"),("vv","vv")]

def production(body_frm, exprs, fname, subname, args, dens, vv=False):
    """body_frm: FORM program text up to (excluding) the final Print/.end;
    exprs: {fortran name: FORM expression of the Locals}."""
    lines = body_frm.rstrip().split("\n")
    lines = [l for l in lines if l.strip() not in ("Print;", ".end") and not l.startswith("Format")]
    dsyms = "".join(f"{a}{b}," for a, b in DOTS if vv or "vv" not in (a + b)).rstrip(",")
    out = lines + [".sort", "Symbols " + dsyms + ";", "Vectors v1,v2,v3,v4;"]
    for n, e in exprs.items():
        out.append(f"Local {n} = {e};")
    out += [".sort", "id e_(v1?,v2?,v3?,v4?) = 0;"]
    for a, b in DOTS:
        if not vv and "vv" in a + b: continue
        out.append(f"id {a}.{b} = {a}{b};")
    out += [".sort", "Drop;"] + [f"NDrop {n};" for n in exprs] + [".sort",
            "Format Fortran;", "Format doublefortran;", "Format O3,stats=on;"]
    for n in exprs:
        out += [f"#optimize {n}", f'#write <{fname}> "%O"', f'#write <{fname}> "      {n} = %e" {n}', "#clearoptimize"]
    out += [".end"]
    return "\n".join(out) + "\n"
