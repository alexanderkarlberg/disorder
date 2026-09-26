"""FORM programs for the O(as^2) tree structures of DISENT (photon check
and helicity-resolved versions). See gen2.py for conventions."""
import sys
from gen2 import string, scalar, lin, slash

LEPTON = "(1/2)*g_(1,kp)*g_(1,mu)*g{l}_(1)*g_(1,k)*g_(1,nu)"
HEADER = ["Vectors k,kp,p1,p2,p3,p4,vv;", "Indices mu,nu,ac,bc,ad,bd,s1,s2,t1,t2;",
          "Symbols Dac,Dad,Dbc,Dbd,Dacd,Dbcd,Dcd,Dab,Dabc,Dabd,Dcab,Ddab,Dcb,Dad2,icd,Dx1,Dx2,Dx3,Dx4,Dx1b,Dx2b;"]
FOOT = ["trace4,1;", "trace4,2;", "trace4,3;", ".sort",
        "id kp = k + p1 - p2 - p3 - p4;",
        "id k.k = 0; id p1.p1 = 0; id p2.p2 = 0; id p3.p3 = 0; id p4.p4 = 0;",
        ".sort", "Format 255;", "Format nospaces;", "Print;", ".end"]

# ---- q qbar g g (generic labels) -------------------------------------
NAV = ("Dcd*(d_(ac,ad)*({pd:s1}-{pc:s1}) - d_(ad,s1)*({pc:ac}+2*{pd:ac})"
       " + d_(s1,ac)*(2*{pc:ad}+{pd:ad}))")
A1 = [("1", [("g","ac"),("S",[(1,"qa"),(1,"pc")],"Dac"),("g","ad"),("S",[(1,"qa"),(1,"pc"),(1,"pd")],"Dacd"),("g","mu")]),
      ("1", [("g","ac"),("S",[(1,"qa"),(1,"pc")],"Dac"),("g","mu"),("S",[(-1,"qb"),(-1,"pd")],"Dbd"),("g","ad")]),
      ("1", [("g","mu"),("S",[(-1,"qb"),(-1,"pc"),(-1,"pd")],"Dbcd"),("g","ac"),("S",[(-1,"qb"),(-1,"pd")],"Dbd"),("g","ad")]),
      (NAV, [("g","s1"),("S",[(1,"qa"),(1,"pc"),(1,"pd")],"Dacd"),("g","mu")]),
      (NAV, [("g","mu"),("S",[(-1,"qb"),(-1,"pc"),(-1,"pd")],"Dbcd"),("g","s1")])]
A2 = [("1", [("g","ad"),("S",[(1,"qa"),(1,"pd")],"Dad"),("g","ac"),("S",[(1,"qa"),(1,"pc"),(1,"pd")],"Dacd"),("g","mu")]),
      ("1", [("g","ad"),("S",[(1,"qa"),(1,"pd")],"Dad"),("g","mu"),("S",[(-1,"qb"),(-1,"pc")],"Dbc"),("g","ac")]),
      ("1", [("g","mu"),("S",[(-1,"qb"),(-1,"pc"),(-1,"pd")],"Dbcd"),("g","ad"),("S",[(-1,"qb"),(-1,"pc")],"Dbc"),("g","ac")]),
      ("-"+NAV, [("g","s1"),("S",[(1,"qa"),(1,"pc"),(1,"pd")],"Dacd"),("g","mu")]),
      ("-"+NAV, [("g","mu"),("S",[(-1,"qb"),(-1,"pc"),(-1,"pd")],"Dbcd"),("g","s1")])]
POL = "(-d_(ac,bc)+({pc:ac}*{pd:bc}+{pd:ac}*{pc:bc})*icd)*(-d_(ad,bd)+({pd:ad}*{pc:bd}+{pc:ad}*{pd:bd})*icd)"

def amp_sum(diags, line, sub, h, conj=False):
    terms = []
    for coef, fac in diags:
        f = fac + [("P", h)] if not conj else [("P", h)] + fac
        # projector sits next to the qb (v) end: for conj the reversed list puts it first
        terms.append(f"({scalar(coef, sub, conj)})*{string(fac, line, sub, conj)}")
    return "\n + ".join(terms)

def qgg(sub, sign, hels):
    out = ["* q qbar g g, labels " + str(sub)] + HEADER
    amps = {"1": A1, "2": A2}
    for tag, l, h in hels:
        for X in "12":
            for Y in "12":
                out.append(f"Local T{X}{Y}h{tag} = ({sign})*" + LEPTON.format(l=l) +
                           f"\n * {slash(lin([(1,'qa')], sub), 2)}*(\n {amp_sum(amps[X], 2, sub, h)}\n )*(1/2)*g{h}_(2)*"
                           f"{slash(lin([(1,'qb')], sub), 2)}*(\n {amp_sum(amps[Y], 2, sub, h, True)}\n ) * {scalar(POL, sub)};")
    out += FOOT[:-3]
    for tag, _, _ in hels:
        out += [f"Local X1h{tag} = T11h{tag} + T22h{tag};", f"Local X2h{tag} = T12h{tag} + T21h{tag};"]
    out += [".sort", "Drop " + ",".join(f"T{X}{Y}h{t}" for t, _, _ in hels for X in "12" for Y in "12") + ";"] + FOOT[-3:]
    return "\n".join(out) + "\n"

# ---- four quarks -------------------------------------------------------
# boson on the (qa,qb) line, gluon -> (qc,qd)
GA = [("Dcd", [("g","s1"),("S",[(1,"qa"),(1,"qc"),(1,"qd")],"Dx1"),("g","mu")]),
      ("Dcd", [("g","mu"),("S",[(-1,"qb"),(-1,"qc"),(-1,"qd")],"Dx2"),("g","s1")])]
GB = [("1", [("g","s1")])]

def fourq_line1(sub, sign, hels):
    """|A|^2 with the boson on the (qa,qb) line; second line (qc,qd) on spin line 3"""
    out = ["* four quarks, boson on (qa,qb), labels " + str(sub)] + HEADER
    for tag, l, h in hels:
        out.append(f"Local Dh{tag} = ({sign})*" + LEPTON.format(l=l) +
                   f"\n * {slash(lin([(1,'qa')], sub), 2)}*(\n {amp_sum(GA, 2, sub, h)}\n )*(1/2)*g{h}_(2)*"
                   f"{slash(lin([(1,'qb')], sub), 2)}*(\n {amp_sum(GA, 2, sub, h, True)}\n )"
                   f"\n * {slash(lin([(1,'qc')], sub), 3)}*g_(3,s1)*{slash(lin([(1,'qd')], sub), 3)}*g_(3,s2);")
    return "\n".join(out + FOOT) + "\n"

# identical quarks qa, qc: direct pairs (qa,qb),(qc,qd); exchanged (qc,qb),(qa,qd).
# Each amplitude: boson on either line; one long trace for A_dir A_exch*.
def e_amp_lines(x, y, u, w):
    """line strings for pairs (x,y) [boson] and (u,w) [gluon-produced]:
    returns (boson-line factors, other-line factors, gluon denominator symbol)"""
    return None

def fourq_exch(sub, sign, hels, attach=("xx", "xy", "yx", "yy")):
    """interference A_dir A_exch^*; attach selects the boson attachments:
    first letter: direct amplitude boson on (qa,qb) [x] or (qc,qd) [y];
    second letter: exchanged amplitude boson on (qc,qb) [x] or (qa,qd) [y]."""
    # boson-line strings with generic end labels
    def bline(o, i, g1, g2, den1, den2):
        # ubar(o) [ g^s S(o+g1+g2) g^mu + g^mu S(-i-g1-g2) g^s ] v(i)
        return [[("g","SI"),("S",[(1,o),(1,g1),(1,g2)],den1),("g","mu")],
                [("g","mu"),("S",[(-1,i),(-1,g1),(-1,g2)],den2),("g","SI")]]
    out = ["* identical-quark interference, labels " + str(sub) + " attach " + str(attach)] + HEADER
    for tag, l, h in hels:
        terms = []
        for att in attach:
            # direct amplitude
            if att[0] == "x":   # boson on (qa,qb), gluon (qc,qd): Dcd
                dir_ab = bline("qa","qb","qc","qd","Dx1","Dx2"); dir_cd = [[("g","SI")]]; gden1 = "Dcd"
            else:               # boson on (qc,qd), gluon (qa,qb): Dab
                dir_ab = [[("g","SI")]]; dir_cd = bline("qc","qd","qa","qb","Dx3","Dx4"); gden1 = "Dab"
            if att[1] == "x":   # exchanged: boson on (qc,qb), gluon (qa,qd)
                ex_cb = bline("qc","qb","qa","qd","Dcab","Dx1b"); ex_ad = [[("g","SI")]]; gden2 = "Dad2"
            else:               # boson on (qa,qd), gluon (qc,qb)
                ex_cb = [[("g","SI")]]; ex_ad = bline("qa","qd","qc","qb","Dabd","Dx2b"); gden2 = "Dcb"
            for fab in dir_ab:
                for fcd in dir_cd:
                    for fcb in ex_cb:
                        for fad in ex_ad:
                            s_ab = string([(f[0], f[1].replace("SI","s1")) if f[0]=="g" else f for f in fab], 2, sub)
                            s_cd = string([(f[0], f[1].replace("SI","s1")) if f[0]=="g" else f for f in fcd], 2, sub)
                            s_cb = string([(f[0], f[1].replace("SI","s2")) if f[0]=="g" else f for f in fcb], 2, sub, True)
                            s_ad = string([(f[0], f[1].replace("SI","s2")) if f[0]=="g" else f for f in fad], 2, sub, True)
                            # conj renames mu->nu, s1->s2 already (s2 stays s2)
                            terms.append(f"{gden1}*{gden2}*{slash(lin([(1,'qa')],sub),2)}*{s_ab}*(1/2)*g{h}_(2)*{slash(lin([(1,'qb')],sub),2)}"
                                         f"*{s_cb}*(1/2)*g{h}_(2)*{slash(lin([(1,'qc')],sub),2)}*{s_cd}*(1/2)*g{h}_(2)*{slash(lin([(1,'qd')],sub),2)}*{s_ad}")
        out.append(f"Local Eh{tag} = ({sign})*" + LEPTON.format(l=l) + " * (\n " + "\n + ".join(terms) + "\n );")
    return "\n".join(out + FOOT) + "\n"

HELS = [("77", 7, 7), ("76", 7, 6)]
if __name__ == "__main__":
    which = sys.argv[1]
    if which == "qgg_g":   # gluon-initiated: g(p1) -> q(p2) qbar(p3) g(p4)
        sub = {"qa": [(1,"p2")], "qb": [(1,"p3")], "pc": [(-1,"p1")], "pd": [(1,"p4")]}
        open("qgg_g.frm","w").write(qgg(sub, "1", HELS))
    elif which == "d_a":   # boson on incoming line 1 -> 2, gluon -> q(p3) qbar(p4)
        sub = {"qa": [(1,"p2")], "qb": [(-1,"p1")], "qc": [(1,"p3")], "qd": [(1,"p4")]}
        open("d_a.frm","w").write(fourq_line1(sub, "-1", HELS))
    elif which == "d_b":   # boson on pair q(p3) qbar(p4), gluon from line 1 -> 2
        sub = {"qa": [(1,"p3")], "qb": [(1,"p4")], "qc": [(1,"p2")], "qd": [(-1,"p1")]}
        open("d_b.frm","w").write(fourq_line1(sub, "-1", HELS))
    elif which == "e":     # identical: incoming q(p1) -> q(p2) q(p3) qbar(p4)
        sub = {"qa": [(1,"p2")], "qb": [(-1,"p1")], "qc": [(1,"p3")], "qd": [(1,"p4")]}
        open("e.frm","w").write(fourq_exch(sub, "-1", HELS))
    elif which == "e2":    # identical: incoming q(p1) -> q(p2) q(p4) qbar(p3)
        sub = {"qa": [(1,"p2")], "qb": [(-1,"p1")], "qc": [(1,"p4")], "qd": [(1,"p3")]}
        open("e2.frm","w").write(fourq_exch(sub, "-1", HELS))

# ---- CONTHR: q qbar g with the gluon polarisation sum replaced by vv vv
C3 = [("1", [("g","ac"),("S",[(1,"qa"),(1,"pc")],"Dac"),("g","mu")]),
      ("1", [("g","mu"),("S",[(-1,"qb"),(-1,"pc")],"Dbc"),("g","ac")])]
def conthr(sub, sign, hels):
    out = ["* spin-correlated q qbar g, labels " + str(sub)] + HEADER
    for tag, l, h in hels:
        out.append(f"Local Ch{tag} = ({sign})*" + LEPTON.format(l=l) +
                   f"\n * {slash(lin([(1,'qa')], sub), 2)}*(\n {amp_sum(C3, 2, sub, h)}\n )*(1/2)*g{h}_(2)*"
                   f"{slash(lin([(1,'qb')], sub), 2)}*(\n {amp_sum(C3, 2, sub, h, True)}\n ) * vv(ac)*vv(bc);")
    foot = [f if "id kp" not in f else "id kp = k + p1 - p2 - p3;" for f in FOOT]
    foot = [f.replace(" id p4.p4 = 0;", "") for f in foot]
    return "\n".join(out + foot) + "\n"

if __name__ == "__main__" and sys.argv[1] == "conthr_q":
    open("conthr_q.frm","w").write(conthr({"qa": [(1,"p2")], "qb": [(-1,"p1")], "pc": [(1,"p3")]}, "-1", HELS))
if __name__ == "__main__" and sys.argv[1] == "conthr_g":
    open("conthr_g.frm","w").write(conthr({"qa": [(1,"p2")], "qb": [(1,"p3")], "pc": [(-1,"p1")]}, "1", HELS))
