"""Wrap FORM's optimised Fortran output into a subroutine.
momenta are passed as (px,py,pz,E) arrays; dens: {symbol: [(sign, name), ...]}
(inverse of the square of the signed sum)."""
import re

def wrap(body_file, subname, momenta, dens, outputs, extra=None, doc=""):
    body = open(body_file).read()
    used = set(re.findall(r"\b([a-z]+\d*[a-z]*\d*)\b", body))
    names = momenta
    lines = [f"C {l}" for l in doc.strip().split("\n")] if doc else []
    lines += [f"      SUBROUTINE {subname}({','.join(names)},{','.join(outputs)})",
              "      IMPLICIT DOUBLE PRECISION (A-Z)",
              f"      DIMENSION {','.join(n + '(4)' for n in names)}",
              "C---dot products (metric +,-,-,-), DISENT layout (px,py,pz,E)"]
    def dotexpr(a, b):
        return f"{a}(4)*{b}(4)-{a}(1)*{b}(1)-{a}(2)*{b}(2)-{a}(3)*{b}(3)"
    for i, a in enumerate(names):
        for b in names[i:]:
            sym = f"{a}{b}"
            if sym in used:
                lines.append(f"      {sym}={dotexpr(a, b)}")
    for d, mom in dens.items():
        v = "+".join(f"{'-' if s < 0 else ''}{n}(I)" for s, n in mom).replace("+-", "-")
        lines += [f"      DO I=1,4", f"        TMP4(I)={v}", f"      ENDDO",
                  f"      {d}=1/(TMP4(4)**2-TMP4(1)**2-TMP4(2)**2-TMP4(3)**2)"]
    if dens:
        j = [i for i, l in enumerate(lines) if "IMPLICIT" in l][0] + 1
        lines.insert(j, "      DIMENSION TMP4(4)")
        lines.insert(j, "      INTEGER I")
    if extra:
        lines += extra
    lines.append(body.rstrip())
    lines += ["      END", ""]
    return "\n".join(lines)
