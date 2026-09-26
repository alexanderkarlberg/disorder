"""One-loop 0 -> qbar(1) q(2) l(3) lbar(4) g(5) helicity amplitudes of
Bern, Dixon, Kosower (Nucl. Phys. B513 (1998) 3), ported from MCFM
(src/W1jet/A51.f, A52.f, A5NLO.f, virt5.f; src/Need/spinoru.f, lnrat.f;
src/Wbb/lfunctions.f), R.K. Ellis and J. Campbell. Momenta are all
outgoing, (px,py,pz,E), negative energy for incoming particles.
Finite parts in MCFM's conventions (dimensional reduction, overall
(4pi)^eps Gamma(1+eps)Gamma(1-eps)^2/Gamma(1-2eps), poles dropped)."""
import cmath, math
from scipy.special import spence
PI = math.pi; PISQO6 = PI**2 / 6

def ddilog(x):
    """real dilogarithm Li2(x) for x <= 1"""
    return spence(1.0 - x)

def spinoru(p):
    """p: list of momenta indexed 1..N (p[0] unused). Returns za, zb, s."""
    N = len(p) - 1
    rt = [0] * (N + 1); c23 = [0] * (N + 1); f = [0] * (N + 1)
    for j in range(1, N + 1):
        if p[j][3] > 0:
            rt[j] = math.sqrt(p[j][3] + p[j][0]); c23[j] = complex(p[j][2], -p[j][1]); f[j] = 1
        else:
            rt[j] = math.sqrt(-p[j][3] - p[j][0]); c23[j] = complex(-p[j][2], p[j][1]); f[j] = 1j
    za = [[0j] * (N + 1) for _ in range(N + 1)]; zb = [[0j] * (N + 1) for _ in range(N + 1)]
    s = [[0.0] * (N + 1) for _ in range(N + 1)]
    for i in range(2, N + 1):
        for j in range(1, i):
            s[i][j] = 2 * (p[i][3]*p[j][3] - p[i][0]*p[j][0] - p[i][1]*p[j][1] - p[i][2]*p[j][2])
            za[i][j] = f[i] * f[j] * (c23[i] * rt[j] / rt[i] - c23[j] * rt[i] / rt[j])
            zb[i][j] = -s[i][j] / za[i][j]
            za[j][i] = -za[i][j]; zb[j][i] = -zb[i][j]; s[j][i] = s[i][j]
    return za, zb, s

def lnrat(x, y):
    """log(x - i ep) - log(y - i ep)"""
    th = lambda v: 1.0 if v > 0 else 0.0
    return complex(math.log(abs(x / y)), -PI * (th(-x) - th(-y)))

def L0(x, y):
    d = 1 - x / y
    if abs(d) < 1e-7: return -1 - d * (0.5 + d / 3)
    return lnrat(x, y) / d

def L1(x, y):
    d = 1 - x / y
    if abs(d) < 1e-7: return -0.5 - d / 3 * (1 + 0.75 * d)
    return (L0(x, y) + 1) / d

def Lsm1(x1, y1, x2, y2):
    r1 = x1 / y1; r2 = x2 / y2; o1 = 1 - r1; o2 = 1 - r2
    d1 = (PISQO6 - ddilog(r1)) - lnrat(x1, y1) * math.log(o1) if o1 > 1 else ddilog(o1)
    d2 = (PISQO6 - ddilog(r2)) - lnrat(x2, y2) * math.log(o2) if o2 > 1 else ddilog(o2)
    return d1 + d2 + lnrat(x1, y1) * lnrat(x2, y2) - PISQO6

class Amp:
    def __init__(self, za, zb, s, musq, epinv=0.0, epinv2=0.0):
        self.za, self.zb, self.s, self.musq, self.ep, self.ep2 = za, zb, s, musq, epinv, epinv2

    def A51(self, j1, j2, j3, j4, j5, za, zb):
        s, musq, ep, ep2 = self.s, self.musq, self.ep, self.ep2
        A5lom = -za[j3][j4]**2 / (za[j1][j2] * za[j2][j3] * za[j4][j5])
        l12 = lnrat(musq, -s[j1][j2]); l23 = lnrat(musq, -s[j2][j3])
        Vcc = -(ep2 + ep*l12 + 0.5*l12**2) - (ep2 + ep*l23 + 0.5*l23**2) - 2*(ep + l23) - 4
        Fcc = za[j3][j4]**2 / (za[j1][j2]*za[j2][j3]*za[j4][j5]) * (
            Lsm1(-s[j1][j2], -s[j4][j5], -s[j2][j3], -s[j4][j5])
            - 2*za[j3][j1]*zb[j1][j5]*za[j5][j4]/za[j3][j4] * L0(-s[j2][j3], -s[j4][j5]) / s[j4][j5])
        Vsc = 0.5*(ep + l23) + 1
        Fsc = (za[j3][j4]*za[j3][j1]*zb[j1][j5]*za[j5][j4] / (za[j1][j2]*za[j2][j3]*za[j4][j5])
               * L0(-s[j2][j3], -s[j4][j5]) / s[j4][j5]
               + 0.5*(za[j3][j1]*zb[j1][j5])**2*za[j4][j5] / (za[j1][j2]*za[j2][j3])
               * L1(-s[j2][j3], -s[j4][j5]) / s[j4][j5]**2)
        return (Vcc + Vsc)*A5lom + Fcc + Fsc

    def A52(self, j1, j2, j3, j4, j5, za, zb):
        s, musq, ep, ep2 = self.s, self.musq, self.ep, self.ep2
        l12 = lnrat(musq, -s[j1][j2]); l45 = lnrat(musq, -s[j4][j5])
        A5lom = za[j2][j4]**2 / (za[j2][j3]*za[j3][j1]*za[j4][j5])
        Vcc = -(ep2 + ep*l12 + 0.5*l12**2) - 2*(ep + l45) - 4
        Fcc = (-za[j2][j4]**2/(za[j2][j3]*za[j3][j1]*za[j4][j5]) * Lsm1(-s[j1][j2], -s[j4][j5], -s[j1][j3], -s[j4][j5])
               + za[j2][j4]*(za[j1][j2]*za[j3][j4] - za[j1][j4]*za[j2][j3]) / (za[j2][j3]*za[j1][j3]**2*za[j4][j5])
               * Lsm1(-s[j1][j2], -s[j4][j5], -s[j2][j3], -s[j4][j5])
               + 2*zb[j1][j3]*za[j1][j4]*za[j2][j4]/(za[j1][j3]*za[j4][j5]) * L0(-s[j2][j3], -s[j4][j5]) / s[j4][j5])
        Vsc = 0.5*(ep + l45) + 0.5
        Fsc = (za[j1][j4]**2*za[j2][j3]/(za[j1][j3]**3*za[j4][j5]) * Lsm1(-s[j1][j2], -s[j4][j5], -s[j2][j3], -s[j4][j5])
               - 0.5*(za[j4][j1]*zb[j1][j3])**2*za[j2][j3]/(za[j1][j3]*za[j4][j5]) * L1(-s[j4][j5], -s[j2][j3]) / s[j2][j3]**2
               + za[j1][j4]**2*za[j2][j3]*zb[j3][j1]/(za[j1][j3]**2*za[j4][j5]) * L0(-s[j4][j5], -s[j2][j3]) / s[j2][j3]
               - za[j2][j1]*zb[j1][j3]*za[j4][j3]*zb[j3][j5]/za[j1][j3] * L1(-s[j4][j5], -s[j1][j2]) / s[j1][j2]**2
               - za[j2][j1]*zb[j1][j3]*za[j3][j4]*za[j1][j4]/(za[j1][j3]**2*za[j4][j5]) * L0(-s[j4][j5], -s[j1][j2]) / s[j1][j2]
               - 0.5*zb[j3][j5]*(zb[j1][j3]*zb[j2][j5] + zb[j2][j3]*zb[j1][j5]) / (zb[j1][j2]*zb[j2][j3]*za[j1][j3]*zb[j4][j5]))
        return (Vcc + Vsc)*A5lom + Fcc + Fsc

    def A5NLO(self, j1, j2, j3, j4, j5, za, zb):
        """returns tree, leading-colour and subleading-colour (coefficient of 1/N^2) one-loop"""
        lo = -zb[j1][j4]**2 / (zb[j2][j5]*zb[j5][j1]*zb[j4][j3])
        return lo, self.A51(j2, j5, j1, j4, j3, zb, za), self.A52(j2, j1, j5, j4, j3, zb, za)

    def virt5(self, ip):
        """tree |A|^2, Re(A_lo^* A51), Re(A_lo^* A52) summed over gluon helicities
        for 0 -> qbar_R(ip1) q_L(ip2) l_L(ip3) a_R(ip4) g(ip5)"""
        za, zb = self.za, self.zb
        i1, i2, i3, i4, i5 = ip
        lom, a51m, a52m = self.A5NLO(i1, i2, i3, i4, i5, za, zb)
        lop, a51p, a52p = self.A5NLO(i2, i1, i4, i3, i5, zb, za)
        tree = abs(lom)**2 + abs(lop)**2
        v51 = (lom.conjugate()*a51m).real + (lop.conjugate()*a51p).real
        v52 = (lom.conjugate()*a52m).real + (lop.conjugate()*a52p).real
        return tree, v51, v52
