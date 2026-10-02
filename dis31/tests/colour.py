# Colour matrices <c_m| T_i.T_k |c_n> for the DIS 3+1 colour bases, all partons
# outgoing (an incoming quark is an outgoing antiquark), CS conventions:
# quark: T^a, antiquark: -T^a^T, gluon: (T^a)_{bc} = -i f_{abc}.
import numpy as np, itertools
from fractions import Fraction
lam = np.zeros((8,3,3), complex)
lam[0][0,1]=lam[0][1,0]=1; lam[1][0,1]=-1j; lam[1][1,0]=1j; lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1; lam[4][0,2]=-1j; lam[4][2,0]=1j; lam[5][1,2]=lam[5][2,1]=1
lam[6][1,2]=-1j; lam[6][2,1]=1j; lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T = lam/2
f = np.zeros((8,8,8))
for a in range(8):
    for b in range(8):
        comm = T[a]@T[b]-T[b]@T[a]
        for c in range(8):
            f[a,b,c] = (-2j*np.trace(comm@T[c])).real   # [Ta,Tb] = i f_abc Tc
Fadj = -1j*f            # (F^a)_{bc} = -i f_abc
def gen(kind):
    if kind=='q': return T
    if kind=='qb': return -np.transpose(T,(0,2,1))
    if kind=='g': return Fadj
def apply(c, kinds, i, a):
    # apply colour generator a of parton i to tensor c (index i of c)
    G = gen(kinds[i])[a]
    return np.moveaxis(np.tensordot(G, c, axes=([1],[i])), 0, i)
def TT(c, kinds, i, k):
    out = np.zeros_like(c)
    for a in range(8):
        out += apply(apply(c, kinds, k, a), kinds, i, a)
    return out
def mat(basis, kinds, i, k):
    n=len(basis); M=np.zeros((n,n))
    for m in range(n):
        for l in range(n):
            M[m,l] = np.vdot(basis[m], TT(basis[l], kinds, i, k)).real
    return M
def frac(x): return str(Fraction(x).limit_denominator(1000))
# q qbar g g: indices (q, qb, g5, g6); c1 = (T^a5 T^a6)_{q qb}, c2 = (T^a6 T^a5)_{q qb}
kinds = ['q','qb','g','g']
c1 = np.einsum('aij,bjk->ikab', T, T)
c2 = np.einsum('bij,ajk->ikab', T, T)
basis=[c1,c2]
G0 = np.array([[np.vdot(x,y).real for y in basis] for x in basis])
print('qqbgg: <c|c> =', [[frac(x) for x in r] for r in G0])
names=['q','qb','g5','g6']
for i,k in itertools.combinations(range(4),2):
    M = mat(basis, kinds, i, k)
    print('qqbgg T_%s.T_%s:' % (names[i],names[k]), [[frac(x) for x in r] for r in M])
# colour conservation check
for i in range(4):
    s = sum(mat(basis,kinds,i,k) for k in range(4) if k!=i)
    Ci = {'q':4/3,'qb':4/3,'g':3}[kinds[i]]
    print('  check sum_k T_%s.T_k + C_i <c|c> =' % names[i], np.abs(s + Ci*G0).max())
# four quarks: (q1, qb2, q3, qb4); direct D = T^a_{q1 qb2} T^a_{q3 qb4}, exchange E = T^a_{q3 qb2} T^a_{q1 qb4}
kinds4 = ['q','qb','q','qb']
D = np.einsum('aij,akl->ijkl', T, T)
E = np.einsum('akj,ail->ijkl', T, T)
basis4=[D,E]
G4 = np.array([[np.vdot(x,y).real for y in basis4] for x in basis4])
print('4q: <c|c> =', [[frac(x) for x in r] for r in G4])
names4=['q1','qb2','q3','qb4']
for i,k in itertools.combinations(range(4),2):
    M = mat(basis4, kinds4, i, k)
    print('4q T_%s.T_%s:' % (names4[i],names4[k]), [[frac(x) for x in r] for r in M])
for i in range(4):
    s = sum(mat(basis4,kinds4,i,k) for k in range(4) if k!=i)
    print('  check', names4[i], np.abs(s + 4/3*G4).max())
