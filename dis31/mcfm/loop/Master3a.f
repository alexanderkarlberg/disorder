!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Master3a(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.7.4
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: Master3a,zab2,zab,zba
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: Lnrat,M3abit1,M3abit2
      real(dp):: t,Delta3,delta12,delta34,delta56
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
c     End statement functions
      delta12=s(j1,j2)-s(j5,j6)-s(j3,j4)
      delta34=s(j3,j4)-s(j5,j6)-s(j1,j2)
      delta56=s(j5,j6)-s(j1,j2)-s(j3,j4)
      Delta3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2*s(j1,j2)*s(j3,j4)-2*s(j3,j4)*s(j5,j6)-2*s(j5,j6)*s(j1,j2)
      M3abit1=
     & 0.5_dp*zb(j2,j4)/(zb(j2,j3)*s(j3,j4)*zab2(j1,j3,j4,j2))
     & *(6*s(j1,j2)*s(j3,j4)/Delta3**2
     & *(t(j1,j3,j4)-t(j2,j3,j4))*zab2(j3,j1,j2,j4)*zab2(j5,j1,j2,j6)

     & +za(j1,j2)*delta12/(zb(j5,j6)*Delta3)
     & *(zab2(j3,j1,j2,j6)*zb(j1,j4)*zb(j2,j6)
     &  +zab2(j3,j2,j4,j6)*zb(j1,j2)*zb(j4,j6)
     &  -za(j1,j3)*zb(j2,j4)*(zb(j1,j6)**2
     &  +zb(j2,j6)**2*zab2(j2,j3,j4,j1)/zab2(j1,j3,j4,j2)))

     & +zb(j1,j2)*delta12/(za(j5,j6)*Delta3)
     & *(zab(j5,j2,j4)*(za(j1,j2)*za(j3,j5)-za(j2,j3)*za(j1,j5))
     & -za(j1,j5)*za(j2,j5)*zab2(j3,j1,j2,j4)
     & +za(j1,j5)**2*za(j1,j3)*zb(j2,j4)*zab2(j2,j3,j4,j1)
     & /zab2(j1,j3,j4,j2))

     & -2*za(j1,j2)*za(j3,j5)/Delta3
     & *(zb(j1,j6)*zb(j2,j4)+zb(j1,j4)*zb(j2,j6))*delta56
     & -(za(j3,j4)*zb(j4,j6)**2/zb(j5,j6)
     & + za(j3,j5)**2*zb(j3,j4)/za(j5,j6))*s(j1,j2)*t(j2,j3,j4)/Delta3

     & +zab(j3,j2,j6)*zb(j4,j6)/zb(j5,j6)

     & -za(j1,j3)*zb(j2,j4)/zab2(j1,j3,j4,j2)
     & *(zab(j5,j1,j2)*za(j2,j5)/za(j5,j6)
     &  +zba(j6,j1,j2)*zb(j2,j6)/zb(j5,j6)))

      M3abit2=
     & +0.5_dp*zb(j2,j4)/(zb(j2,j3)*zab2(j1,j3,j4,j2)*Delta3)
     & *(3._dp*delta56*(t(j1,j3,j4)-t(j2,j3,j4))
     & *zab2(j3,j1,j2,j4)*zab2(j5,j1,j2,j6)/Delta3

     & -(zab(j3,j2,j6)*zb(j4,j6)/zb(j5,j6)
     &  -zab(j5,j1,j4)*za(j3,j5)/za(j5,j6)
     & +2*za(j1,j5)*zb(j2,j4)/(za(j5,j6)*zab2(j1,j3,j4,j2))
     & *(zab2(j3,j2,j4,j1)*za(j1,j5)-zab(j3,j2,j6)*za(j6,j5)))*delta34

     & +(zab2(j3,j1,j2,j4)
     & +4*za(j1,j3)*zb(j2,j4)*s(j5,j6)/zab2(j1,j3,j4,j2))
     & *(zab(j5,j1,j2)*za(j2,j5)/za(j5,j6)
     &  -zba(j6,j1,j2)*zb(j2,j6)/zb(j5,j6))

     & +8*s(j1,j2)*za(j1,j5)*zb(j2,j4)*zab(j3,j5,j6)/zab2(j1,j3,j4,j2)

     & -2*(zab(j3,j1,j4)-zab(j3,j2,j4))*zab2(j5,j1,j2,j6)
     & +4*(zba(j4,j1,j2)*zb(j2,j6)*za(j3,j5)
     &    +zab(j5,j2,j1)*za(j1,j3)*zb(j4,j6)))

      Master3a=
     & +M3abit1*lnrat(-s(j1,j2),-s(j5,j6))
     & +M3abit2*lnrat(-s(j3,j4),-s(j5,j6))

      return
      end
