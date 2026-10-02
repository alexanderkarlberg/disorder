!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Master2a(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.7.2
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: Master2a,zab2,zba2,zab,zba,Lnrat,M2abit1,M2abit2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t,delta12,delta34,delta56,Delta3
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
c     End statement functions

      delta12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      delta34=s(j3,j4)-s(j5,j6)-s(j1,j2)
      delta56=s(j5,j6)-s(j1,j2)-s(j3,j4)
      Delta3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2*s(j1,j2)*s(j3,j4)-2*s(j3,j4)*s(j5,j6)-2*s(j5,j6)*s(j1,j2)
      M2abit1=
     & 0.5_dp*zb(j1,j2)/(zb(j2,j3)*zab2(j4,j1,j2,j3))
     & *(6*za(j1,j2)*zab2(j3,j1,j2,j4)*zab2(j5,j1,j2,j6)
     & *(zb(j1,j3)*delta56-2*zba(j1,j2,j4)*zb(j4,j3))/Delta3**2

     & +1._dp/Delta3*(
     & zab(j2,j1,j4)*za(j3,j5)*zab2(j5,j1,j2,j3)/za(j5,j6)
     & -za(j1,j2)*za(j4,j5)/(za(j3,j4)*za(j5,j6))
     & *(zb(j1,j3)*za(j4,j5)*zab2(j3,j1,j2,j4)*delta12/zab2(j4,j1,j2,j3)
     & -za(j3,j5)*(2*zb(j1,j4)*delta12-zba(j1,j2,j3)*zb(j3,j4)))

     & -za(j2,j4)*za(j3,j5)*delta12/(za(j3,j4)*za(j5,j6))
     & *(zab2(j5,j1,j2,j4)-zab(j5,j6,j4))

     & +zb(j4,j6)/zb(j5,j6)
     & *(2*(zab(j2,j3,j6)*delta34-zab(j2,j5,j6) *delta56 )
     & +zab(j3,j4,j6)*zab2(j2,j1,j4,j3))

     & -zab(j2,j1,j3)*(zab(j5,j3,j6)-zab(j5,j4,j6))*delta12
     & /zab2(j4,j1,j2,j3)

     & -4*za(j1,j2)*zb(j3,j6)*(za(j3,j5)*zb(j1,j4)
     & +za(j3,j4)*za(j5,j6)*zb(j1,j3)*zb(j4,j6)/zab2(j4,j1,j2,j3))

     & +zb(j4,j6)*(zab(j2,j4,j6)*za(j6,j5)-zab(j2,j1,j3)*za(j3,j5)))


     & +zab2(j5,j1,j2,j3)
     & *(za(j2,j5)*za(j3,j4)-za(j2,j3)*za(j4,j5))
     & /(za(j3,j4)*za(j5,j6)*zab2(j4,j1,j2,j3)))

      M2abit2=
     & -0.5_dp/(zb(j2,j3)*zab2(j4,j1,j2,j3)*Delta3)
     & *(
     & -6*zb(j1,j2)*zab2(j5,j1,j2,j6)*zb(j4,j3)*zab2(j3,j1,j2,j4)
     & *(za(j2,j4)*delta56-2*zab(j2,j1,j3)*za(j3,j4))/Delta3

     & +zb(j1,j3)/zab2(j4,j1,j2,j3)
     & *(delta34*(t(j1,j2,j3)-t(j1,j2,j4))-Delta3)
     & *(zab(j5,j3,j4)*za(j4,j5)/za(j5,j6)
     &  +zba(j6,j3,j4)*zb(j4,j6)/zb(j5,j6))

     & +zab(j5,j3,j4)/za(j5,j6)
     & *(2*(zab(j5,j4,j1)*delta34-zab(j5,j6,j1)*delta56)
     & +zab2(j4,j2,j3,j1)*zab2(j5,j1,j2,j4))

     & -zb(j4,j6)/zb(j5,j6)
     & *(2*zab(j3,j2,j1)
     & *(zb(j3,j6)*delta12-2*zba(j3,j4,j5)*zb(j5,j6))
     & +zba2(j1,j2,j3,j4)*zb(j4,j3)*zab2(j3,j1,j2,j6)
     & -3*zba(j6,j4,j3)*zb(j3,j1)*delta34))

      Master2a=
     & +M2abit1*lnrat(-s(j1,j2),-s(j5,j6))
     & +M2abit2*lnrat(-s(j3,j4),-s(j5,j6))

      return
      end
