!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFPMccT(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.10.17
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FFPMccT,zab2,zab,zba
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t,delta12,delta34
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
c     End statement functions
      delta12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      delta34=s(j3,j4)-s(j5,j6)-s(j1,j2)
      FFPMccT=
     & -za(j2,j4)**2*zb(j1,j6)*zab2(j5,j2,j4,j3)
     & /(za(j2,j3)*t(j2,j3,j4)**2)
     & +0.5_dp*zb(j1,j3)*za(j2,j4)*za(j2,j5)*zab2(j5,j2,j4,j3)
     & /(za(j2,j3)*zb(j3,j4)*za(j5,j6)*t(j2,j3,j4))
     & +0.5_dp*za(j2,j4)**2*zb(j1,j6)*zab2(j4,j1,j3,j6)
     & /(za(j2,j3)*za(j3,j4)*zb(j5,j6)*t(j2,j3,j4))

     & -za(j4,j5)*zab(j2,j1,j3)*zab2(j2,j1,j3,j6)
     & /(za(j2,j3)*t(j1,j2,j3)*zab2(j1,j2,j3,j4))

     & +0.5_dp/(zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2))
     & *(4*zab(j5,j1,j6)*zab(j4,j2,j3)*zab2(j1,j2,j4,j3)/t(j2,j3,j4)

     & +zab2(j1,j2,j4,j3)*zab2(j5,j2,j4,j3)**2
     & *(2*s(j3,j4)*s(j5,j6)+delta12*t(j2,j3,j4))
     & /(zb(j3,j4)*za(j5,j6)*t(j2,j3,j4)**2)
     & +zab(j5,j1,j3)**2
     & *(za(j1,j2)*delta12+zab2(j1,j2,j3,j4)*za(j4,j2))
     & /(za(j2,j3)*zb(j3,j4)*za(j5,j6))

     & +za(j2,j4)**3*zb(j2,j6)**2*zab2(j1,j2,j3,j4)
     & /(za(j2,j3)*za(j3,j4)*zb(j5,j6))
     & -za(j2,j4)*za(j4,j5)
     & *zab2(j1,j2,j3,j4)*(zab(j4,j1,j6)-zab(j4,j2,j6))
     & /(za(j2,j3)*za(j3,j4))
     & +za(j1,j2)*zb(j3,j6)*(zab(j5,j1,j3)-zab(j5,j2,j3))
     & *(s(j2,j3)+t(j1,j3,j4))/(za(j2,j3)*zb(j3,j4))

     & -za(j1,j2)*za(j4,j5)/za(j2,j3)
     & *(zb(j3,j6)*(2*s(j2,j4)+delta34)-2*zba(j3,j1,j5)*zb(j5,j6)))


     & +0.5_dp/(zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4))
     & *(zab(j1,j4,j6)
     & *(zab(j5,j1,j3)*s(j1,j3)-zab(j5,j2,j3)*s(j2,j3))
     & +5*zab(j1,j2,j4)*zab(j4,j1,j3)*zab(j5,j4,j6)

     & +zab(j5,j4,j3)*zab2(j1,j2,j3,j6)*(s(j1,j3)-s(j2,j4))
     & +zab(j1,j2,j3)
     & *(zab(j5,j2,j6)*s(j3,j4)
     &  +zab(j5,j4,j6)*s(j1,j3)
     & +2*zab(j2,j1,j4)*za(j4,j3)*zab2(j5,j1,j3,j6)/za(j2,j3)

     & +3*(za(j5,j4)*zb(j4,j3)*za(j3,j1)*zb(j1,j6)
     & -zab(j5,j3,j6)*s(j2,j4))
     & -2*(zab(j5,j1,j6)*(s(j1,j4)+s(j2,j4))
     & +za(j5,j1)*zb(j1,j4)*za(j4,j3)*zb(j3,j6))))

      return
      end
