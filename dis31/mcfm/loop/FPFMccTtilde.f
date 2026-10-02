!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPFMccTtilde(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq 9.12
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FPFMccTtilde,zab2,zba2,zab
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t,deltat23
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      deltat23=s(j2,j3)-s(j5,j6)-s(j1,j4)
      FPFMccTtilde=
     & 0.5_dp*(t(j1,j2,j4)*deltat23+2*s(j1,j4)*s(j5,j6))/t(j1,j2,j4)**2
     & *(zab2(j4,j1,j2,j6)**2/(za(j1,j2)*zb(j5,j6)*zab2(j2,j1,j4,j3))
     &+(zb(j1,j2)*za(j3,j5))**2/(zb(j1,j4)*za(j5,j6)*zab2(j3,j1,j2,j4)))

     & -2*zab(j4,j1,j2)*za(j2,j4)*zab(j5,j3,j6)
     & /(za(j1,j2)*t(j1,j2,j4)*zab2(j2,j1,j4,j3))
     & -zab(j3,j1,j2)*za(j4,j5)*zab2(j3,j1,j2,j6)
     & /(za(j1,j2)*t(j1,j2,j3)*zab2(j3,j1,j2,j4))

     & +1._dp/(zab2(j2,j1,j4,j3)*zab2(j3,j1,j2,j4))
     & *(zab(j4,j1,j2)*za(j2,j3)*zab(j5,j3,j6)/za(j1,j2)
     & +zab(j3,j1,j2)*zab(j5,j4,j6)/(za(j1,j2)*zab2(j1,j2,j3,j4))
     & *(za(j1,j2)*(s(j1,j3)-s(j2,j4))+za(j1,j3)*zba2(j3,j1,j4,j2))

     & -0.5_dp*zab2(j4,j1,j2,j6)/(za(j1,j2)*zb(j5,j6))
     & *(zab2(j3,j1,j2,j6)*deltat23+2*zab(j3,j1,j6)*s(j5,j6))
     & -0.5_dp*zb(j1,j2)*za(j3,j5)*zab2(j5,j3,j4,j1)*deltat23
     & /(zb(j1,j4)*za(j5,j6))

     & +zb(j2,j6)*(zab(j4,j2,j1)*za(j3,j5)+zab(j3,j4,j1)*za(j4,j5))
     & -za(j3,j4)*zb(j1,j6)*zab2(j5,j1,j3,j2))

      return
      end
