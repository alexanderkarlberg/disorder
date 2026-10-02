!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFPMscTtilde(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.10.22
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: Deltat3,deltat14,deltat56
      complex(dp):: FFPMscTtilde,zab2,zab
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions
      deltat14=s(j1,j4)-s(j2,j3)-s(j5,j6)
      deltat56=s(j5,j6)-s(j1,j4)-s(j2,j3)
      Deltat3=s(j1,j4)**2+s(j2,j3)**2+s(j5,j6)**2
     & -2*s(j1,j4)*s(j2,j3)-2*s(j2,j3)*s(j5,j6)-2*s(j5,j6)*s(j1,j4)
      FFPMscTtilde=3*za(j1,j4)*zb(j2,j3)*zab2(j2,j3,j4,j1)
     & /(zab2(j3,j1,j2,j4)*Deltat3**2)
     & *((deltat14*zab(j2,j5,j6)-s(j5,j6)*zab(j2,j3,j6))*zab(j5,j4,j1)
     & -0.5_dp*deltat56*zab(j2,j5,j6)*zab(j5,j6,j1))

     & -zab2(j2,j3,j4,j1)/(zab2(j3,j1,j2,j4)*Deltat3)
     & *(0.5_dp*zab(j4,j5,j6)*(zab(j5,j6,j3)+2*zab(j5,j2,j3))
     & -zab(j5,j2,j3)*zab(j4,j1,j6))

      return
      end
