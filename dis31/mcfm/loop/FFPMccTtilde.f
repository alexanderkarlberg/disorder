!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFPMccTtilde(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.10.18
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMccTtilde,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPMccTtilde=1._dp/(zab2(j1,j2,j3,j4)*zab2(j3,j1,j4,j2))
     & *(0.5_dp*(zab(j4,j1,j3)-zab(j4,j2,j3))*zab2(j5,j1,j4,j6)
     & +zab(j5,j1,j3)*zab(j4,j2,j6)

     & +zab2(j4,j1,j2,j3)/zab2(j3,j1,j2,j4)
     & *(zab(j5,j2,j4)*zab(j3,j1,j6)
     & -zab(j3,j1,j4)*(zab(j5,j2,j6)-zab(j5,j3,j6))))
     & -zb(j1,j3)*za(j4,j5)*zab2(j2,j1,j3,j6)
     & /(zab2(j3,j1,j2,j4)*t(j1,j2,j3))

      return
      end
