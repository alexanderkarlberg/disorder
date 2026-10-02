!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPFMccT(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.9.11
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMccT,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FPFMccT=
     & 2*za(j3,j4)**2*zab2(j5,j3,j4,j2)*zb(j1,j6)
     & /(za(j2,j3)*t(j2,j3,j4)**2)
     & +za(j1,j3)/(za(j2,j3)*t(j1,j2,j3)*t(j2,j3,j4))
     & *(zab2(j4,j2,j3,j1)*zab2(j5,j3,j4,j2)*zab2(j3,j1,j2,j6)
     & /zab2(j1,j2,j3,j4)
     & -zb(j1,j2)*za(j3,j4)*za(j4,j5)*zb(j1,j6))
      return
      end
