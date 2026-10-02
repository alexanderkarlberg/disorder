!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFMPsc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.12
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp)::FFMPsc,FFMPsc_unsym
      integer:: j1,j2,j3,j4,j5,j6
      FFMPsc=FFMPsc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
     &      +FFMPsc_unsym(j2,j1,j4,j3,j6,j5,zb,za)
      return
      end function FFMPsc

      function FFMPsc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FFMPsc_unsym,Master2,Master3,zab2,Lnrat,FFMPsc1
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

c      FFMPsc1=
c     & +0.5_dp*za(j4,j5)*zab2(j3,j1,j2,j4)*zab2(j5,j2,j3,j4)
c     & /(s(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2)*zab2(j4,j1,j2,j3))

      FFMPsc_unsym=Master2(j1,j2,j3,j4,j5,j6,za,zb)
     &            +Master3(j1,j2,j3,j4,j5,j6,za,zb)
c      ffmp_scalar_lzero_1
     & +FFMPsc1(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s(j1,j2),-s(j5,j6))
c      ffmp_poly_1
     & -0.5_dp*za(j3,j5)*(za(j2,j3)*za(j4,j5)-za(j2,j5)*za(j3,j4))
     & /(za(j1,j4)*za(j3,j4)*za(j5,j6)*zab2(j4,j1,j2,j3))
     & +0.5_dp*zb(j2,j4)*za(j3,j5)
     & /(zb(j2,j3)*za(j5,j6)*zab2(j1,j3,j4,j2))
     & *(za(j3,j5)/za(j3,j4)+zab2(j5,j2,j3,j4)/t(j2,j3,j4))
c      +flip2
      return
      end function FFMPsc_unsym

      function FFMPsc1(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: zab2,FFMPsc1
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFMPsc1=
     & +0.5_dp*za(j4,j5)*zab2(j3,j1,j2,j4)*zab2(j5,j2,j3,j4)
     & /(s(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2)*zab2(j4,j1,j2,j3))

      return
      end function FFMPsc1
