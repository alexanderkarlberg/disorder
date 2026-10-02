!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPPFcc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq8.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FPPFcc,Atree,zba2,zab,za12a
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,Lsm1,Lsm1_2me
      real(dp):: t
c     Statement functions
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      za12a(j1,j2,j3,j4,j5)=za(j1,j2)*zba2(j2,j3,j4,j5)
c     End statement functions

      Atree=-za(j4,j5)**2/(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j5,j6))
      FPPFcc=Atree
     & *(-Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))
     &   -Lsm1(-s(j2,j3),-t(j2,j3,j4),-s(j3,j4),-t(j2,j3,j4))
     &   -Lsm1_2me(t(j1,j2,j3),t(j2,j3,j4),s(j2,j3),s(j5,j6)))

     & +2*za(j4,j5)*zab(j5,j2,j3)/(za(j1,j2)*za(j2,j3)*za(j5,j6))
     & *L0(-t(j2,j3,j4),-s(j3,j4))/s(j3,j4)
     & +2*za(j4,j5)*za12a(j5,j1,j2,j3,j4)
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j5,j6))
     & *L0(-s(j5,j6),-t(j2,j3,j4))/t(j2,j3,j4)

      return
      end
