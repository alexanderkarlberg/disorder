!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPMFcc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq8.12
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp)::FPMFcc,FPMFcc_unsym
      integer:: j1,j2,j3,j4,j5,j6
      FPMFcc=FPMFcc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
     &      +FPMFcc_unsym(j4,j3,j2,j1,j6,j5,zb,za)
      return
      end

      function FPMFcc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FPMFcc_unsym,zab2,zba2,zab,zb21b
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,Lsm1_2mht,I3m,Lsm1
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zb21b(j1,j2,j3,j4,j5)=zba2(j1,j2,j3,j4)*zb(j4,j5)
c     End statement functions

      FPMFcc_unsym=(za(j1,j3)*zab2(j3,j1,j2,j6)**2
     & /(za(j1,j2)*za(j2,j3)*zb(j5,j6)*t(j1,j2,j3)*zab2(j1,j2,j3,j4))
     & +zb(j1,j2)**3*za(j4,j5)**2
     & /(zb(j2,j3)*zb(j1,j3)*za(j5,j6)*t(j1,j2,j3)*zab2(j4,j2,j3,j1)))
     & *Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))

     & +(za(j1,j3)*zab2(j3,j1,j2,j6)**2
     & /(za(j1,j2)*za(j2,j3)*zb(j5,j6)*t(j1,j2,j3)*zab2(j1,j2,j3,j4))
     & +(zb(j1,j2)*za(j4,j5))**2*zab2(j4,j1,j3,j2)
     & /(zb(j2,j3)*za(j5,j6)*t(j1,j2,j3)
     & *zab2(j4,j2,j3,j1)*zab2(j4,j1,j2,j3)))
     & *Lsm1_2mht(s(j3,j4),t(j1,j2,j3),s(j5,j6),s(j1,j2))

     & +0.5_dp*zb(j1,j2)
     & *((zab2(j4,j1,j2,j3)*za(j3,j5)+zab2(j4,j1,j2,j4)*za(j4,j5))**2
     & -s(j1,j2)*s(j3,j4)*za(j4,j5)**2)
     & /(za(j1,j2)*zb(j3,j4)*za(j5,j6)*zab2(j4,j2,j3,j1)
     & *zab2(j4,j1,j2,j3))*I3m(s(j1,j2),s(j3,j4),s(j5,j6))

     & -2*za(j1,j3)*zab2(j3,j1,j2,j6)
     & /(za(j1,j2)*zb(j5,j6)*zab2(j1,j2,j3,j4))
     & *(zb21b(j6,j2,j3,j1,j2)/t(j1,j2,j3)
     & *L0(-s(j2,j3),-t(j1,j2,j3))/t(j1,j2,j3)
     & +zab(j3,j4,j6)/za(j2,j3)
     & *L0(-s(j5,j6),-t(j1,j2,j3))/t(j1,j2,j3))

c      +flip1
      return
      end
