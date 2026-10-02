!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPMFsc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq8.14
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp)::FPMFsc,FPMFsc_unsym
      integer:: j1,j2,j3,j4,j5,j6
      FPMFsc=FPMFsc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
     &      +FPMFsc_unsym(j4,j3,j2,j1,j6,j5,zb,za)
      return
      end

      function FPMFsc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FPMFsc_unsym,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L1
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FPMFsc_unsym=-0.5_dp*za(j1,j3)
     & /(za(j1,j2)*za(j2,j3)*zb(j5,j6)*t(j1,j2,j3)*zab2(j1,j2,j3,j4))
     & *((za(j3,j2)*zb(j2,j1)*zab2(j1,j2,j3,j6))**2
     & *L1(-t(j1,j2,j3),-s(j2,j3))/s(j2,j3)**2
     & +zab(j3,j4,j6)**2*L1(-s(j5,j6),-t(j1,j2,j3)))
     & +0.5_dp*zb(j6,j2)**2/(za(j1,j2)*zb(j2,j3)*zb(j3,j4)*zb(j5,j6))
c     +flip1
      return
      end
