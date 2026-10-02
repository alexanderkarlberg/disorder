!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPPFsc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq8.8
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FPPFsc,Atree,Atrees,zab2,zba2,zab,zb12b
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L1
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zb12b(j1,j2,j3,j4,j5)=zb(j1,j2)*zab2(j2,j3,j4,j5)
c     End statement functions

      Atree=-za(j4,j5)**2/(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j5,j6))
      Atrees=1._dp/(3._dp*za(j2,j3)**2*s(j5,j6))
     & *(-za(j4,j5)*zba2(j6,j1,j2,j3)*zb(j3,j1)/t(j1,j2,j3)
     &   +zb(j1,j6)*zab2(j5,j4,j2,j3)*za(j3,j4)/t(j2,j3,j4))
      FPPFsc=Atree
     & *(-0.5_dp*(zab(j4,j3,j2)*za(j2,j5)/za(j4,j5))**2
     & *L1(-s(j3,j4),-t(j2,j3,j4))/t(j2,j3,j4)**2
     & +0.5_dp*(zab2(j4,j2,j3,j1)*za(j1,j5)/za(j4,j5))**2
     & *L1(-s(j5,j6),-t(j2,j3,j4))/t(j2,j3,j4)**2)
     & +0.5_dp*(
     & -zab(j5,j2,j3)*za(j5,j4)
     & /(za(j1,j2)*za(j2,j3)*t(j2,j3,j4)*za(j5,j6))
     & -zb(j2,j3)*za(j4,j5)*zab2(j5,j2,j4,j3)
     & /(za(j1,j2)*t(j1,j2,j3)*t(j2,j3,j4)*za(j5,j6))

     & +zab2(j4,j2,j3,j6)**2
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*t(j2,j3,j4)*zb(j5,j6))
     & +zab(j4,j2,j3)*zb12b(j6,j1,j2,j3,j6)
     & /(za(j1,j2)*za(j2,j3)*t(j1,j2,j3)*t(j2,j3,j4)*zb(j5,j6)))

     & +Atrees
      return
      end
