!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPFPcc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.4
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FPFPcc
      complex(dp):: FPFPcc1,FPFPcc2,FPFPcc3
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,Lsm1_2me,Lsm1
      real(dp):: t,s234

c      FPFPcc1=
c     & -2*za(j3,j5)*zab(j5,j4,j2)/(za(j4,j1)*za(j5,j6)*za(j2,j4))
c      FPFPcc2=
c     & -2*za(j3,j5)*zab(j5,j2,j4)/(za(j1,j2)*za(j5,j6)*za(j2,j4))
c      FPFPcc3=
c     & +2*za(j1,j3)*za(j5,j1)*za(j5,j3)*zab2(j3,j2,j4,j1)
c     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))

      s234=t(j2,j3,j4)
      FPFPcc=za(j3,j5)*(za(j2,j3)*za(j4,j5)+za(j3,j4)*za(j5,j2))
     & /(za(j1,j2)*za(j3,j4)*za(j5,j6)*za(j2,j4)**2)
     & *Lsm1(-s(j2,j3),-s234,-s(j3,j4),-s234)

     & -za(j1,j3)*za(j3,j5)**2
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))
     & *Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))

     & +za(j1,j3)*za(j3,j5)*(za(j1,j3)*za(j4,j5)+za(j3,j4)*za(j5,j1))
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1)**2*za(j5,j6))
     & *Lsm1_2me(t(j1,j2,j3),s234,s(j2,j3),s(j5,j6))

     & -za(j1,j3)*za(j3,j5)**2
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))
     & *Lsm1_2me(t(j1,j2,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))

     & -za(j3,j5)**2/(za(j2,j3)*za(j4,j1)*za(j5,j6)*za(j2,j4))
     & *Lsm1(-s(j1,j4),-t(j1,j2,j4),-s(j1,j2),-t(j1,j2,j4))

     & +FPFPcc1(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j2,j3))/s(j2,j3)
     & +FPFPcc2(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j3,j4))/s(j3,j4)
     & +FPFPcc3(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j5,j6))/s(j5,j6)

      return
      end


      function FPFPcc1(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.4
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPcc1,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FPFPcc1=
     & -2*za(j3,j5)*zab(j5,j4,j2)/(za(j4,j1)*za(j5,j6)*za(j2,j4))

      return
      end


      function FPFPcc2(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.4
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPcc2,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FPFPcc2=
     & -2*za(j3,j5)*zab(j5,j2,j4)/(za(j1,j2)*za(j5,j6)*za(j2,j4))

      return
      end


      function FPFPcc3(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.4
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPcc3,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FPFPcc3=
     & +2*za(j1,j3)*za(j5,j1)*za(j5,j3)*zab2(j3,j2,j4,j1)
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))

      return
      end



