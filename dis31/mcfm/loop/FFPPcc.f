!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFPPcc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.4
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FFPPcc
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,Lsm1_2me,Lsm1,FFPPcc1,FFPPcc2
      real(dp):: t,s123,s234

c      FFPPcc1=
c     & -2*za(j5,j2)*zab(j5,j4,j3)
c     & /(za(j3,j4)*za(j4,j1)*za(j5,j6))
c      FFPPcc2=
c     & +2*za(j2,j5)*za(j1,j5)*zab2(j2,j3,j4,j1)
c     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))

      s123=t(j1,j2,j3)
      s234=t(j2,j3,j4)
      FFPPcc=
     & za(j2,j5)*(za(j1,j2)*za(j4,j5)-za(j2,j4)*za(j1,j5))
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)**2*za(j5,j6))
     & *Lsm1_2me(s123,s234,s(j2,j3),s(j5,j6))

     & +za(j2,j5)*(za(j2,j3)*za(j1,j5)-za(j1,j2)*za(j3,j5))
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j1,j3)*za(j5,j6))
     & *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)

     & -za(j2,j5)*(za(j2,j3)*za(j4,j5)+za(j2,j4)*za(j3,j5))
     & /(za(j2,j3)*za(j3,j4)**2*za(j4,j1)*za(j5,j6))
     & *Lsm1_2me(t(j1,j2,j4),s123,s(j1,j2),s(j5,j6))

     & -za(j5,j2)**2/(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))
     & *(Lsm1(-s(j1,j4),-t(j1,j2,j4),-s(j1,j2),-t(j1,j2,j4))
     & + Lsm1_2me(t(j1,j3,j4),t(j1,j2,j4),s(j1,j4),s(j5,j6)))

     & +FFPPcc1(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j2,j3))/s(j2,j3)
     & +FFPPcc2(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j5,j6))/s(j5,j6)

      return
      end


      function FFPPcc1(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.4
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPcc1,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPPcc1=
     & -2*za(j5,j2)*zab(j5,j4,j3)
     & /(za(j3,j4)*za(j4,j1)*za(j5,j6))

      return
      end


      function FFPPcc2(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.4
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPcc2,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPPcc2=
     & +2*za(j2,j5)*za(j1,j5)*zab2(j2,j3,j4,j1)
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))

      return
      end

