!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFPMcc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.16
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp)::FFPMcc,FFPMcc_unsym
      integer:: j1,j2,j3,j4,j5,j6
      FFPMcc=FFPMcc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
     &      +FFPMcc_unsym(j2,j1,j4,j3,j6,j5,zb,za)
      return
      end

      function FFPMcc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: FFPMcc_unsym,zab2,FFPMccTtilde,FFPMccT
      complex(dp):: Master1,L0,Lsm1_2mh,I3m,Lsm1
      complex(dp):: FFPMcc1,FFPMcc2,FFPMcc3
      real(dp):: t,s123,s234
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

c      FFPMcc1=
c     & +2*zab(j5,j1,j3)*zab2(j1,j2,j4,j3)*zab2(j5,j2,j4,j3)
c     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2))
c      FFPMcc2=
c     & +2*zab(j4,j2,j3)*zab2(j5,j2,j3,j4)*zab2(j5,j2,j4,j3)
c     & /(za(j5,j6)*t(j2,j3,j4)*zb(j2,j4)*zab2(j1,j2,j3,j4))
c      FFPMcc3=
c     & +2*zab(j4,j2,j3)*zab2(j5,j2,j4,j3)*zab2(j5,j3,j4,j2)
c     & /(za(j5,j6)*t(j2,j3,j4)*zb(j2,j4)*zab2(j1,j3,j4,j2))
     
      s123=t(j1,j2,j3)
      s234=t(j2,j3,j4)

      FFPMcc_unsym=((za(j1,j2)*zb(j4,j6)*s123)**2
     & /(za(j2,j3)*zb(j5,j6)*zab2(j1,j2,j3,j4)**3*zab2(j3,j1,j2,j4))
     & -(zab2(j1,j2,j3,j6)*zab2(j2,j1,j3,j4))**2
     & /(za(j2,j3)*zb(j5,j6)*zab2(j1,j2,j3,j4)**3*zab2(j3,j1,j2,j4)))
     & *Lsm1_2mh(s(j1,j4),s123,s(j2,j3),s(j5,j6))

     & +((zab2(j2,j1,j3,j4)*zab2(j3,j1,j2,j6))**2
     & /(za(j2,j3)*zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)**3)
     & -za(j2,j3)*(zb(j4,j6)*s123)**2
     & /(zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)**3))
     & *Lsm1_2mh(s(j3,j4),s123,s(j1,j2),s(j5,j6))

     & +((zb(j2,j3)*za(j1,j5))**2*s234*zab2(j1,j2,j4,j3)
     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2)**3)
     & -zab2(j1,j2,j4,j3)**3*zab2(j5,j3,j4,j2)**2
     & /(zb(j3,j4)*za(j5,j6)*s234
     & *zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2)**3)
     & -za(j2,j4)**3*zb(j1,j6)**2
     & /(za(j2,j3)*za(j3,j4)*zb(j5,j6)*s234*zab2(j2,j3,j4,j1)))
     & *Lsm1_2mh(s(j1,j2),s234,s(j3,j4),s(j5,j6))

     & +((za(j1,j2)*zab2(j3,j1,j2,j6))**2
     & /(za(j2,j3)*zb(j5,j6)*za(j1,j3)**2
     & *zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4))
     & -za(j2,j3)*zab2(j1,j2,j3,j6)**2
     & /(zb(j5,j6)*za(j1,j3)**2*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)))
     & *Lsm1(-s(j1,j2),-s123,-s(j2,j3),-s123)

     & +((zab2(j5,j2,j3,j4)*zb(j2,j3))**2
     & /(za(j5,j6)*s234*zb(j2,j4)**3*zab2(j1,j2,j3,j4))
     & -(zab2(j5,j3,j4,j2)*zb(j3,j4))**2
     & /(za(j5,j6)*s234*zb(j2,j4)**3*zab2(j1,j2,j3,j4)))
     & *Lsm1(-s(j2,j3),-s234,-s(j3,j4),-s234)

     & +FFPMccT(j1,j2,j3,j4,j5,j6,za,zb)
     & *I3m(s(j1,j2),s(j3,j4),s(j5,j6))
     & +FFPMccTtilde(j1,j2,j3,j4,j5,j6,za,zb)
     & *I3m(s(j1,j4),s(j2,j3),s(j5,j6))


     & +FFPMcc1(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j5,j6))/s(j5,j6)
     & +FFPMcc2(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j2,j3))/s(j2,j3)
     & +FFPMcc3(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j3,j4))/s(j3,j4)

     & + Master1(j1,j4,j3,j2,j5,j6,za,zb)
c     +flip2
      return
      end


      function FFPMcc1(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: FFPMcc1,zab2,zab
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPMcc1=
     & +2*zab(j5,j1,j3)*zab2(j1,j2,j4,j3)*zab2(j5,j2,j4,j3)
     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2))

      return
      end

      function FFPMcc2(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: FFPMcc2,zab2,zab
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPMcc2=
     & +2*zab(j4,j2,j3)*zab2(j5,j2,j3,j4)*zab2(j5,j2,j4,j3)
     & /(za(j5,j6)*t(j2,j3,j4)*zb(j2,j4)*zab2(j1,j2,j3,j4))

      return
      end

      function FFPMcc3(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: FFPMcc3,zab2,zab
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPMcc3=
     & +2*zab(j4,j2,j3)*zab2(j5,j2,j4,j3)*zab2(j5,j3,j4,j2)
     & /(za(j5,j6)*t(j2,j3,j4)*zb(j2,j4)*zab2(j1,j3,j4,j2))

      return
      end

