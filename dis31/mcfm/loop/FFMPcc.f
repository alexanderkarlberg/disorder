!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFMPcc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.10
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp)::FFMPcc,FFMPcc_unsym
      integer:: j1,j2,j3,j4,j5,j6
      FFMPcc=FFMPcc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
     &      +FFMPcc_unsym(j2,j1,j4,j3,j6,j5,zb,za)
      return
      end

      function FFMPcc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.10
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FFMPcc_unsym,zab2,zba2,zab,za22a,zb21b,za21a
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,Lsm1_2mh,Lsm1_2mht,I3m,Lsm1,FFMPcc1,FFMPcc2
      real(dp):: t,delta12
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      za22a(j1,j2,j3,j4,j5,j6)=zab2(j1,j2,j3,j4)*za(j4,j6)
     &                        +zab2(j1,j2,j3,j5)*za(j5,j6)
      zb21b(j1,j2,j3,j4,j5)=zba2(j1,j2,j3,j4)*zb(j4,j5)
      za21a(j1,j2,j3,j4,j5)=zab2(j1,j2,j3,j4)*za(j4,j5)
c     End statement functions
      delta12=s(j1,j2)-s(j3,j4)-s(j5,j6)

c      FFMPcc1=
c     & +2*za(j2,j3)*zb(j2,j4)*zab2(j5,j2,j3,j4)*zab2(j5,j3,j4,j2)
c     & /(zb(j2,j3)*za(j5,j6)*t(j2,j3,j4)*zab2(j1,j3,j4,j2))
c      FFMPcc2=
c     & -2*zb(j1,j4)*zb(j2,j4)*za(j1,j5)*zab2(j5,j2,j3,j4)
c     & /(zb(j2,j3)*zb(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2))

      FFMPcc_unsym=zab2(j5,j2,j3,j1)**2
     & /(zb(j2,j3)*za(j5,j6)*zab2(j4,j1,j2,j3)*zab2(j4,j2,j3,j1))
     & *(Lsm1_2mht(s(j1,j4),t(j1,j2,j3),s(j2,j3),s(j5,j6))
     & +Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3)))

     & +(zb(j2,j4)**3*za(j1,j5)**2*t(j2,j3,j4)
     & /(zb(j2,j3)*zb(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2)**3)
     & -zb(j2,j4)*(zab2(j1,j2,j3,j4)*zab2(j5,j3,j4,j2))**2
     &/(zb(j2,j3)*zb(j3,j4)*za(j5,j6)*t(j2,j3,j4)*zab2(j1,j3,j4,j2)**3)

     & -(za(j2,j3)*zb(j1,j6))**2*zab2(j3,j2,j4,j1)
     & /(za(j3,j4)*zb(j5,j6)*t(j2,j3,j4)
     & *zab2(j2,j3,j4,j1)*zab2(j4,j2,j3,j1)))
     & *Lsm1_2mh(s(j1,j2),t(j2,j3,j4),s(j3,j4),s(j5,j6))

     & +((zb(j1,j3)*za(j4,j5)*t(j1,j2,j3))**2
     & /(zb(j2,j3)*za(j5,j6)*zab2(j4,j1,j2,j3)**3*zab2(j4,j2,j3,j1))
     & -zab2(j4,j2,j3,j1)*zab2(j5,j1,j2,j3)**2
     & /(zb(j2,j3)*za(j5,j6)*zab2(j4,j1,j2,j3)**3))
     & *Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))

     & -(za(j2,j3)*zb(j1,j6))**2
     & /(zb(j5,j6)*t(j2,j3,j4)*za(j2,j4)*zab2(j4,j2,j3,j1))
     & *Lsm1(-s(j2,j3),-t(j2,j3,j4),-s(j3,j4),-t(j2,j3,j4))

     & -0.5_dp*zb(j1,j4)*(za22a(j2,j1,j4,j2,j3,j5)**2
     &  -za(j2,j5)**2*s(j1,j4)*s(j2,j3))
     & /(za(j1,j4)*zb(j2,j3)*za(j5,j6)
     & *zab2(j2,j1,j4,j3)*zab2(j2,j3,j4,j1))
     & *I3m(s(j1,j4),s(j2,j3),s(j5,j6))

     & +(0.5_dp*zb(j2,j4)*zab2(j5,j2,j3,j4)**2
     & *(2*s(j3,j4)*s(j5,j6)+delta12*t(j2,j3,j4))
     & /(zb(j2,j3)*zb(j3,j4)*za(j5,j6)*t(j2,j3,j4)**2*zab2(j1,j3,j4,j2))

     & +za(j2,j3)*zb(j1,j6)*zab2(j3,j2,j4,j1)*zab2(j5,j2,j3,j4)
     & /(t(j2,j3,j4)**2*zab2(j4,j2,j3,j1))

     & +0.5_dp*za(j2,j3)*zb(j1,j6)*zab2(j3,j1,j4,j6)*zab2(j3,j2,j4,j1)
     & /(za(j3,j4)*zb(j5,j6)*t(j2,j3,j4)*zab2(j4,j2,j3,j1))
     & +0.5_dp*zb(j1,j4)*za(j2,j5)*zb(j1,j6)*zab(j3,j2,j4)
     & /(zb(j3,j4)*t(j2,j3,j4)*zab2(j4,j2,j3,j1))
     & +2._dp*za(j2,j3)*zb(j2,j4)**2*zab(j5,j1,j6)
     & /(zb(j2,j3)*t(j2,j3,j4)*zab2(j1,j3,j4,j2))

     & -zab(j3,j2,j1)*zab2(j5,j2,j3,j1)*zb(j4,j6)
     & /(zb(j2,j3)*t(j1,j2,j3)*zab2(j4,j2,j3,j1))

     & +0.5_dp/(zab2(j1,j3,j4,j2)*zab2(j4,j2,j3,j1))

     & *(zb(j1,j4)*zab(j3,j1,j6)*zab(j5,j2,j4)/zb(j3,j4)
     & +zb(j1,j4)*zab(j5,j4,j2)*zab2(j5,j2,j3,j4)*delta12
     & /(zb(j2,j3)*zb(j3,j4)*za(j5,j6))

     & -za(j2,j3)*zb(j1,j6)/(za(j3,j4)*zb(j5,j6))
     & *(zab(j3,j2,j6)*zab2(j3,j1,j4,j2)+zab(j3,j1,j2)*zab(j3,j5,j6))

     & +2*zb(j1,j2)*zab(j3,j5,j6)*zab2(j5,j2,j3,j4)/zb(j2,j3)

     & +zb(j2,j4)*zb(j1,j6)/(zb(j2,j3)*zb(j3,j4))
     & *(zab(j5,j2,j4)*(s(j5,j6)-t(j1,j2,j3))
     &  +zab(j5,j3,j4)*(s(j5,j6)+t(j2,j3,j4)))

     & +zb(j1,j4)*za(j3,j5)*(2*zab(j3,j4,j6)+zab(j3,j2,j6)))

     & +0.5_dp/(zab2(j1,j3,j4,j2)*zab2(j4,j1,j2,j3)*zab2(j4,j2,j3,j1))
     & *(za(j4,j1)*zb(j1,j2)*za(j2,j3)*zb(j3,j4)/zb(j2,j3)
     & *((zab(j5,j1,j2)-zab(j5,j4,j2))*zb(j1,j6)
     & -2*zb(j1,j2)*zab(j5,j4,j6))

     & +zb21b(j1,j3,j4,j2,j4)
     & *(za(j3,j5)*zab(j4,j1,j6)+za(j4,j5)*zab(j3,j2,j6))
     & +zab(j5,j3,j6)*(zab(j3,j2,j1)*s(j1,j3)+zab(j3,j4,j1)*s(j1,j2))

     & +zb(j1,j3)*(zab(j3,j1,j6)-zab(j3,j2,j6))*za21a(j3,j1,j2,j4,j5)
     & +zab2(j5,j3,j4,j1)*zab(j3,j4,j6)*s(j2,j3)
     & +zab(j5,j4,j1)*zab2(j3,j1,j2,j6)*s(j1,j3)

     & +(3*zab(j3,j2,j1)+zab(j3,j4,j1))*zab(j5,j2,j6)*s(j1,j4)
     & -zab(j5,j3,j1)*zab2(j1,j3,j4,j6)*(zab(j3,j2,j1)-zab(j3,j4,j1))

     & +2*zab(j3,j4,j1)*zb(j1,j6)
     & *(za(j5,j1)*(s(j1,j2)+t(j1,j2,j3))-zab(j5,j2,j3)*za(j3,j1))
     & -zab(j3,j2,j1)*zab(j5,j4,j6)*s(j3,j4)

     & +4*za(j5,j1)*zb(j1,j2)*za(j2,j3)*zb(j3,j1)*zab(j3,j4,j6)

     & -za(j2,j3)*zab2(j4,j2,j3,j1)
     & *(2._dp*(zab(j5,j2,j4)-zab(j5,j3,j4))*zb(j2,j6)
     & -zab(j5,j1,j6)*zb(j2,j4))))
     &  *I3m(s(j1,j2), s(j3,j4), s(j5,j6))

     & +FFMPcc1(j1,j2,j3,j4,j5,j6,za,zb)
     & *L0(-t(j2,j3,j4),-s(j3,j4))/s(j3,j4)
     & +FFMPcc2(j1,j2,j3,j4,j5,j6,za,zb)
     & *L0(-t(j2,j3,j4),-s(j5,j6))/s(j5,j6)

      return
      end

      function FFMPcc1(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFMPcc1,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFMPcc1=
     & +2*za(j2,j3)*zb(j2,j4)*zab2(j5,j2,j3,j4)*zab2(j5,j3,j4,j2)
     & /(zb(j2,j3)*za(j5,j6)*t(j2,j3,j4)*zab2(j1,j3,j4,j2))

      return
      end

      function FFMPcc2(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFMPcc2,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFMPcc2=
     & -2*zb(j1,j4)*zb(j2,j4)*za(j1,j5)*zab2(j5,j2,j3,j4)
     & /(zb(j2,j3)*zb(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2))

      return
      end

