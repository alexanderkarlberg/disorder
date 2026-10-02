!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FMPFcc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq8.18
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp)::FMPFcc,FMPFcc_unsym
      integer:: j1,j2,j3,j4,j5,j6
      FMPFcc=FMPFcc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
     &      +FMPFcc_unsym(j4,j3,j2,j1,j6,j5,zb,za)
      return
      end

      function FMPFcc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FMPFcc_unsym,zab2,zab3,zab,Master1,zb12b
      complex(dp):: L0,I3m,Lsm1,Lsm1_2mh
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t,delta34,delta56
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab3(j1,j2,j3,j4,j5)=+za(j1,j2)*zb(j2,j5)
     & +za(j1,j3)*zb(j3,j5)+za(j1,j4)*zb(j4,j5)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zb12b(j1,j2,j3,j4,j5)=zb(j1,j2)*zab2(j2,j3,j4,j5)
c     End statement functions

      delta34=s(j3,j4)-s(j5,j6)-s(j1,j2)
      delta56=s(j5,j6)-s(j1,j2)-s(j3,j4)

      FMPFcc_unsym=(zb(j1,j3)**3*za(j4,j5)**2
     & /(zb(j1,j2)*zb(j2,j3)*za(j5,j6)*t(j1,j2,j3)*zab2(j4,j2,j3,j1))
     & +za(j1,j2)**3*zab2(j3,j1,j2,j6)**2
     & /(za(j2,j3)*zb(j5,j6)*za(j1,j3)**3*t(j1,j2,j3)*zab2(j1,j2,j3,j4))
     & -za(j1,j2)*za(j2,j3)*zab2(j1,j2,j3,j6)**2
     & /(zb(j5,j6)*za(j1,j3)**3*t(j1,j2,j3)*zab2(j1,j2,j3,j4)))
     & *Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))
c       pmpm_gf_b2
     & +(zb(j1,j3)**3*za(j4,j5)**2
     & /(zb(j1,j2)*zb(j2,j3)*za(j5,j6)*t(j1,j2,j3)*zab2(j4,j2,j3,j1))
     & +zab2(j3,j1,j2,j6)**2*zab2(j2,j1,j3,j4)**3
     & /(za(j2,j3)*zb(j5,j6)*t(j1,j2,j3)
     & *zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)**3)
     & -za(j2,j3)*zb(j4,j6)**2*t(j1,j2,j3)*zab2(j2,j1,j3,j4)
     & /(zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)**3))
     & *Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))
c       Three-mass triangle coefficient (D=6 form): pmpm_gf_3m_nodel_1
     & +(-2*zab(j2,j1,j3)*zab(j5,j4,j6)*zab2(j2,j1,j3,j4)
     & /(t(j1,j2,j3)*zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4))
     & -0.5_dp*(t(j1,j2,j3)*delta34+2*s(j1,j2)*s(j5,j6))/t(j1,j2,j3)**2
     & *(zb(j1,j3)**3*za(j4,j5)**2
     & /(zb(j1,j2)*zb(j2,j3)*za(j5,j6)*zab2(j4,j2,j3,j1))
     & +zab2(j2,j1,j3,j6)**2*zab2(j2,j1,j3,j4)
     & /(za(j2,j3)*zb(j5,j6)*zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4)))
     & +0.5_dp/(zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4))
     & *(
     & zb(j3,j6)*zab(j2,j4,j6)*zab2(j2,j1,j3,j4)*zab2(j1,j3,j4,j2)
     & /zb(j5,j6)
     & -zab(j2,j1,j6)*zab(j2,j4,j6)*zab3(j4,j1,j2,j3,j4)
     & *zab2(j1,j3,j4,j2)/(za(j3,j4)*zb(j5,j6))
     & -(zab(j1,j4,j6)*zab(j2,j3,j4)*zab(j5,j6,j1)*za(j1,j4)/za(j3,j4))
     & +za(j1,j5)*zab(j2,j3,j4)
     & *(zb(j3,j6)*t(j2,j3,j4)-zb12b(j3,j4,j1,j3,j6))
     & -zab(j2,j1,j6)*zab(j1,j4,j3)*zab(j5,j2,j4)
     & -0.5_dp*s(j1,j4)*za(j4,j5)*zb(j1,j6)*delta56*zab2(j1,j2,j3,j4)
     & /(zb(j1,j2)*za(j3,j4))
     & +0.5_dp*(s(j1,j4)-s(j2,j3))
     & *za(j1,j2)*zb(j3,j4)*zab2(j5,j2,j3,j6)))
     & *I3m(s(j1,j2),s(j3,j4),s(j5,j6))
c       L_0 terms:
     & -2*zab(j2,j1,j3)*zab2(j2,j1,j3,j6)
     & /(zb(j5,j6)*za(j1,j3)*t(j1,j2,j3))
     & *(zab2(j3,j1,j2,j6)/zab2(j3,j1,j2,j4)
     & *L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)
     & +zab2(j1,j2,j3,j6)/zab2(j1,j2,j3,j4)
     & *L0(-t(j1,j2,j3),-s(j2,j3))/s(j2,j3))
     & -2*zab(j2,j4,j6)*zab2(j2,j1,j3,j6)*zab2(j2,j1,j3,j4)
     & /(za(j2,j3)*zb(j5,j6)*zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4))
     & *L0(-t(j1,j2,j3),-s(j5,j6))/s(j5,j6)
     & +Master1(j1,j2,j3,j4,j5,j6,za,zb)
c      +flip1

      return
      end
