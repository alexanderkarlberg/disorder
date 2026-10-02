!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FMPFsc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq8.20
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp)::FMPFsc,FMPFsc_unsym
      integer:: j1,j2,j3,j4,j5,j6
      FMPFsc=FMPFsc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
     &      +FMPFsc_unsym(j4,j3,j2,j1,j6,j5,zb,za)
      return
      end

      function FMPFsc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FMPFsc_unsym,zab2,zab,zba,zb12b
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lsm1_2mh,I3m,Lnrat,Ls1
      real(dp):: t,delta12,delta34,delta56,Delta3
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
      zb12b(j1,j2,j3,j4,j5)=zb(j1,j2)*zab2(j2,j3,j4,j5)
c     End statement functions
      delta12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      delta34=s(j3,j4)-s(j5,j6)-s(j1,j2)
      delta56=s(j5,j6)-s(j1,j2)-s(j3,j4)
      Delta3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2*s(j1,j2)*s(j3,j4)-2*s(j3,j4)*s(j5,j6)-2*s(j5,j6)*s(j1,j2)

      FMPFsc_unsym=
     &  za(j1,j2)*za(j2,j3)*zb(j1,j3)**2*zab2(j1,j2,j3,j6)**2
     & /(zb(j5,j6)*za(j1,j3)*t(j1,j2,j3)*zab2(j1,j2,j3,j4))
     & *(Ls1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))
     & /t(j1,j2,j3)**2
     & -0.5_dp*L1(-t(j1,j2,j3),-s(j2,j3))/s(j2,j3)**2)
c       pmpm_scalar_b2
     & +za(j2,j3)*zb(j4,j6)**2*t(j1,j2,j3)*zab2(j2,j1,j3,j4)
     & /(zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)**3)
     & *Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))
c       Scalar three-mass triangle (D=6 form) pmpm_scalar_3m_1
     & +za(j1,j2)*zb(j4,j6)/(zab2(j1,j2,j3,j4)*Delta3)
     & *(3*zb(j3,j4)*za(j5,j6)*zab2(j2,j3,j4,j1)
     & *(zab(j3,j5,j6)*delta56-zab(j3,j4,j6)*delta34)
     & *zab2(j4,j1,j2,j3)/(zab2(j3,j1,j2,j4)*Delta3)
     & +(3*zab(j5,j6,j4)*zab(j2,j3,j1)
     & -zab(j5,j3,j4)*zab2(j2,j3,j4,j1))
     & *zab2(j4,j1,j2,j3)/zab2(j3,j1,j2,j4)
     & -zb(j1,j3)*za(j2,j4)*zb(j3,j6)*za(j5,j6)
     & +zb(j1,j4)*za(j2,j3)*zab(j5,j6,j4)*(t(j1,j2,j3)-t(j1,j2,j4))
     & *zab2(j4,j1,j2,j3)/zab2(j3,j1,j2,j4)**2)
     & *I3m(s(j1,j2),s(j3,j4),s(j5,j6))
c       L_0 terms: pmpm_scalar_lzero_1
     & +za(j2,j4)*zb(j4,j6)**2*zab2(j2,j1,j3,j4)*t(j1,j2,j3)
     & /(zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4))
     & *(-0.5_dp*za(j2,j4)/za(j2,j3)*L1(-s(j5,j6),-t(j1,j2,j3))
     & /t(j1,j2,j3)**2
     & +1._dp/zab2(j3,j1,j2,j4)
     & *L0(-s(j5,j6),-t(j1,j2,j3))/t(j1,j2,j3))
     & +zab(j2,j1,j3)*za(j3,j2)*zb(j4,j6)**2*t(j1,j2,j3)
     &  /(zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)**2)
     & *L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)
     & +0.5_dp*(zab(j2,j1,j3)*zab2(j3,j1,j2,j6))**2
     & /(zb(j5,j6)*za(j1,j3)*t(j1,j2,j3)*zab2(j3,j1,j2,j4))
     & *L1(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)**2
c       pmpm_scalar_lns12_1
     & +za(j1,j2)*zb(j4,j6)/zab2(j1,j2,j3,j4)*(
     &  3*(zab(j5,j3,j4)*delta34-zab(j5,j6,j4)*delta56)
     & *zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3)
     & /(zab2(j3,j1,j2,j4)*Delta3**2)
     & +0.5_dp/(za(j3,j4)*zb(j5,j6)*zab2(j3,j1,j2,j4)*Delta3)
     & *(-za(j2,j4)*delta12
     & *(zba(j6,j5,j3)*zb(j3,j1)+zb12b(j6,j4,j2,j3,j1))
     & +zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3)
     & *(zab(j3,j4,j6)-zab(j3,j5,j6))
     & -zab2(j2,j3,j4,j1)*zab(j3,j4,j6)*delta34
     & *(t(j1,j2,j3)-t(j1,j2,j4))/zab2(j3,j1,j2,j4))
     & +0.5_dp*zb(j4,j6)*zab2(j2,j3,j4,j1)
     & /(zb(j5,j6)*zab2(j3,j1,j2,j4)**2))*lnrat(-s(j1,j2),-s(j5,j6))
c      pmpm_scalar_lns34_1
     & +zb(j4,j6)/(zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4))
     & *(-3*za(j1,j2)*zb(j3,j4)*zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3)
     & *(za(j3,j5)*delta12-2*zab(j3,j4,j6)*za(j6,j5))/Delta3**2
     & -0.5_dp/(zb(j1,j2)*zb(j5,j6)*Delta3)
     & *(zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3)
     & *(zb(j4,j3)*zab2(j3,j1,j2,j6)+zb(j4,j6)*(delta56-2*s(j1,j2)))
     & -delta34*zab(j2,j4,j3)
     & *(zba(j6,j5,j3)*zb(j3,j1)+zb12b(j6,j4,j2,j3,j1))
     & -2*zb(j4,j6)*t(j1,j2,j3)*(2*(zab(j2,j3,j1)-zab(j2,j4,j1))
     & *zab2(j4,j1,j2,j3)
     & +(t(j1,j2,j3)-t(j1,j2,j4))*(zb(j1,j3)*za(j2,j4)
     & +zb(j1,j4)*za(j2,j3)*zab2(j4,j1,j2,j3)/zab2(j3,j1,j2,j4))))
     & -0.5_dp*zb(j1,j3)*zab(j2,j4,j6)/(zb(j1,j2)*zb(j5,j6)))
     & *lnrat(-s(j3,j4),-s(j5,j6))
c       pmpm_scalar_poly_1
     & +0.5_dp*zb(j4,j6)*zab2(j2,j3,j4,j1)*zab2(j4,j1,j2,j3)
     & *(zab(j3,j5,j6)*delta56-zab(j3,j4,j6)*delta34)
     & /(zb(j1,j2)*za(j3,j4)*zb(j5,j6)*zab2(j1,j2,j3,j4)
     & *zab2(j3,j1,j2,j4)*Delta3)
     & -0.5_dp*zab(j2,j4,j6)
     & *(zba(j6,j5,j3)*zb(j3,j1)+zb12b(j6,j4,j2,j3,j1))
     & /(zb(j1,j2)*za(j3,j4)*zb(j5,j6)*zab2(j1,j2,j3,j4)
     & *zab2(j3,j1,j2,j4))
     & -0.5_dp*(zb(j1,j3)*zab2(j1,j2,j3,j6))**2
     & /(zb(j1,j2)*zb(j2,j3)*zb(j5,j6)*za(j1,j3)
     & *t(j1,j2,j3)*zab2(j1,j2,j3,j4))
c      + flip1
      return
      end
