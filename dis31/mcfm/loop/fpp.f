!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function fpp(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: fpp
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lsm1,Lsm1_2mht
      complex(dp):: t0,s1,s2,s3
      real(dp):: t
      s2 = -2*(-za(j1,j3)*zb(j1,j6)-za(j2,j3)*zb(j2,j6))/(za(j1,j2)*zb(j
     &2,j4)+za(j1,j3)*zb(j3,j4))/zb(j5,j6)*(L0(-s(j2,j3),-t(j1,j2,j3))*z
     &b(j1,j2)*(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))/t(j1,j2,j3)**2+
     &L0(-s(j5,j6),-t(j1,j2,j3))/za(j2,j3)*za(j3,j4)*zb(j4,j6)/t(j1,j2,j
     &3))
      s3 = (Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))+Lsm1_2mh
     &t(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6)))*(-za(j4,j5)**2/za(j5,j6
     &)*zb(j1,j2)**2/(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))/zb(j2,j3)
     &/t(j1,j2,j3)+1/za(j2,j3)*(-za(j1,j3)*zb(j1,j6)-za(j2,j3)*zb(j2,j6)
     &)**2/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))/zb(j5,j6)/t(j1,j2,j
     &3))-1/za(j2,j3)/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*(L1(-t(j
     &1,j2,j3),-s(j2,j3))/s(j2,j3)**2*za(j2,j3)**2*zb(j1,j2)**2*(za(j1,j
     &2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))**2+L1(-s(j5,j6),-t(j1,j2,j3))*za
     &(j3,j4)**2*zb(j4,j6)**2)/zb(j5,j6)/t(j1,j2,j3)/2
      s1 = s2+s3
      s2 = s1-2/za(j5,j6)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))/(za(
     &j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*(-L0(-s(j2,j3),-t(j2,j3,j4))
     &*za(j3,j4)*(-za(j2,j5)*zb(j2,j4)-za(j3,j5)*zb(j3,j4))/t(j2,j3,j4)*
     &*2+L0(-s(j5,j6),-t(j2,j3,j4))*za(j1,j5)*zb(j1,j2)/zb(j2,j3)/t(j2,j
     &3,j4))
      t0 = s2-(Lsm1(-s(j3,j4),-t(j2,j3,j4),-s(j2,j3),-t(j2,j3,j4))+Lsm1_
     &2mht(s(j1,j2),t(j2,j3,j4),s(j3,j4),s(j5,j6)))*(1/za(j5,j6)/zb(j2,j
     &3)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))**2/(za(j1,j2)*zb(j2,j
     &4)+za(j1,j3)*zb(j3,j4))/t(j2,j3,j4)-1/za(j2,j3)*za(j3,j4)**2/(za(j
     &2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j1,j6)**2/zb(j5,j6)/t(j2,j
     &3,j4))+1/za(j5,j6)/zb(j2,j3)/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,
     &j4))*(L1(-s(j5,j6),-t(j2,j3,j4))*za(j1,j5)**2*zb(j1,j2)**2+L1(-t(j
     &2,j3,j4),-s(j2,j3))/s(j2,j3)**2*za(j3,j4)**2*zb(j2,j3)**2*(-za(j2,
     &j5)*zb(j2,j4)-za(j3,j5)*zb(j3,j4))**2)/t(j2,j3,j4)/2
      fpp=t0
      end

