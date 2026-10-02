!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function fsl(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: fsl

      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lnrat,i3m,Lsm1_2mh
      complex(dp):: t0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,
     & s11,s12,s13,s14,s15
      real(dp):: t
      s4 = Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))*(za(j4,j5)**
     &2/za(j5,j6)/zb(j1,j2)*zb(j1,j3)**2/(-za(j1,j4)*zb(j1,j3)-za(j2,j4)
     &*zb(j2,j3))/t(j1,j2,j3)-1/za(j1,j2)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3
     &)*zb(j2,j4))**3*(-za(j1,j3)*zb(j1,j6)-za(j2,j3)*zb(j2,j6))**2*(-za
     &(j1,j2)*zb(j1,j4)+za(j2,j3)*zb(j3,j4))**2/zb(j5,j6)/t(j1,j2,j3))-3
     &._dp/4._dp*(Lnrat(-s(j3,j4),-s(j5,j6))+Lnrat(-t(j1,j2,j3),-s(j5,j6))
     &)/za(j1,j2)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(-za(j1,j2)
     &*zb(j1,j6)+za(j2,j3)*zb(j3,j6))**2/zb(j5,j6)/t(j1,j2,j3)
      s3 = s4-L0(-t(j1,j2,j3),-s(j5,j6))/s(j5,j6)/za(j1,j2)*za(j2,j3)*za
     &(j2,j4)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))**2*zb(j4,j6)**2
     &/zb(j5,j6)*t(j1,j2,j3)+L1(-s(j5,j6),-t(j1,j2,j3))/za(j1,j2)*za(j2,
     &j4)**2/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*zb(j4,j6)**2/zb(
     &j5,j6)/t(j1,j2,j3)/2
      s4 = s3+(Lnrat(-s(j3,j4),-s(j5,j6))+Lnrat(-t(j1,j2,j3),-s(j5,j6)))
     &/za(j1,j2)*za(j2,j3)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))**2
     &*(-za(j1,j2)*zb(j1,j6)+za(j2,j3)*zb(j3,j6))*zb(j4,j6)/zb(j5,j6)/2
      s2 = s4+2*L0(-t(j1,j2,j3),-s(j5,j6))/s(j5,j6)/za(j1,j2)*za(j2,j4)/
     &(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(-za(j1,j2)*zb(j1,j6)+z
     &a(j2,j3)*zb(j3,j6))*zb(j4,j6)/zb(j5,j6)+3._dp/4._dp*(Lnrat(-s(j3,j4)
     &,-s(j5,j6))+Lnrat(-t(j1,j2,j4),-s(j5,j6)))/za(j5,j6)/zb(j1,j2)*(za
     &(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))**2/(-za(j1,j3)*zb(j1,j4)-za
     &(j2,j3)*zb(j2,j4))/t(j1,j2,j4)-(Lnrat(-s(j3,j4),-s(j5,j6))+Lnrat(-
     &t(j1,j2,j4),-s(j5,j6)))*za(j3,j5)/za(j5,j6)/zb(j1,j2)*zb(j1,j4)*(z
     &a(j2,j5)*zb(j1,j2)+za(j4,j5)*zb(j1,j4))/(-za(j1,j3)*zb(j1,j4)-za(j
     &2,j3)*zb(j2,j4))**2/2
      s3 = s2+L0(-t(j1,j2,j4),-s(j5,j6))/s(j5,j6)*za(j3,j5)**2/za(j5,j6)
     &/zb(j1,j2)*zb(j1,j3)*zb(j1,j4)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(
     &j2,j4))**2*t(j1,j2,j4)-L1(-s(j5,j6),-t(j1,j2,j4))*za(j3,j5)**2/za(
     &j5,j6)/zb(j1,j2)*zb(j1,j3)**2/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j
     &2,j4))/t(j1,j2,j4)/2-2*L0(-t(j1,j2,j4),-s(j5,j6))/s(j5,j6)*za(j3,j
     &5)/za(j5,j6)/zb(j1,j2)*zb(j1,j3)*(za(j2,j5)*zb(j1,j2)+za(j4,j5)*zb
     &(j1,j4))/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))
      s5 = s3
      s8 = i3m(s(j1,j2),s(j3,j4),s(j5,j6))
      s10 = za(j2,j3)*(za(j2,j5)*zb(j1,j2)+za(j3,j5)*zb(j1,j3))/(-za(j1,
     &j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))**2*zb(j4,j6)+(-s(j1,j2)+s(j3,j4
     &)-s(j5,j6))/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2
     &)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j2,j3)*((s(j1,j2)-s
     &(j3,j4)-s(j5,j6))*za(j2,j5)*zb(j1,j2)+(-s(j1,j2)-s(j3,j4)+s(j5,j6)
     &)*za(j5,j6)*zb(j1,j6))/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*
     &*2*zb(j4,j6)-1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1
     &,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j5,j6)*zb(j1,j2)
     &/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(za(j2,j3)*zb(j3,j6)+z
     &a(j2,j4)*zb(j4,j6))**2
      s11 = s10-2*za(j4,j5)*zb(j1,j3)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb
     &(j2,j4))*(-za(j1,j3)*zb(j1,j6)-za(j2,j3)*zb(j2,j6))*(-za(j1,j2)*zb
     &(j1,j4)+za(j2,j3)*zb(j3,j4))/t(j1,j2,j3)**2
      s9 = s11+s(j3,j4)/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s
     &(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j2,j5)*zb(j1,
     &j6)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(t(j1,j2,j3)-t(j1,j
     &2,j4))/2+(-s(j1,j2)**2+2*s(j1,j2)*s(j3,j4)-s(j3,j4)**2+3*s(j3,j4)*
     &(-s(j1,j2)+s(j3,j4)-s(j5,j6))+2*s(j1,j2)*s(j5,j6)+2*s(j3,j4)*s(j5,
     &j6)-s(j5,j6)**2)/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(
     &j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*(-za(j2,j3)*zb
     &(j1,j3)-za(j2,j4)*zb(j1,j4))/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2
     &,j4))*(-za(j3,j5)*zb(j3,j6)-za(j4,j5)*zb(j4,j6))*(t(j1,j2,j3)-t(j1
     &,j2,j4))/2
      s7 = s8*s9
      s9 = -cone
      s11 = Lnrat(-s(j1,j2),-s(j3,j4))
      s13 = -1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s
     &(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)*(-za(j3,j5)*zb(
     &j1,j3)+za(j4,j5)*zb(j1,j4))*zb(j1,j6)/(-za(j1,j3)*zb(j1,j4)-za(j2,
     &j3)*zb(j2,j4))+3._dp/2._dp*(-s(j1,j2)-s(j3,j4)+s(j5,j6))/(s(j1,j2)**
     &2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s
     &(j5,j6)+s(j5,j6)**2)**2*(-za(j2,j3)*zb(j1,j3)-za(j2,j4)*zb(j1,j4))
     &/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(-za(j3,j5)*zb(j3,j6)-
     &za(j4,j5)*zb(j4,j6))*(-t(j1,j2,j3)+t(j1,j2,j4))
      s12 = s13-1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2
     &)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j2,j5)/za(j5,j6)/(-
     &za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(-za(j3,j5)*zb(j1,j3)*t(j
     &1,j2,j3)+za(j4,j5)*zb(j1,j4)*t(j1,j2,j4)+(-s(j1,j2)+s(j3,j4)-s(j5,
     &j6))*za(j3,j5)*zb(j1,j4)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4)
     &)*t(j1,j2,j4))+1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(
     &j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)*za(j3,j
     &5)*zb(j1,j4)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))**2*(-2*za(
     &j2,j5)*zb(j1,j2)*zb(j5,j6)+zb(j1,j6)*(-t(j1,j2,j3)+t(j1,j2,j4)))
      s10 = s11*s12
      s8 = s9*s10
      s6 = s7+s8
      s4 = s5+s6
      s5 = s4+1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*
     &s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)/za(j1,j2)*za(j2,j5)/za(j
     &5,j6)*((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j5)-2*za(j1,j2)*za(j5,j
     &6)*zb(j1,j6))/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(-t(j1,j2
     &,j3)+t(j1,j2,j4))/4
      s6 = s5
      s9 = Lnrat(-s(j1,j2),-s(j3,j4))
      s11 = 1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s(
     &j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j2,j5)*zb(j1,j2)/(-za(j
     &1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(za(j2,j3)*zb(j3,j6)-za(j2,j4
     &)*zb(j4,j6))-1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1
     &,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j2,j3)*zb(j1,j2)
     &/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))**2*zb(j4,j6)*(-2*za(j1
     &,j2)*za(j5,j6)*zb(j1,j6)+za(j2,j5)*(t(j1,j2,j3)-t(j1,j2,j4)))
      s12 = s11
      s14 = 3._dp/2._dp*(-s(j1,j2)-s(j3,j4)+s(j5,j6))/(s(j1,j2)**2-2*s(j1,
     &j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s
     &(j5,j6)**2)**2*(-za(j2,j3)*zb(j1,j3)-za(j2,j4)*zb(j1,j4))/(-za(j1,
     &j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(-za(j3,j5)*zb(j3,j6)-za(j4,j5)
     &*zb(j4,j6))*(t(j1,j2,j3)-t(j1,j2,j4))
      s15 = 1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s(
     &j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j1,j6)/(-za(j1,j3)*zb(j
     &1,j4)-za(j2,j3)*zb(j2,j4))/zb(j5,j6)*(za(j2,j3)*zb(j3,j6)*t(j1,j2,
     &j3)+(-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j3)/(-za(j1,j3)*zb(j1,j4)-
     &za(j2,j3)*zb(j2,j4))*zb(j4,j6)*t(j1,j2,j3)-za(j2,j4)*zb(j4,j6)*t(j
     &1,j2,j4))
      s13 = s14+s15
      s10 = s12+s13
      s8 = s9*s10
      s9 = -1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s(
     &j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)/zb(j1,j2)*zb(j1,j6)/(-za(j
     &1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))/zb(j5,j6)*((-s(j1,j2)+s(j3,j4
     &)-s(j5,j6))*zb(j1,j6)-2*za(j2,j5)*zb(j1,j2)*zb(j5,j6))*(t(j1,j2,j3
     &)-t(j1,j2,j4))/4
      s7 = s8+s9
      s1 = s6+s7
      s3 = s1-L0(-t(j3,j5,j6),-s(j1,j2))/s(j1,j2)*za(j3,j5)*za(j4,j5)/za
     &(j5,j6)/zb(j1,j2)*zb(j1,j4)**2/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(
     &j4,j6))**2*t(j3,j5,j6)+L1(-s(j1,j2),-t(j3,j5,j6))*za(j4,j5)**2/za(
     &j5,j6)/zb(j1,j2)*zb(j1,j4)**2/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j
     &4,j6))/t(j3,j5,j6)/2-3._dp/4._dp*(Lnrat(-s(j3,j4),-s(j1,j2))+Lnrat(-
     &t(j3,j5,j6),-s(j1,j2)))/za(j5,j6)/zb(j1,j2)*(za(j3,j5)*zb(j1,j3)-z
     &a(j5,j6)*zb(j1,j6))**2/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))/
     &t(j3,j5,j6)
      s4 = s3+(Lnrat(-s(j3,j4),-s(j1,j2))+Lnrat(-t(j3,j5,j6),-s(j1,j2)))
     &*za(j3,j5)/za(j5,j6)/zb(j1,j2)*zb(j1,j4)*(za(j3,j5)*zb(j1,j3)-za(j
     &5,j6)*zb(j1,j6))/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))**2/2
      s5 = s4+2*L0(-t(j3,j5,j6),-s(j1,j2))/s(j1,j2)*za(j4,j5)/za(j5,j6)/
     &zb(j1,j2)*zb(j1,j4)*(za(j3,j5)*zb(j1,j3)-za(j5,j6)*zb(j1,j6))/(-za
     &(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))
      s2 = s5+Lsm1_2mh(s(j3,j4),t(j3,j5,j6),s(j5,j6),s(j1,j2))*(-1/za(j5
     &,j6)/zb(j1,j2)*(-za(j3,j5)*zb(j1,j5)-za(j3,j6)*zb(j1,j6))**2/(-za(
     &j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))**3*(-za(j3,j5)*zb(j3,j4)-za(
     &j5,j6)*zb(j4,j6))**2/t(j3,j5,j6)+1/za(j1,j2)*za(j2,j4)**2*zb(j3,j6
     &)**2/(-za(j4,j5)*zb(j3,j5)-za(j4,j6)*zb(j3,j6))/zb(j5,j6)/t(j3,j5,
     &j6))+3._dp/4._dp*(Lnrat(-s(j3,j4),-s(j1,j2))+Lnrat(-t(j4,j5,j6),-s(j
     &1,j2)))/za(j1,j2)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))/zb(j5
     &,j6)*(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6))**2/t(j4,j5,j6)
      s4 = s2-(Lnrat(-s(j3,j4),-s(j1,j2))+Lnrat(-t(j4,j5,j6),-s(j1,j2)))
     &/za(j1,j2)*za(j2,j3)*zb(j4,j6)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(
     &j4,j6))**2/zb(j5,j6)*(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6))/2
      s3 = s4+L0(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)/za(j1,j2)*za(j2,j3)**2
     &*zb(j3,j6)*zb(j4,j6)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))**2
     &/zb(j5,j6)*t(j4,j5,j6)-L1(-s(j1,j2),-t(j4,j5,j6))/za(j1,j2)*za(j2,
     &j3)**2*zb(j3,j6)**2/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))/zb(
     &j5,j6)/t(j4,j5,j6)/2-2*L0(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)/za(j1,j
     &2)*za(j2,j3)*zb(j3,j6)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))/
     &zb(j5,j6)*(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6))
      s5 = s3
      s8 = -cone
      s10 = Lnrat(-s(j5,j6),-s(j3,j4))
      s12 = -1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s
     &(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j5,j6)*zb(j1,j6)*(-za(
     &j2,j3)*zb(j3,j6)+za(j2,j4)*zb(j4,j6))/(-za(j3,j5)*zb(j4,j5)-za(j3,
     &j6)*zb(j4,j6))+3._dp/2._dp*(s(j1,j2)-s(j3,j4)-s(j5,j6))/(s(j1,j2)**2
     &-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(
     &j5,j6)+s(j5,j6)**2)**2*(-za(j2,j3)*zb(j1,j3)-za(j2,j4)*zb(j1,j4))/
     &(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))*(-za(j3,j5)*zb(j3,j6)-z
     &a(j4,j5)*zb(j4,j6))*(-t(j3,j5,j6)+t(j4,j5,j6))
      s11 = s12-1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2
     &)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)/za(j1,j2)*za(j2,j5)/(-
     &za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))*(-za(j2,j3)*zb(j3,j6)*t(j
     &3,j5,j6)+za(j2,j4)*zb(j4,j6)*t(j4,j5,j6)+(-s(j1,j2)+s(j3,j4)-s(j5,
     &j6))*za(j2,j3)*zb(j4,j6)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6)
     &)*t(j4,j5,j6))-1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(
     &j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j2,j3)*za(j5,j
     &6)*zb(j4,j6)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))**2*(2*za(j
     &2,j5)*zb(j1,j2)*zb(j5,j6)-zb(j1,j6)*(-t(j3,j5,j6)+t(j4,j5,j6)))
      s9 = s10*s11
      s7 = s8*s9
      s9 = Lnrat(-s(j5,j6),-s(j3,j4))
      s11 = 1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s(
     &j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j2,j5)*(za(j3,j5)*zb(j1
     &,j3)-za(j4,j5)*zb(j1,j4))/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6
     &))*zb(j5,j6)+1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1
     &,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j5)*zb(j1,j4)
     &/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))**2*zb(j5,j6)*(2*za(j1,
     &j2)*za(j5,j6)*zb(j1,j6)-za(j2,j5)*(t(j3,j5,j6)-t(j4,j5,j6)))
      s12 = s11
      s14 = 3._dp/2._dp*(s(j1,j2)-s(j3,j4)-s(j5,j6))/(s(j1,j2)**2-2*s(j1,j
     &2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(
     &j5,j6)**2)**2*(-za(j2,j3)*zb(j1,j3)-za(j2,j4)*zb(j1,j4))/(-za(j3,j
     &5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))*(-za(j3,j5)*zb(j3,j6)-za(j4,j5)*
     &zb(j4,j6))*(t(j3,j5,j6)-t(j4,j5,j6))
      s15 = 1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s(
     &j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)/zb(j1,j2)*zb(j1,j6)/(-za(j
     &3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))*(za(j3,j5)*zb(j1,j3)*t(j3,j5,
     &j6)+(-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j3,j5)*zb(j1,j4)/(-za(j3,j5)*
     &zb(j4,j5)-za(j3,j6)*zb(j4,j6))*t(j3,j5,j6)-za(j4,j5)*zb(j1,j4)*t(j
     &4,j5,j6))
      s13 = s14+s15
      s10 = s12+s13
      s8 = s9*s10
      s6 = s7+s8
      s4 = s5+s6
      s6 = s4
      s8 = i3m(s(j5,j6),s(j3,j4),s(j1,j2))
      s10 = -1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*s
     &(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j1,j2)*(za(j3,j5)*zb(j
     &1,j3)+za(j4,j5)*zb(j1,j4))**2/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j
     &4,j6))*zb(j5,j6)+za(j3,j5)*zb(j1,j4)/(-za(j3,j5)*zb(j4,j5)-za(j3,j
     &6)*zb(j4,j6))**2*(za(j2,j3)*zb(j3,j6)+za(j2,j5)*zb(j5,j6))+(-s(j1,
     &j2)+s(j3,j4)-s(j5,j6))/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**
     &2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j5)*z
     &b(j1,j4)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))**2*((s(j1,j2)-
     &s(j3,j4)-s(j5,j6))*za(j1,j2)*zb(j1,j6)+(-s(j1,j2)-s(j3,j4)+s(j5,j6
     &))*za(j2,j5)*zb(j5,j6))
      s11 = s10-2*za(j2,j4)*(-za(j3,j5)*zb(j1,j5)-za(j3,j6)*zb(j1,j6))*z
     &b(j3,j6)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))*(-za(j3,j5)*zb
     &(j3,j4)-za(j5,j6)*zb(j4,j6))/t(j3,j5,j6)**2
      s9 = s11+s(j3,j4)/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s
     &(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j2,j5)*zb(j1,
     &j6)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))*(t(j3,j5,j6)-t(j4,j
     &5,j6))/2+(-s(j1,j2)**2+2*s(j1,j2)*s(j3,j4)-s(j3,j4)**2+3*s(j3,j4)*
     &(-s(j1,j2)+s(j3,j4)-s(j5,j6))+2*s(j1,j2)*s(j5,j6)+2*s(j3,j4)*s(j5,
     &j6)-s(j5,j6)**2)/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(
     &j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*(-za(j2,j3)*zb
     &(j1,j3)-za(j2,j4)*zb(j1,j4))/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4
     &,j6))*(-za(j3,j5)*zb(j3,j6)-za(j4,j5)*zb(j4,j6))*(t(j3,j5,j6)-t(j4
     &,j5,j6))/2
      s7 = s8*s9
      s5 = s6+s7
      t0 = s5-1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2*s(j1,j2)*
     &s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)/za(j1,j2)*za(j2,j5)/za(j
     &5,j6)*(-(-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j5)+2*za(j1,j2)*za(j5,
     &j6)*zb(j1,j6))/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))*(-t(j3,j
     &5,j6)+t(j4,j5,j6))/4+1/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s(j3,j4)**
     &2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)/zb(j1,j2)*z
     &b(j1,j6)/(-za(j3,j5)*zb(j4,j5)-za(j3,j6)*zb(j4,j6))/zb(j5,j6)*(-(-
     &s(j1,j2)+s(j3,j4)-s(j5,j6))*zb(j1,j6)+2*za(j2,j5)*zb(j1,j2)*zb(j5,
     &j6))*(t(j3,j5,j6)-t(j4,j5,j6))/4

      fsl=t0
      return
      end
