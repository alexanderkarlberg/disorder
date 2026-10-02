!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function fpm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: fpm
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lnrat,Lsm1,i3m,Lsm1_2mh
      complex(dp):: t0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,
     & s11,s12,s13,s14,s15,s16,s17,s18
      real(dp):: t
      s3 = -2*L0(-t(j2,j3,j4),-s(j5,j6))/s(j5,j6)*za(j1,j5)/za(j5,j6)*zb
     &(j1,j3)/zb(j2,j3)/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*(-za(j
     &2,j5)*zb(j2,j3)+za(j4,j5)*zb(j3,j4))-2*L0(-t(j1,j2,j3),-s(j5,j6))/
     &s(j5,j6)/za(j2,j3)*za(j2,j4)/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,
     &j4))*(-za(j1,j2)*zb(j1,j6)+za(j2,j3)*zb(j3,j6))*zb(j4,j6)/zb(j5,j6
     &)
      s4 = s3+Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))*(-za(j
     &4,j5)**2/za(j5,j6)*zb(j1,j3)**2/(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(
     &j1,j3))/zb(j2,j3)/t(j1,j2,j3)+za(j1,j2)**2/za(j1,j3)**2/za(j2,j3)*
     &(-za(j1,j3)*zb(j1,j6)-za(j2,j3)*zb(j2,j6))**2/(za(j1,j2)*zb(j2,j4)
     &+za(j1,j3)*zb(j3,j4))/zb(j5,j6)/t(j1,j2,j3))
      s2 = s4+Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))*(-za(j4,j
     &5)**2/za(j5,j6)*zb(j1,j3)**2/(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,
     &j3))/zb(j2,j3)/t(j1,j2,j3)+1/za(j2,j3)/(-za(j1,j3)*zb(j1,j4)-za(j2
     &,j3)*zb(j2,j4))**2*(-za(j1,j3)*zb(j1,j6)-za(j2,j3)*zb(j2,j6))**2/(
     &za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*(-za(j1,j2)*zb(j1,j4)+za(
     &j2,j3)*zb(j3,j4))**2/zb(j5,j6)/t(j1,j2,j3))+L0(-t(j1,j2,j3),-s(j1,
     &j2))/s(j1,j2)*za(j1,j2)/za(j1,j3)*zb(j1,j3)/(-za(j1,j3)*zb(j1,j4)-
     &za(j2,j3)*zb(j2,j4))*(-za(j1,j3)*zb(j1,j6)-za(j2,j3)*zb(j2,j6))**2
     &/zb(j5,j6)/t(j1,j2,j3)
      s3 = s2+L0(-t(j1,j2,j3),-s(j2,j3))/s(j2,j3)*za(j1,j2)/za(j1,j3)*zb
     &(j1,j3)*(-za(j1,j3)*zb(j1,j6)-za(j2,j3)*zb(j2,j6))/(za(j1,j2)*zb(j
     &2,j4)+za(j1,j3)*zb(j3,j4))*(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6
     &))/zb(j5,j6)/t(j1,j2,j3)-L1(-t(j1,j2,j3),-s(j2,j3))/s(j2,j3)**2*za
     &(j2,j3)*zb(j1,j3)**2/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*(za
     &(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))**2/zb(j5,j6)/t(j1,j2,j3)/2
      s4 = s3+L0(-t(j1,j2,j3),-s(j2,j3))/s(j2,j3)*zb(j1,j3)/(za(j1,j2)*z
     &b(j2,j4)+za(j1,j3)*zb(j3,j4))*(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3
     &,j6))*(-za(j1,j2)*zb(j1,j6)+za(j2,j3)*zb(j3,j6))/zb(j5,j6)/t(j1,j2
     &,j3)
      s5 = s4-L1(-s(j5,j6),-t(j1,j2,j3))/za(j2,j3)*za(j2,j4)**2/(za(j1,j
     &2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j4,j6)**2/zb(j5,j6)/t(j1,j2,j
     &3)/2
      s6 = s5
      s8 = L0(-t(j1,j2,j3),-s(j5,j6))/s(j5,j6)*za(j2,j4)/(-za(j1,j3)*zb(
     &j1,j4)-za(j2,j3)*zb(j2,j4))/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j
     &4))*zb(j4,j6)**2/zb(j5,j6)*t(j1,j2,j3)
      s10 = i3m(s(j1,j2),s(j3,j4),s(j5,j6))
      s12 = -1/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))/(za(j1,j2)*zb(
     &j2,j4)+za(j1,j3)*zb(j3,j4))*(za(j2,j3)*zb(j3,j4)*(-za(j3,j5)*zb(j3
     &,j6)-za(j4,j5)*zb(j4,j6))-za(j1,j2)*zb(j1,j4)*(-za(j3,j5)*zb(j3,j6
     &)+za(j4,j5)*zb(j4,j6)))/2-2*za(j1,j2)*za(j4,j5)*zb(j1,j3)/(za(j1,j
     &2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j4,j6)/t(j1,j2,j3)
      s13 = s12
      s15 = (za(j4,j5)**2/za(j5,j6)*zb(j1,j3)**2/(za(j2,j4)*zb(j1,j2)+za
     &(j3,j4)*zb(j1,j3))/zb(j2,j3)-1/za(j2,j3)/(za(j1,j2)*zb(j2,j4)+za(j
     &1,j3)*zb(j3,j4))*(-za(j1,j2)*zb(j1,j6)+za(j2,j3)*zb(j3,j6))**2/zb(
     &j5,j6))/t(j1,j2,j3)**2*(2*s(j1,j2)*s(j5,j6)+(-s(j1,j2)+s(j3,j4)-s(
     &j5,j6))*t(j1,j2,j3))/2
      s16 = (-s(j1,j2)-s(j3,j4)+s(j5,j6))/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j
     &4)+s(j3,j4)**2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2
     &)/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(-za(j1,j5)*zb(j1,j6)
     &-za(j2,j5)*zb(j2,j6))/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*(-
     &s(j5,j6)*za(j1,j2)*zb(j1,j4)-(-za(j1,j2)*zb(j1,j4)+za(j2,j3)*zb(j3
     &,j4))*t(j1,j2,j3))/2
      s14 = s15+s16
      s11 = s13+s14
      s9 = s10*s11
      s7 = s8+s9
      s1 = s6+s7
      s4 = s1
      s6 = Lnrat(-s(j1,j2),-s(j5,j6))/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+s
     &(j3,j4)**2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)/(-
     &za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*(2*s(j1,j2)*(-za(j2,j5)*z
     &b(j1,j6)+za(j1,j5)*(-za(j2,j3)*zb(j1,j3)-za(j2,j4)*zb(j1,j4))/(-za
     &(j1,j3)*zb(j2,j3)-za(j1,j4)*zb(j2,j4))*zb(j2,j6))+za(j1,j2)*(zb(j1
     &,j6)**2+(-za(j2,j3)*zb(j1,j3)-za(j2,j4)*zb(j1,j4))/(-za(j1,j3)*zb(
     &j2,j3)-za(j1,j4)*zb(j2,j4))*zb(j2,j6)**2)/zb(j5,j6)*t(j1,j2,j3)+1/
     &za(j5,j6)*zb(j1,j2)*(za(j2,j5)**2+za(j1,j5)**2*(-za(j2,j3)*zb(j1,j
     &3)-za(j2,j4)*zb(j1,j4))/(-za(j1,j3)*zb(j2,j3)-za(j1,j4)*zb(j2,j4))
     &)*t(j1,j2,j4))
      s7 = -Lsm1(-s(j3,j4),-t(j2,j3,j4),-s(j2,j3),-t(j2,j3,j4))*(1/za(j5
     &,j6)/zb(j2,j3)/zb(j2,j4)**2*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j
     &4))**2*zb(j3,j4)**2/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))/t(j2
     &,j3,j4)-1/za(j2,j3)*za(j2,j4)**2/(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb
     &(j1,j3))*zb(j1,j6)**2/zb(j5,j6)/t(j2,j3,j4))
      s5 = s6+s7
      s3 = s4+s5
      s4 = s3-Lsm1_2mh(s(j1,j2),t(j2,j3,j4),s(j3,j4),s(j5,j6))*(1/za(j5,
     &j6)/zb(j2,j3)/(-za(j1,j3)*zb(j2,j3)-za(j1,j4)*zb(j2,j4))**2*(za(j3
     &,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))**2/(za(j1,j2)*zb(j2,j4)+za(j1,
     &j3)*zb(j3,j4))*(za(j1,j2)*zb(j2,j3)-za(j1,j4)*zb(j3,j4))**2/t(j2,j
     &3,j4)-1/za(j2,j3)*za(j2,j4)**2/(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j
     &1,j3))*zb(j1,j6)**2/zb(j5,j6)/t(j2,j3,j4))
      s2 = s4-L0(-t(j2,j3,j4),-s(j3,j4))/s(j3,j4)*za(j2,j4)/za(j5,j6)/zb
     &(j2,j4)/(-za(j1,j3)*zb(j2,j3)-za(j1,j4)*zb(j2,j4))*(za(j3,j5)*zb(j
     &2,j3)+za(j4,j5)*zb(j2,j4))**2*zb(j3,j4)/t(j2,j3,j4)+L1(-s(j5,j6),-
     &t(j2,j3,j4))*za(j1,j5)**2/za(j5,j6)*zb(j1,j3)**2/zb(j2,j3)/(za(j1,
     &j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))/t(j2,j3,j4)/2
      s3 = s2-L0(-t(j2,j3,j4),-s(j2,j3))/s(j2,j3)*za(j2,j4)/za(j5,j6)/zb
     &(j2,j4)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))*zb(j3,j4)/(za(j1
     &,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*(-za(j2,j5)*zb(j2,j4)-za(j3,j5
     &)*zb(j3,j4))/t(j2,j3,j4)+L1(-t(j2,j3,j4),-s(j2,j3))/s(j2,j3)**2*za
     &(j2,j4)**2/za(j5,j6)*zb(j2,j3)/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j
     &3,j4))*(-za(j2,j5)*zb(j2,j4)-za(j3,j5)*zb(j3,j4))**2/t(j2,j3,j4)/2
      s4 = s3-L0(-t(j2,j3,j4),-s(j2,j3))/s(j2,j3)*za(j2,j4)/za(j5,j6)/(z
     &a(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*(-za(j2,j5)*zb(j2,j4)-za(j
     &3,j5)*zb(j3,j4))*(-za(j2,j5)*zb(j2,j3)+za(j4,j5)*zb(j3,j4))/t(j2,j
     &3,j4)
      s5 = s4-L0(-t(j2,j3,j4),-s(j5,j6))/s(j5,j6)*za(j1,j5)**2/za(j5,j6)
     &*zb(j1,j3)/(-za(j1,j3)*zb(j2,j3)-za(j1,j4)*zb(j2,j4))/(za(j1,j2)*z
     &b(j2,j4)+za(j1,j3)*zb(j3,j4))*t(j2,j3,j4)
      s6 = s5
      s8 = -Lnrat(-s(j3,j4),-s(j5,j6))/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j4)+
     &s(j3,j4)**2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)/(
     &-za(j1,j3)*zb(j2,j3)-za(j1,j4)*zb(j2,j4))*(2*s(j3,j4)*(-za(j4,j5)*
     &zb(j3,j6)+za(j3,j5)*(-za(j1,j4)*zb(j1,j3)-za(j2,j4)*zb(j2,j3))/(-z
     &a(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4))*zb(j4,j6))+za(j3,j4)*(zb(j
     &3,j6)**2+(-za(j1,j4)*zb(j1,j3)-za(j2,j4)*zb(j2,j3))/(-za(j1,j3)*zb
     &(j1,j4)-za(j2,j3)*zb(j2,j4))*zb(j4,j6)**2)/zb(j5,j6)*t(j1,j3,j4)+1
     &/za(j5,j6)*(za(j4,j5)**2+za(j3,j5)**2*(-za(j1,j4)*zb(j1,j3)-za(j2,
     &j4)*zb(j2,j3))/(-za(j1,j3)*zb(j1,j4)-za(j2,j3)*zb(j2,j4)))*zb(j3,j
     &4)*t(j2,j3,j4))
      s10 = -cone
      s12 = i3m(s(j3,j4),s(j1,j2),s(j5,j6))
      s14 = -1/(-za(j1,j3)*zb(j2,j3)-za(j1,j4)*zb(j2,j4))/(za(j1,j2)*zb(
     &j2,j4)+za(j1,j3)*zb(j3,j4))*(za(j1,j2)*zb(j2,j3)*(-za(j1,j5)*zb(j1
     &,j6)-za(j2,j5)*zb(j2,j6))-za(j1,j4)*(za(j1,j5)*zb(j1,j6)-za(j2,j5)
     &*zb(j2,j6))*zb(j3,j4))/2-2*za(j1,j5)*za(j2,j4)*zb(j1,j6)*zb(j3,j4)
     &/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))/t(j2,j3,j4)
      s15 = s14
      s17 = (-1/za(j5,j6)/zb(j2,j3)/(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3
     &,j4))*(-za(j2,j5)*zb(j2,j3)+za(j4,j5)*zb(j3,j4))**2+1/za(j2,j3)*za
     &(j2,j4)**2/(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j1,j6)**2/
     &zb(j5,j6))/t(j2,j3,j4)**2*(2*s(j3,j4)*s(j5,j6)+(s(j1,j2)-s(j3,j4)-
     &s(j5,j6))*t(j2,j3,j4))/2
      s18 = (-s(j1,j2)-s(j3,j4)+s(j5,j6))/(s(j1,j2)**2-2*s(j1,j2)*s(j3,j
     &4)+s(j3,j4)**2-2*s(j1,j2)*s(j5,j6)-2*s(j3,j4)*s(j5,j6)+s(j5,j6)**2
     &)/(-za(j1,j3)*zb(j2,j3)-za(j1,j4)*zb(j2,j4))/(za(j1,j2)*zb(j2,j4)+
     &za(j1,j3)*zb(j3,j4))*(-za(j3,j5)*zb(j3,j6)-za(j4,j5)*zb(j4,j6))*(-
     &s(j5,j6)*za(j1,j4)*zb(j3,j4)-(za(j1,j2)*zb(j2,j3)-za(j1,j4)*zb(j3,
     &j4))*t(j2,j3,j4))/2
      s16 = s17+s18
      s13 = s15+s16
      s11 = s12*s13
      s9 = s10*s11
      s7 = s8+s9
      t0 = s6+s7
      fpm=t0
      return
      end
