!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Master3(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.7.5
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: Master3,Master3a,zab2,zab,zba,zb12b,zb22b
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lsm1_2mh,I3m,M3bit1,M3bit2,M3bit3,M3bit4
      real(dp):: t,delta12,delta34,delta56,Delta3
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
      zb12b(j1,j2,j3,j4,j5)=zb(j1,j2)*zab2(j2,j3,j4,j5)
      zb22b(j1,j2,j3,j4,j5,j6)=
     & zb(j1,j2)*zab2(j2,j4,j5,j6)+zb(j1,j3)*zab2(j3,j4,j5,j6)
c     End statement functions

      delta12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      delta34=s(j3,j4)-s(j5,j6)-s(j1,j2)
      delta56=s(j5,j6)-s(j1,j2)-s(j3,j4)
      Delta3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2*s(j1,j2)*s(j3,j4)-2*s(j3,j4)*s(j5,j6)-2*s(j5,j6)*s(j1,j2)

      M3bit1=+0.5_dp*za(j1,j5)**2*zb(j1,j4)*zb(j2,j4)*t(j2,j3,j4)
     & /(zb(j2,j3)*zb(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2))
     & *(-2*zb(j2,j4)/zab2(j1,j3,j4,j2))

      M3bit2=+0.5_dp*za(j1,j5)**2*zb(j1,j4)*zb(j2,j4)*t(j2,j3,j4)
     & /(zb(j2,j3)*zb(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2))
     & *zb(j1,j4)

      M3bit3=-0.5_dp*za(j2,j3)*zb(j2,j4)*zab2(j5,j3,j4,j2)
     & /(zb(j2,j3)*za(j5,j6)*zab2(j1,j3,j4,j2))
     & *((zab2(j5,j2,j3,j4)/t(j2,j3,j4)
     & -2*za(j1,j5)*zb(j2,j4)/zab2(j1,j3,j4,j2)))

      M3bit4=-0.5_dp*za(j2,j3)*zb(j2,j4)*zab2(j5,j3,j4,j2)
     & /(zb(j2,j3)*za(j5,j6)*zab2(j1,j3,j4,j2))
     & *(zab(j5,j2,j4))

      Master3=Master3a(j1,j2,j3,j4,j5,j6,za,zb)
     & -za(j1,j5)**2*zb(j2,j4)**3*t(j2,j3,j4)
     & /(zb(j2,j3)*zb(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2)**3)
     & *Lsm1_2mh(s(j1,j2),t(j2,j3,j4),s(j3,j4),s(j5,j6))

     & +zb(j2,j4)/(zb(j2,j3)*zab2(j1,j3,j4,j2)*Delta3)
     & *(1.5_dp*s(j1,j2)*delta12*(t(j1,j3,j4)-t(j2,j3,j4))
     & *zab2(j3,j1,j2,j4)*zab2(j5,j1,j2,j6)/Delta3

     & +0.5_dp*s(j1,j2)*zab2(j3,j1,j2,j4)*zab2(j5,j1,j2,j6)

     & -za(j1,j3)*za(j1,j5)*zb(j2,j4)*zab2(j2,j3,j4,j1)
     & /zab2(j1,j3,j4,j2)
     & *(zb12b(j2,j1,j3,j4,j6)-zb22b(j2,j3,j4,j1,j2,j6))

     & +zab(j2,j1,j4)*(zab2(j3,j1,j4,j2)*zab2(j5,j1,j2,j6)
     & +zab(j5,j1,j2)*zab(j3,j4,j6) )

     & +za(j1,j5)*za(j2,j3)*(zb(j2,j4)
     & *(zb12b(j1,j2,j3,j4,j6)-zb22b(j1,j3,j4,j1,j2,j6))

     & -zb(j1,j2)*zb22b(j4,j1,j2,j3,j4,j6)))
     & *I3m(s(j1,j2),s(j3,j4),s(j5,j6))

     & +M3bit1*L0(-s(j5,j6),-t(j2,j3,j4))/t(j2,j3,j4)
     & +M3bit2*L1(-s(j5,j6),-t(j2,j3,j4))/t(j2,j3,j4)**2
     & +M3bit3*L0(-t(j2,j3,j4),-s(j3,j4))/s(j3,j4)
     & +M3bit4*L1(-s(j3,j4),-t(j2,j3,j4))/t(j2,j3,j4)**2

     & +0.5_dp*zb(j2,j4)/(zb(j2,j3)*zab2(j1,j3,j4,j2)*Delta3)
     & *(-za(j3,j5)/za(j3,j4)*(za(j2,j5)/za(j5,j6)
     & *(zab(j3,j1,j2)*delta12-zab(j3,j4,j2)*delta34)
     & -zb(j1,j6)*(za(j1,j3)*delta56-2*zab(j1,j2,j4)*za(j4,j3)))

     & +zb(j4,j6)/zb(j3,j4)*(zb(j1,j6)/zb(j5,j6)
     & *(zab(j1,j2,j4)*delta12-zab(j1,j3,j4)*delta34)
     & +za(j2,j5)*(zb(j2,j4)*delta56-2*zba(j2,j1,j3)*zb(j3,j4))))

      return
      end



