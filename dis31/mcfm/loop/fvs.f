!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Fvs(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Fvs
c--- This function returns the result of BDK
c--- Published in Nucl.Phys.B513:3-86,1998.
c--- e-Print: hep-ph/9708239
c--- Eqs.(11.1) and Eqs.(11.6) as appropriate dependent on choice of "st"

      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'heldefs.f'
      integer st
      complex(dp):: L0,L1,Lsm1_2me,Lnrat
      complex(dp):: Brackpm,Brackpp,Brackppa
      real(dp):: t
c---statement functions
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)

c---This is the square bracket in Eq.(11.1) subject to exchange 16,25
      Brackppa(j1,j2,j3,j4,j5,j6)=
     & za(j1,j2)*(za(j5,j3)*zb(j3,j1))**2/za(j5,j6)/za(j3,j4)**2
     & *L1(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)**2
     & +za(j5,j2)*za(j5,j3)*zb(j3,j1)/za(j5,j6)/za(j3,j4)**2
     & *L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)

c---This is the whole expression in Eq.(11.1) subject to exchange 34
      Brackpp(j1,j2,j3,j4,j5,j6)=
     & -za(j2,j3)*za(j2,j4)*za(j3,j5)*za(j4,j5)
     & /za(j1,j2)/za(j5,j6)/za(j3,j4)**4
     & *Lsm1_2me(t(j1,j2,j3),t(j1,j2,j4),s(j1,j2),s(j5,j6))
     &  -(za(j2,j4)*za(j3,j5)+za(j2,j3)*za(j4,j5))
     & /za(j1,j2)/za(j5,j6)/za(j3,j4)**3
     & *(za(j5,j3)*zb(j3,j1)*za(j1,j2)
     & *L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)
     & +za(j5,j6)*zb(j6,j4)*za(j4,j2)
     & *L0(-t(j1,j2,j3),-s(j5,j6))/s(j5,j6)
     & +0.5_dp*za(j5,j2)*(Lnrat(-t(j1,j2,j3),-t(j1,j2,j4))))
     & -Brackppa(j1,j2,j3,j4,j5,j6)-Brackppa(j6,j5,j3,j4,j2,j1)
     & +0.5_dp/za(j3,j4)**2
     & *(za(j2,j5)**2/za(j1,j2)/za(j5,j6)
     & -zb(j1,j6)**2/zb(j1,j2)/zb(j5,j6))

c---end statement functions

      if(st==hqpqbmgpgm) then
c---flip2:( 1<-->2, 3<-->4, 5<-->6, za<-->zb)
      Fvs=
     & +Brackpm(j1,j2,j3,j4,j5,j6,za,zb)
     & +Brackpm(j2,j1,j4,j3,j6,j5,zb,za)

      elseif(st==hqpqbmgpgp) then

c---exch34:(3<-->4)
      Fvs=
     & +Brackpp(j1,j2,j3,j4,j5,j6)
     & +Brackpp(j1,j2,j4,j3,j5,j6)

      endif
      return
      end

      function Brackpm(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Brackpm
c---This is the whole expression in Eq.(11.6) subject to exchange flip_2
c---ie  flip2:( 1<-->2, 3<-->4, 5<-->6, za<-->zb)

      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: Lsm1_2mh,I3m,Lnrat,zab2,Brackpma
      real(dp):: t,delta,IDelta

c---statement functions
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      delta(j1,j2,j3,j4,j5,j6)=s(j1,j2)-s(j3,j4)-s(j5,j6)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      IDelta=1._dp/(s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2._dp*(+s(j1,j2)*s(j3,j4)+s(j1,j2)*s(j5,j6)+s(j5,j6)*s(j3,j4)))

      Brackpm=
     & -2._dp*za(j2,j3)*zb(j4,j6)*zab2(j3,j1,j2,j6)
     & *zab2(j2,j5,j6,j4)*t(j1,j2,j3)
     & /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**4
     & *Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))
     & +(2._dp*za(j2,j3)*zb(j4,j6)
     & *(zb(j1,j2)*za(j5,j6)*za(j2,j3)*zb(j4,j6)
     & -zab2(j5,j2,j4,j1)*zab2(j3,j1,j2,j4))*zab2(j4,j1,j2,j3)
     & *(t(j1,j2,j3)-t(j1,j2,j4))/zab2(j3,j1,j2,j4)**3*IDelta
     & -3._dp*(s(j3,j4)*delta(j3,j4,j1,j2,j5,j6)*zab2(j2,j3,j4,j1)
     &  *zab2(j5,j3,j4,j6)*IDelta
     & -za(j2,j3)*zb(j3,j1)*za(j5,j4)*zb(j4,j6)
     & -za(j2,j4)*zb(j4,j1)*za(j5,j3)*zb(j3,j6))
     & *zab2(j4,j1,j2,j3)/zab2(j3,j1,j2,j4)*IDelta
     & -zb(j1,j4)*za(j3,j5)*za(j2,j3)*zb(j4,j6)
     & *zab2(j4,j1,j2,j3)**2/zab2(j3,j1,j2,j4)**2*IDelta
     & -zb(j1,j3)*za(j4,j5)*za(j2,j4)*zb(j3,j6)*IDelta)
     & *I3m(s(j1,j2),s(j3,j4),s(j5,j6))
     & +(2._dp*za(j2,j3)*zb(j4,j6)*t(j1,j2,j3)*zab2(j2,j1,j4,j6)
     & /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**3
     & -(zab2(j2,j1,j4,j6)**2
     &  +2._dp*za(j2,j3)*zb(j3,j6)*za(j2,j4)*zb(j4,j6))
     & /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**2)
     & *Lnrat(-t(j1,j2,j3),-s(j3,j4))
c---ie  flip3:( 1<-->5, 2<-->6, 3<-->4, za<-->zb)
     &  +Brackpma(j1,j2,j3,j4,j5,j6,za,zb)
     &  +Brackpma(j5,j6,j4,j3,j1,j2,zb,za)
     &  +zb(j1,j6)*zab2(j4,j1,j2,j3)
     &  *(zb(j1,j6)*delta(j3,j4,j5,j6,j1,j2)
     & -2._dp*zb(j1,j2)*za(j2,j5)*zb(j5,j6))
     & /zb(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)*IDelta
     &  +(zab2(j2,j1,j4,j6)**2+za(j2,j1)*zb(j1,j6)*za(j2,j5)*zb(j5,j6))
     & /za(j1,j2)/zb(j5,j6)/zab2(j3,j1,j2,j4)**2

       return
       end


      function Brackpma(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Brackpma
c---This is the curly braces expression in Eq.(11.6) subject to flip_3
c---ie  flip3:( 1<-->5, 2<-->6, 3<-->4, za<-->zb)

      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: L0,L1,Lnrat,zab2
      real(dp):: t,delta,IDelta

c---statement functions
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j3,j1)
      delta(j1,j2,j3,j4,j5,j6)=s(j1,j2)-s(j3,j4)-s(j5,j6)
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      IDelta=1._dp/(s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2._dp*(s(j1,j2)*s(j3,j4)+s(j1,j2)*s(j5,j6)+s(j5,j6)*s(j3,j4)))

       Brackpma=
     & 2._dp*zb(j1,j2)*za(j2,j3)*zb(j3,j6)
     & /zb(j5,j6)/zab2(j3,j1,j2,j4)**2
     & *(zab2(j3,j1,j2,j6)*zab2(j2,j1,j3,j4)/zab2(j3,j1,j2,j4)
     & *L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)
     & -za(j2,j3)*zb(j3,j6)
     & *(L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)
     & -1/2._dp*L1(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)))
     & -(3._dp*delta(j5,j6,j3,j4,j1,j2)*zab2(j2,j3,j4,j1)
     & *zab2(j5,j3,j4,j6)*zab2(j4,j1,j2,j3)
     & /zab2(j3,j1,j2,j4)*IDelta**2
     & +2._dp*za(j3,j2)*zb(j2,j1)*zb(j4,j6)*(t(j1,j2,j3)-t(j1,j2,j4))
     & *(za(j2,j5)*t(j1,j2,j3)+za(j2,j1)*zb(j1,j6)*za(j6,j5))
     & /zab2(j3,j1,j2,j4)**3*IDelta
     & -(2._dp*za(j2,j3)*zb(j3,j6)*(t(j1,j2,j3)-t(j1,j2,j4))
     & +delta(j1,j2,j3,j4,j5,j6)
     & *(za(j2,j5)*t(j1,j2,j3)
     & +za(j2,j1)*zb(j1,j6)*za(j6,j5))/za(j5,j6))
     & *za(j5,j2)*zb(j2,j1)/zab2(j3,j1,j2,j4)**2*IDelta
     & +0.5_dp*(za(j1,j2)*zb(j1,j6)**2/zb(j5,j6)
     &  +zb(j1,j2)*za(j2,j5)**2/za(j5,j6)
     & -2._dp*za(j2,j5)*zb(j1,j6))*zab2(j4,j1,j2,j3)
     & /zab2(j3,j1,j2,j4)*IDelta)
     & *(Lnrat(-s(j1,j2),-s(j3,j4)))

      return
      end


