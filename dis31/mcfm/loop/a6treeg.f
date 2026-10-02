!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6treeg

c----hep-ph/9708239, Eqs(8.4,8.8,8.14,9.2,9.7,10.1) Multiplied by (-i)
      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'heldefs.f'
      integer st
      real(dp):: t

      if(st==hqpgmgpqbm) then
      a6treeg=
     &((za(j2,j4)*za(j4,j5)*zb(j1,j3)*zb(j1,j6))/
     &(s(j2,j3)*s(j5,j6)*za(j3,j4)*zb(j1,j2))-
     &(za(j4,j5)*zb(j1,j3)**2*(-(za(j1,j2)*zb(j1,j6))+za(j2,j3)*zb(j3,j6
     &)))/
     &(s(j2,j3)*s(j5,j6)*zb(j1,j2)*t(j1,j2,j3))+
     &(za(j2,j4)**2*zb(j1,j6)*(-(za(j2,j5)*zb(j2,j3))+za(j4,j5)*zb(j3,j4
     &)))/
     &(s(j2,j3)*s(j5,j6)*za(j3,j4)*t(j2,j3,j4)))
      elseif(st==hqpgpgmqbm) then
      a6treeg=
     &(-(((za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))*
     &(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6)))/
     &(s(j2,j3)*s(j5,j6)*za(j1,j2)*zb(j3,j4)))-
     &(za(j1,j3)*za(j4,j5)*zb(j1,j2)*
     &(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6)))/
     &(s(j2,j3)*s(j5,j6)*za(j1,j2)*t(j1,j2,j3))+
     &(za(j3,j4)*zb(j1,j6)*zb(j2,j4)*
     &(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4)))/
     &(s(j2,j3)*s(j5,j6)*zb(j3,j4)*t(j2,j3,j4)))
      elseif(st==hqpgpgpqbm) then
c---This amplitude corresponds to
c    q(4)+l(5) --> q_R(1)+l_R(6)+g_R(2)+g_R(3)
      a6treeg=
     &(-za(j4,j5)**2)/(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j5,j6))
      elseif(st==hqpgpqbmgm) then
      a6treeg=
     &(((za(j3,j5)*zb(j1,j3)+za(j4,j5)*zb(j1,j4))*
     &(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6)))/
     &(s(j5,j6)*za(j1,j2)*za(j2,j3)*zb(j1,j4)*zb(j3,j4))+
     &(za(j3,j5)*zb(j1,j2)*(-(za(j1,j4)*zb(j1,j6))-za(j2,j4)*zb(j2,j6)))
     &/
     &(s(j5,j6)*za(j1,j2)*zb(j1,j4)*t(j1,j2,j4))-
     &(za(j3,j4)*zb(j1,j6)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4)))/
     &(s(j5,j6)*za(j2,j3)*zb(j3,j4)*t(j2,j3,j4)))
      elseif(st==hqpgpqbmgp) then
      a6treeg=
     &(-za(j1,j3)*za(j3,j5)**2)/
     &(za(j1,j2)*za(j1,j4)*za(j2,j3)*za(j3,j4)*za(j5,j6))
      elseif(st==hqpqbmgmgp) then
      a6treeg=
     &(((-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4))*
     &(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6)))/
     &(s(j3,j4)*s(j5,j6)*za(j1,j4)*zb(j2,j3))-
     &(za(j1,j3)*za(j2,j5)*zb(j1,j4)*
     &(-(za(j1,j3)*zb(j1,j6))+za(j3,j4)*zb(j4,j6)))/
     &(s(j3,j4)*s(j5,j6)*za(j1,j4)*t(j1,j3,j4))-
     &(za(j2,j3)*zb(j1,j6)*zb(j2,j4)*
     &(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4)))/
     &(s(j3,j4)*s(j5,j6)*zb(j2,j3)*t(j2,j3,j4)))
      elseif(st==hqpqbmgpgm) then
      a6treeg=
     &-(-((za(j2,j4)*za(j2,j5)*zb(j1,j3)*zb(j1,j6))/
     &(s(j3,j4)*s(j5,j6)*za(j2,j3)*zb(j1,j4)))+
     &(za(j2,j5)*zb(j1,j3)**2*(-(za(j1,j4)*zb(j1,j6))-za(j3,j4)*zb(j3,j6
     &)))/
     &(s(j3,j4)*s(j5,j6)*zb(j1,j4)*t(j1,j3,j4))+
     &(za(j2,j4)**2*zb(j1,j6)*(-(za(j2,j5)*zb(j2,j3))+za(j4,j5)*zb(j3,j4
     &)))/
     &(s(j3,j4)*s(j5,j6)*za(j2,j3)*t(j2,j3,j4)))
      elseif(st==hqpqbmgpgp) then
      a6treeg=
     &(-za(j2,j5)**2)/(za(j1,j4)*za(j2,j3)*za(j3,j4)*za(j5,j6))
      else
      write(6,*) 'unimplemented st',st
      stop
      endif
      return
      end
