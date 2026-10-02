!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function atree(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: atree

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     July, 1998.                                                      *
c***********************************************************************
c---- the basic process is atreepp which is the amplitude for
c---- q-(p4)+Q(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
c---- All outgoing particles are right-handed,
c---- except Q(p3) which is left-handed
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      character(len=2):: st
      real(dp):: t
      complex(dp):: tree
c-----Note all tree amplitudes divided by im wrt to BDKW paper
c-----Statement function for pp (others derived by crossing)
c-----This is an implementation of BDKW Eq.(3.10)
      tree(j1,j2,j3,j4,j5,j6)=
     & +zb(j1,j2)*za(j5,j4)*(za(j3,j1)*zb(j1,j6)+za(j3,j2)*zb(j2,j6))
     &      /(s(j2,j3)*s(j5,j6)*t(j1,j2,j3))
     & +za(j3,j4)*zb(j6,j1)*(za(j5,j3)*zb(j3,j2)+za(j5,j4)*zb(j4,j2))
     &       /(s(j2,j3)*s(j5,j6)*t(j2,j3,j4))

      if (st == 'pp') then
         atree=+tree(j1,j2,j3,j4,j5,j6)
      elseif (st == 'pm') then
         atree=-tree(j1,j3,j2,j4,j5,j6)
      elseif (st == 'sl') then
c---this is not a real tree amplitude but it is an auxiliary quantity.
         atree=-tree(j3,j1,j2,j4,j5,j6)
      endif

      return
      end
