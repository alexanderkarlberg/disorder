!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function atreez(i1,i2,i3,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: atreez
c--- this function should be passed the (true) helicities of the particles,
c--- namely with i1 incoming and i2,i3 outgoing (1=left, 2=right)
c--- to compare with BDKW paper, first flip i1 to outgoing and then
c--- note that left corresponds to + and right to -
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,nhel
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: atree

      nhel=(i1-1)+(i2-1)*2+(i3-1)*4

c--- note that reversing all helicities swaps za,zb
c--- and switching 3rd helicity swaps j5 and j6
      if     (nhel == 0) then
        atreez=+atree('pm',j1,j2,j3,j4,j5,j6,zb,za)
      elseif (nhel == 1) then
        atreez=+atree('pp',j1,j2,j3,j4,j6,j5,za,zb)
      elseif (nhel == 2) then
        atreez=+atree('pp',j1,j2,j3,j4,j5,j6,zb,za)
      elseif (nhel == 3) then
        atreez=+atree('pm',j1,j2,j3,j4,j6,j5,za,zb)
      elseif (nhel == 4) then
        atreez=+atree('pm',j1,j2,j3,j4,j6,j5,zb,za)
      elseif (nhel == 5) then
        atreez=+atree('pp',j1,j2,j3,j4,j5,j6,za,zb)
      elseif (nhel == 6) then
        atreez=+atree('pp',j1,j2,j3,j4,j6,j5,zb,za)
      elseif (nhel == 7) then
        atreez=+atree('pm',j1,j2,j3,j4,j5,j6,za,zb)
      else
        write(*,*) 'Unsupported helicity assignment'
        stop
      endif

      return
      end
