!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function a6(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     July, 1998.                                                      *
c     st is a string to choose between pp,pm or sl
c     implementation of Eq. (3.3) of BDKW hep-ph/9610370
c     character string st can take the value pp,pm or sl
c***********************************************************************
      include 'mxpart.f'
      include 'zprods_decl.f'
      character(len=2):: st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: atree,vv,switchyard

      a6=atree(st,j1,j2,j3,j4,j5,j6,za,zb)*vv(st,j1,j2,j3,j4,j5,j6)
      a6=a6+switchyard(st,j1,j2,j3,j4,j5,j6,za,zb)

      end

      function switchyard(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: switchyard
c-----switchyard function to direct to pm,pp or st
      include 'mxpart.f'
      include 'zprods_decl.f'
      character(len=2):: st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: fpp,fpm,fsl
      if     (st == 'pp') then
           switchyard=fpp(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st == 'pm') then
           switchyard=fpm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st == 'sl') then
           switchyard=fsl(j1,j2,j3,j4,j5,j6,za,zb)
      endif
      end



