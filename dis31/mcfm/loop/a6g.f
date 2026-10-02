!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function a6g(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6g

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     August, 1999.                                                    *
c     implementation of Eq. (6.1) of BDKW hep-ph/9708239               *
c     character string st can take various values                      *
c***********************************************************************
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6treeg,vvg,fcc,fsc

      a6g=a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)*vvg(st,j1,j2,j3,j4,j5,j6)
      a6g=a6g+fcc(st,j1,j2,j3,j4,j5,j6,za,zb)
     &       +fsc(st,j1,j2,j3,j4,j5,j6,za,zb)

      end

