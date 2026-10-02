!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function a61(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a61
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     July, 1998.                                                      *
c***********************************************************************
c     implementation of Eqs. (2.7) and (2.8) of BDKW hep-ph/9610370
c     with ns=0
c     character string st can take the value pp or pm
c     a61(pp,j1,j2,j3,j4,j5,j6,za,zb) corresponds to
c     q(j1,+)+Q(j3,-)+e(j6,+)+q~(j4)+Q~(j2)+e~(j5)
c     a61(pm,j1,j2,j3,j4,j5,j6,za,zb) corresponds to
c     q(j1,+)+Q(j3,+)+e(j6,+)+q~(j4)+Q~(j2)+e~(j5)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      character(len=2):: st
      complex(dp):: a6,aa6sf,aa6tp,aa6uv

c----Includes ultraviolet subtraction aa6uv
      call a6routine(st,j1,j2,j3,j4,j5,j6,za,zb,aa6sf,aa6tp,aa6uv)
      a61=(one-two/xnsq)*a6(st,j1,j2,j3,j4,j5,j6,za,zb)
     & -(real(nf,dp)*aa6sf-aa6tp)/xn-aa6uv

      if (st == 'pp') then
      a61=a61+(-two*a6('pm',j1,j3,j2,j4,j5,j6,za,zb)
     &             +a6('sl',j2,j3,j1,j4,j5,j6,za,zb))/xnsq
      elseif (st == 'pm') then
      a61=a61+(-two*a6('pp',j1,j3,j2,j4,j5,j6,za,zb)
     &             -a6('sl',j3,j2,j1,j4,j5,j6,za,zb))/xnsq
      else
      write(6,*) 'Unimplemented st in a61',st
      stop
      endif
      return
      end

