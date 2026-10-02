!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function vv(st,j1,j2,j3,j4,j5,j6)
      implicit none
      include 'types.f'
      complex(dp):: vv
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     July, 1998.                                                      *
c***********************************************************************
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: Lnrat,xl12,xl34,xl23,xl56
      character(len=2):: st


      xl12=Lnrat(musq,-s(j1,j2))
      xl34=Lnrat(musq,-s(j3,j4))

      if ((st == 'pp') .or. (st == 'pm')) then

c   Vpm =Vpp=
c     - 2*epinv^2   + epinv* ( 2/3 - xlog12 - xlog34 )
c     + 10/9 - 1/2*xlog12^2 - 1/2*xlog34^2 + 13/6*xlog23 - 3/2*xlog56
c   Vsl =
c    -2*epinv^2 + epinv* ( - 3 - xlog12 - xlog34 )
c    -15/2 - 3/2*xlog12 - 1/2*xlog12^2 - 3/2*xlog34 - 1/2*xlog34^2

      xl23=Lnrat(musq,-s(j2,j3))
      xl56=Lnrat(musq,-s(j5,j6))
      vv =-epinv*(two*epinv2-two/three+xl12+xl34)
     & +10._dp/9._dp-0.5_dp*(xl12**2+xl34**2)+13._dp/6._dp*xl23-1.5_dp*xl56


      elseif (st == 'sl') then
      vv=
     & -epinv*(two*epinv2+three+xl12+xl34)
     & -7.5_dp-1.5_dp*(xl12+xl34)-0.5_dp*(xl12**2+xl34**2)


      endif

      end




