!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function lnrat(x,y)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: lnrat
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     August, 1998.                                                    *
c     lnrat(x,y)=log(x-i*ep)-log(y-i*ep)                               *
c     this function is hard-wired for sign of epsilon we must adjust   *
c     sign of x and y to get the right sign for epsilon                *
c***********************************************************************
      include 'constants.f'
       real(dp):: x,y,htheta
c--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+half*sign(one,x)
c---  end statement function
      lnrat=cplx1(log(abs(x/y)))
     & -impi*(htheta(-x)-htheta(-y))
      return
      end

