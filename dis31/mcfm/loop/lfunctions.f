!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     July, 1998/July 2005/June 2007                                   *
c***********************************************************************

      function L0(x,y)
      implicit none
      include 'types.f'
      complex(dp):: L0
      include 'constants.f'
      include 'cplx.h'
      complex(dp):: Lnrat
      real(dp):: x,y,denom
      denom=one-x/y
      if (abs(denom) < 1.e-7_dp) then
      L0=-cone-cplx1(denom*(half+denom/three))
      else
      L0=Lnrat(x,y)/cplx1(denom)
      endif
      return
      end

      function L1(x,y)
      implicit none
      include 'types.f'
      complex(dp):: L1
      include 'constants.f'
      include 'cplx.h'
      real(dp):: x,y,denom
      complex(dp):: L0
      denom=one-x/y
      if (abs(denom) < 1.e-7_dp) then
      L1=-half*cone-cplx1(denom/three*(one+0.75_dp*denom))
      else
      L1=(L0(x,y)+cone)/cplx1(denom)
      endif
      return
      end

      function L2(x,y)
      implicit none
      include 'types.f'
      complex(dp):: L2
      include 'constants.f'
      include 'cplx.h'
      complex(dp):: Lnrat
      real(dp):: x,y,r,denom
      r=x/y
      denom=one-r
      if (abs(denom) < 1.e-7_dp) then
      L2=(cplx1(10._dp)+denom*(cplx1(15._dp)+cplx1(18._dp)*denom))
     & /cplx1(60._dp)
      else
      L2=(Lnrat(x,y)-cplx1(half*(r-one/r)))/cplx1(denom)**3
      endif
      return
      end

      function L0old(x,y)
      implicit none
      include 'types.f'
      complex(dp):: L0old
      include 'constants.f'
      include 'cplx.h'
      complex(dp):: Lnrat
      real(dp):: x,y,denom
      denom=one-x/y
      L0old=Lnrat(x,y)/cplx1(denom)
      return
      end


      function L1old(x,y)
      implicit none
      include 'types.f'
      complex(dp):: L1old
      include 'constants.f'
      include 'cplx.h'
      real(dp):: x,y,denom
      complex(dp):: L0old
      denom=one-x/y
      L1old=(L0old(x,y)+cone)/cplx1(denom)
      return
      end

      function Ls0(x1,y1,x2,y2)
      implicit none
      include 'types.f'
      complex(dp):: Ls0
      include 'constants.f'
      include 'cplx.h'
      real(dp):: x1,x2,y1,y2,r1,r2
      complex(dp):: Lsm1
      r1=x1/y1
      r2=x2/y2
      Ls0=Lsm1(x1,y1,x2,y2)/cplx1(one-r1-r2)
      return
      end

      function Ls1(x1,y1,x2,y2)
      implicit none
      include 'types.f'
      complex(dp):: Ls1
      include 'constants.f'
      include 'cplx.h'
      real(dp):: x1,x2,y1,y2,r1,r2
      complex(dp):: Ls0,L0
      r1=x1/y1
      r2=x2/y2
      Ls1=(Ls0(x1,y1,x2,y2)+L0(x1,y1)+L0(x2,y2))/cplx1(one-r1-r2)
      return
      end



      function Lsm1(x1,y1,x2,y2)
      implicit none
      include 'types.f'
      complex(dp):: Lsm1

      include 'constants.f'
      include 'cplx.h'
      real(dp):: x1,x2,y1,y2,r1,r2,omr1,omr2,ddilog
      complex(dp):: dilog1,dilog2,Lnrat
      r1=x1/y1
      r2=x2/y2
      omr1=one-r1
      omr2=one-r2
      if (omr1 > one) then
         dilog1=cplx1(pisqo6-ddilog(r1))-Lnrat(x1,y1)*cplx1(log(omr1))
      else
         dilog1=cplx1(ddilog(omr1))
      endif
      if (omr2 > one) then
         dilog2=cplx1(pisqo6-ddilog(r2))-Lnrat(x2,y2)*cplx1(log(omr2))
      else
          dilog2=cplx1(ddilog(omr2))
      endif
      lsm1=dilog1+dilog2+Lnrat(x1,y1)*Lnrat(x2,y2)-cplx1(pisqo6)
      return
      end

      function Lsm1_2mh(s,t,m1sq,m2sq)
      implicit none
      include 'types.f'
      complex(dp):: Lsm1_2mh
      include 'constants.f'
      real(dp):: s,t,m1sq,m2sq
      complex(dp):: lsm1_2mht,I3m
      Lsm1_2mh=Lsm1_2mht(s,t,m1sq,m2sq)
     & +(half*(s-m1sq-m2sq)+m1sq*m2sq/t)*I3m(s,m1sq,m2sq)
      return
      end

      function Lsm1_2mht(s,t,m1sq,m2sq)
      implicit none
      include 'types.f'
      complex(dp):: Lsm1_2mht
      include 'constants.f'
      include 'cplx.h'
      real(dp):: s,t,m1sq,m2sq,ddilog,r1,r2,omr1,omr2
      complex(dp):: Lnrat,dilog1,dilog2
      r1=m1sq/t
      r2=m2sq/t
      omr1=one-r1
      omr2=one-r2

      if (omr1 > one) then
      dilog1=cplx1(pisqo6-ddilog(r1))-Lnrat(-m1sq,-t)*cplx1(log(omr1))
      else
      dilog1=cplx1(ddilog(omr1))
      endif
      if (omr2 > one) then
      dilog2=cplx1(pisqo6-ddilog(r2))-Lnrat(-m2sq,-t)*cplx1(log(omr2))
      else
      dilog2=cplx1(ddilog(omr2))
      endif
      lsm1_2mht=-dilog1-dilog2
     & +half*(Lnrat(-s,-m1sq)*Lnrat(-s,-m2sq)-Lnrat(-s,-t)**2)
      return
      end


      function Lsm1_2me(s,t,m1sq,m3sq)
      implicit none
      include 'types.f'
      complex(dp):: Lsm1_2me

c---- formula taken from G.~Duplancic and B~Nizic [arXiv:hep-ph/0006249 v2]
c---- Eq 71
c---- Lsm1_2me notation follows from
c----  Z.~Bern, L.~J.~Dixon and D.~A.~Kosower,
c----  %``Dimensionally regulated pentagon integrals,''
c----  Nucl.\ Phys.\ B {\bf 412}, 751 (1994)
c----  [arXiv:hep-ph/9306240].
c----  %%CITATION = HEP-PH 9306240;%%
c----  Eqs. (I.13)
c---- analytic continuation has been checked by calculating numerically.
      include 'constants.f'
      include 'cplx.h'
      integer:: j
      real(dp):: s,t,m1sq,m3sq,ddilog,arg(4),omarg(4),f2me,htheta
      complex(dp):: Li2(4),wlog(4)
c--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(s)=half+half*sign(one,s)

      f2me=(s+t-m1sq-m3sq)/(s*t-m1sq*m3sq)

      arg(1)=f2me*s
      arg(2)=f2me*t
      arg(3)=f2me*m1sq
      arg(4)=f2me*m3sq

      do j=1,4
         omarg(j)=one-arg(j)
         wlog(j)=log(abs(arg(j)))
     &     +impi*cplx1(htheta(-arg(j))*sign(one,f2me))
         if (omarg(j) > one) then
             Li2(j)=cplx1(pisqo6-ddilog(arg(j)))
     &       -wlog(j)*cplx1(log(omarg(j)))
          else
             Li2(j)=cplx1(ddilog(omarg(j)))
         endif
      enddo
      Lsm1_2me=Li2(1)+Li2(2)-Li2(3)-Li2(4)

      return
      end

c      function Lsm1_2me(s,t,m1sq,m3sq)
c      implicit none
c      include 'types.f'
c      complex(dp):: Lsm1_2me

c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      integer:: j
c      real(dp):: s,t,m1sq,m3sq,ddilog,arg(5),omarg(5)
c      complex(dp):: Lnrat,Li2(5),wlog(5)
c      arg(1)=m1sq/s
c      wlog(1)=Lnrat(-m1sq,-s)

c      arg(2)=m1sq/t
c      wlog(2)=Lnrat(-m1sq,-t)

c      arg(3)=m3sq/s
c      wlog(3)=Lnrat(-m3sq,-s)

c      arg(4)=m3sq/t
c      wlog(4)=Lnrat(-m3sq,-t)

c      arg(5)=arg(1)*arg(4)
c      wlog(5)=Lnrat(-m1sq,-s)+Lnrat(-m3sq,-t)

c      do j=1,5
c         omarg(j)=one-arg(j)
c         if (omarg(j) > one) then
c             Li2(j)=cplx1(pisqo6-ddilog(arg(j)))-wlog(j)*log(omarg(j))
c          else
c             Li2(j)=cplx1(ddilog(omarg(j)))
c         endif
c      enddo
c      Lsm1_2me=Li2(5)-Li2(1)-Li2(2)-Li2(3)-Li2(4)-half*Lnrat(-s,-t)**2
c      return
c      end

