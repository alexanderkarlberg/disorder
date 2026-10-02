!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine checkndotp(p,n,ig)
      implicit none
      include 'types.f'
c--- routine to check that the vector "n" contracted with the LO
c--- matrix elements, is properly defined such that n.p(ig)=0
      include 'mxpart.f'
      real(dp):: p(mxpart,4),n(4),test,tolerance
      integer:: ig
      parameter (tolerance=1.e-2_dp)
      logical:: failedndp
      common/failedndp/failedndp
!$omp threadprivate(/failedndp/)

      failedndp=.false.
      test=n(4)*p(ig,4)-n(1)*p(ig,1)-n(2)*p(ig,2)-n(3)*p(ig,3)
      test=test/max(abs(n(4)),abs(n(1)),abs(n(2)),abs(n(3)))
     &         /max(abs(p(ig,4)),abs(p(ig,1)),abs(p(ig,2)),abs(p(ig,3)))

      if (test > tolerance) then
        write(6,*) 'warning: tolerance for n.p check in gvec routine'
     & //' exceeded: ',tolerance,test
c        write(*,*) 'tolerance: ',tolerance
c        write(6,*) 'test: ',test
c        write(*,*) 'index of gluon: ',ig
c        write(6,*) 'n: ',n(1),n(2),n(3),n(4)
c        write(6,*) 'p: ',p(ig,1),p(ig,2),p(ig,3),p(ig,4)
c        write(6,*)
c        call writeout(p)
        call flush(6)
        failedndp=.true.
c        stop
      endif

      return
      end

