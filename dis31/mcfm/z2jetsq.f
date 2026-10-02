!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine z2jetsq(i1,i2,i3,i4,i5,i6,za,zb,msq)
      implicit none
      include 'types.f'
c-----Author R.K. Ellis
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> Gamma^* +g(p5) +g(p6)
c                          |
c                          --> l(p3)+a(p4)

c--all momenta incoming
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'mmsq_cs.f'
      complex(dp):: qcd1LL(-1:1,-1:1),qcd2LL(-1:1,-1:1)
      complex(dp):: qcd1LR(-1:1,-1:1),qcd2LR(-1:1,-1:1)
c      complex(dp):: qcd1RL(-1:1,-1:1),qcd2RL(-1:1,-1:1)
c      complex(dp):: qcd1RR(-1:1,-1:1),qcd2RR(-1:1,-1:1)
      complex(dp):: qedLL(-1:1,-1:1),qedLR(-1:1,-1:1)
c      complex(dp):: qedRL(-1:1,-1:1),qedRR(-1:1,-1:1)
      real(dp):: msq1(2,2),msq2(2,2),msqq(2,2),msq(2,2)
      integer:: i1,i2,i3,i4,i5,i6,j,k,pq,pl
      integer,parameter::pol(2)=(/-1,1/)

      call subqcd(i1,i2,i3,i4,i5,i6,za,zb,qcd1LL)
      call subqcd(i1,i2,i3,i4,i6,i5,za,zb,qcd2LL)


c      call subqcd(i2,i1,i4,i3,i5,i6,za,zb,qcd1RR)
c      call subqcd(i2,i1,i4,i3,i6,i5,za,zb,qcd2RR)

      call subqcd(i1,i2,i4,i3,i5,i6,za,zb,qcd1LR)
      call subqcd(i1,i2,i4,i3,i6,i5,za,zb,qcd2LR)

c      call subqcd(i2,i1,i3,i4,i5,i6,za,zb,qcd1RL)
c      call subqcd(i2,i1,i3,i4,i6,i5,za,zb,qcd2RL)


      do j=1,2
      do k=1,2
      qedLL(pol(j),pol(k))=qcd1LL(pol(j),pol(k))+qcd2LL(pol(k),pol(j))
      qedLR(pol(j),pol(k))=qcd1LR(pol(j),pol(k))+qcd2LR(pol(k),pol(j))
c      qedRL(pol(j),pol(k))=qcd1RL(pol(j),pol(k))+qcd2RL(pol(k),pol(j))
c      qedRR(pol(j),pol(k))=qcd1RR(pol(j),pol(k))+qcd2RR(pol(k),pol(j))
      enddo
      enddo

      do pq=1,1
      do pl=1,2
      msq1(pq,pl)=0._dp
      msq2(pq,pl)=0._dp
      msqq(pq,pl)=0._dp
      enddo
      enddo

c---sum over gluon polarizations
      do j=1,2
      do k=1,2
      msq1(1,1)=msq1(1,1)+real(qcd1LL(pol(j),pol(k))
     &                  *conjg(qcd1LL(pol(j),pol(k))),dp)
      msq2(1,1)=msq2(1,1)+real(qcd2LL(pol(j),pol(k))
     &                  *conjg(qcd2LL(pol(j),pol(k))),dp)
      msqq(1,1)=msqq(1,1)+real(qedLL(pol(j),pol(k))
     &                  *conjg(qedLL(pol(j),pol(k))),dp)

      msq1(1,2)=msq1(1,2)+real(qcd1LR(pol(j),pol(k))
     &                  *conjg(qcd1LR(pol(j),pol(k))),dp)
      msq2(1,2)=msq2(1,2)+real(qcd2LR(pol(j),pol(k))
     &                  *conjg(qcd2LR(pol(j),pol(k))),dp)
      msqq(1,2)=msqq(1,2)+real(qedLR(pol(j),pol(k))
     &                  *conjg(qedLR(pol(j),pol(k))),dp)

c      msq1(2,1)=msq1(2,1)+abs(qcd1RL(pol(j),pol(k)))**2
c      msq2(2,1)=msq2(2,1)+abs(qcd2RL(pol(j),pol(k)))**2
c      msqq(2,1)=msqq(2,1)+abs(qedRL(pol(j),pol(k)))**2

c      msq1(2,2)=msq1(2,2)+abs(qcd1RR(pol(j),pol(k)))**2
c      msq2(2,2)=msq2(2,2)+abs(qcd2RR(pol(j),pol(k)))**2
c      msqq(2,2)=msqq(2,2)+abs(qedRR(pol(j),pol(k)))**2

      enddo
      enddo


      msq1(2,2)=msq1(1,1)
      msq1(2,1)=msq1(1,2)

      msq2(2,2)=msq2(1,1)
      msq2(2,1)=msq2(1,2)

      msqq(2,2)=msqq(1,1)
      msqq(2,1)=msqq(1,2)

      do pq=1,2
      do pl=1,2
      mmsq_cs(0,pq,pl)=-ninth*msqq(pq,pl)
      mmsq_cs(1,pq,pl)=msq1(pq,pl)
      mmsq_cs(2,pq,pl)=msq2(pq,pl)
      msq(pq,pl)=mmsq_cs(2,pq,pl)+mmsq_cs(1,pq,pl)+mmsq_cs(0,pq,pl)
      enddo
      enddo
      return
      end
