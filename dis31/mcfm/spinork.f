!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine spinork(N,p,zabikj,zbaikj,k)
      implicit none
      include 'types.f'
c---Calculate spinor products dotted in with a vector k
c---extended to deal with negative energies ie with all momenta outgoing
c   zabikj=<i-|k|j-> zbaikj=<i+|k|j+>
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,
c---za(i,j)*zb(j,i)=s(i,j)
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),rt(mxpart),k(4),kp,km,flip(mxpart)
      complex(dp):: pr(mxpart),pl(mxpart),f(mxpart),kr,kl,
     & zabikj(mxpart,mxpart),zbaikj(mxpart,mxpart)
      integer:: i,j,N

c--setup components for vector which is contracted in
      kp=+k(4)+k(1)
      km=+k(4)-k(1)
      kr=cplx2(+k(3),-k(2))
      kl=cplx2(+k(3),+k(2))

c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
            zabikj(j,j)=czip
            zbaikj(j,j)=czip

c-----positive energy case
            if (p(j,4) > 0._dp) then
            flip(j)=1._dp
            f(j)=cone
            else
            flip(j)=-1._dp
            f(j)=im
            endif
            rt(j)=flip(j)*(p(j,4)+p(j,1))
!-----  for numerical safety
            if (rt(j) > 0._dp) then
              rt(j)=sqrt(rt(j))
            else
              rt(j)=1.e-32_dp
            endif
            pr(j)=cplx2(flip(j)*p(j,3),-flip(j)*p(j,2))
            pl(j)=conjg(pr(j))
      enddo
      do i=1,N
         do j=1,i
         zabikj(i,j)=f(i)*f(j)
     & *(pr(i)*pl(j)*cplx1(kp/(rt(i)*rt(j)))
     &    -pr(i)*kl*cplx1(rt(j)/rt(i))
     &    -cplx1(rt(i)/rt(j))*kr*pl(j)+cplx1(rt(i)*rt(j)*km))
         zbaikj(j,i)=zabikj(i,j)
         zabikj(j,i)=flip(i)*flip(j)*conjg(zabikj(i,j))
         zbaikj(i,j)=flip(i)*flip(j)*conjg(zbaikj(j,i))

      enddo
      enddo

      return
      end
