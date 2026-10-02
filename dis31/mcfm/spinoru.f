!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine spinoru(N,p,za,zb)
      implicit none
      include 'types.f'
c---Calculate spinor products
c---extended to deal with negative energies ie with all momenta outgoing
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,
c---za(i,j)*zb(j,i)=s(i,j)
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'

      integer, intent(in) :: N
      real(dp), intent(in) :: p(mxpart,4)
      complex(dp), intent(out) :: za(mxpart,mxpart),zb(mxpart,mxpart)

      real(dp):: rt(mxpart)
      complex(dp):: c23(mxpart),f(mxpart)
      integer:: i,j

c---if one of the vectors happens to be zero this routine fails.
      do j=1,N
         za(j,j)=czip
         zb(j,j)=za(j,j)

c-----positive energy case
         if (p(j,4) > zip) then
            rt(j)=p(j,4)+p(j,1)
!-----  for numerical safety
            if (rt(j) > 0._dp) then
              rt(j)=sqrt(rt(j))
            else
              rt(j)=1.e-32_dp
            endif
            c23(j)=cmplx(p(j,3),-p(j,2),kind=dp)
            f(j)=cone
         else
c-----negative energy case
            rt(j)=-p(j,4)-p(j,1)
!-----  for numerical safety
            if (rt(j) > 0._dp) then
              rt(j)=sqrt(rt(j))
            else
              rt(j)=1.e-32_dp
            endif
            c23(j)=cmplx(-p(j,3),p(j,2),kind=dp)
            f(j)=im
         endif
      enddo
      do i=2,N
         do j=1,i-1
c disorder (2026-10-02): Minkowski product inline instead of MCFM's
c function dot, whose name clashes with DISENT's DOT
         s(i,j)=2*(p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)
     &            -p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)
     &   *(c23(i)*rt(j)/rt(i)-c23(j)*rt(i)/rt(j))
! Check for numerical safety -- this routine is sometimes called with
! exactly collinear momenta, e.g. auxiliary vectors for massive vectors,
! and we do not want it to create a FPE
         if (abs(za(i,j)) == 0._dp) then
           zb(i,j)=0._dp
         else

c Disabled this check August 2018 since it can cause poor cancellation
c Reinstated with a much tighter check for 4f single top calculation
c         if (abs(s(i,j))<1.e-10_dp) then
c           zb(i,j)=-(f(i)*f(j))**2*conjg(za(i,j))
c         else
           zb(i,j)=-s(i,j)/za(i,j)
c         endif
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

      return
      end
