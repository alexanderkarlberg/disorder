!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c      function aqqb_zbb(i1,i2,i3,i4,i5,i6)
c      implicit none
c      include 'types.f'
c      complex(dp):: aqqb_zbb
c--- Note that this is the amplitude for particle labels
c--- q1, qb2, Q5, Qb6, l3, lb4
c--- This corresponds to A++(1,6,5,2) of eq. (12.3) in BDK
c--- and we note that    A++(1,6,5,2) = A++(1,4,3,2)

c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'sprods_com.f'
c      include 'zprods_com.f'
c      integer:: i1,i2,i3,i4,i5,i6
c      complex(dp):: t2a
c      real(dp):: s234,s256,prop
c--- statement function
c      t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

c      s234=s(i2,i3)+s(i2,i4)+s(i3,i4)
c      s256=s(i2,i6)+s(i2,i5)+s(i5,i6)
c      prop=s(i5,i6)*s(i3,i4)
c      aqqb_zbb=
c     & +zb(i1,i4)*za(i5,i2)*t2a(i3,i1,i4,i6)/(prop*s256)
c     & +za(i3,i2)*zb(i6,i1)*t2a(i5,i2,i3,i4)/(prop*s234)

c      return
c      end

      function aqqb_zbb_new(i1,i2,i3,i4,i5,i6)
      implicit none
      include 'types.f'
      complex(dp):: aqqb_zbb_new
c--- This corresponds to A++(1,2,3,4) of eq. (12.3) in BDK
c    The notation of BDK calculates the following amplitude

c     q3(L)----<----------q2            q3(L)------<--------q2
c                 0                             0
c                 0                             0
c                 0                             0
c     q1(R)------<--------q4            q1(R)------<--------q4
c             )                                         )
c            (                                         (
c             )                                         )
c     l5(L)-------<-------l6            l5(L)-------<-------l6

c     Note that this function has the property
c     Conjg(aqqb_zbb_new(i1,i2,i3,i4,i5,i6))=
c          -aqqb_zbb_new(i4,i3,i2,i1,i6,i5)
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: zab2
      real(dp):: s123,s234,prop
c--- statement function
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      s123=s(i1,i2)+s(i1,i3)+s(i2,i3)
      s234=s(i2,i3)+s(i2,i4)+s(i3,i4)
      prop=s(i2,i3)*s(i5,i6)

      aqqb_zbb_new=
     & +zb(i1,i2)*za(i5,i4)*zab2(i3,i1,i2,i6)/(prop*s123)
     & +za(i3,i4)*zb(i6,i1)*zab2(i5,i3,i4,i2)/(prop*s234)

      return
      end
