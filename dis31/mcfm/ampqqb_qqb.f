!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine ampqqb_qqb(i1,i2,i5,i6,qqbA,qqbB)
      implicit none
      include 'types.f'

      integer:: i1,i2,i5,i6,j,k
      complex(dp):: aqqb_zbb_new
      complex(dp):: qqbA(2,2,2),qqbB(2,2,2)
      integer,parameter::swap(2)=(/2,1/)
c--- also include diagrams where the Z is attached to b-bbar line
c--- notation: (qqb hel, bbb hel, outgoing lepton helicity is left-handed

c--- Z to qqb, L L L
      qqbA(1,1,1)=+aqqb_zbb_new(i1,i6,i5,i2,3,4)
c--- Z to bbb, L L
      qqbB(1,1,1)=+aqqb_zbb_new(i6,i1,i2,i5,3,4)

c--- Z to qqb, R R
      qqbA(2,2,1)=-aqqb_zbb_new(i2,i5,i6,i1,3,4)
c--- Z to bbb, R R
      qqbB(2,2,1)=-aqqb_zbb_new(i5,i2,i1,i6,3,4)

c--- Z to qqb, L R
      qqbA(1,2,1)=+aqqb_zbb_new(i1,i5,i6,i2,3,4)
c--- Z to bbb, L R
      qqbB(1,2,1)=-aqqb_zbb_new(i5,i1,i2,i6,3,4)

c--- Z to qqb, R L
      qqbA(2,1,1)=-aqqb_zbb_new(i2,i6,i5,i1,3,4)
c--- Z to bbb, R L
      qqbB(2,1,1)=+aqqb_zbb_new(i6,i2,i1,i5,3,4)

      do j=1,2
      do k=1,2
      qqbA(j,k,2)=conjg(qqbA(swap(j),swap(k),1))
      qqbB(j,k,2)=conjg(qqbB(swap(j),swap(k),1))
      enddo
      enddo
      return
      end


