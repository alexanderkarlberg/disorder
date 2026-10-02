!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine storecsz(mcs)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the Z2jet matrix elements into separate arrays for each
c-- incoming parton case

      include 'mmsq_cs.f'
      integer:: i,j,k
      real(dp):: mcs(0:2,2,2)

      do i=0,2
        do j=1,2
          do k=1,2
        mcs(i,j,k)=mmsq_cs(i,j,k)
          enddo
        enddo
      enddo
      return
      end
