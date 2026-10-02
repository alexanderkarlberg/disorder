!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function t(j1,j2,j3)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      real(dp):: t
      integer:: j1,j2,j3
      t=s(j1,j2)+s(j2,j3)+s(j1,j3)
      end
