!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function ampqqbgll(p1,p2,p3,p4,p5,al,be,ga,de,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
c     This routine calculates the amplitude for a right-handed quark
c     0--> q^+(1)+qb^1(2)+g^+(3)+l^-(4)+lb^+(5)
c     according to Eq.22 of 1309.3245v3
      integer:: p1,p2,p3,p4,p5
      complex(dp)::ampqqbgll
      complex(dp)::al,be,ga,de
      ampqqbgll=
     & +al*za(p1,p2)*zb(p1,p5)*za(p4,p2)/(za(p1,p3)*za(p3,p2))
     & +be*zb(p3,p5)*za(p4,p2)/za(p1,p3)
     & +ga*zb(p3,p1)*zb(p3,p5)*za(p4,p1)/(za(p1,p3)*zb(p2,p3))
     & +de*zb(p3,p1)/(za(p1,p3)*zb(p2,p1))
     & *(zb(p1,p5)*za(p4,p1)+zb(p2,p5)*za(p4,p2)+zb(p3,p5)*za(p4,p3))

      return
      end
