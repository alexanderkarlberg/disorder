!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Master1(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq 7.1
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t,delta12,delta56,Delta3
      complex(dp):: Master1,zab2,zba2,zab,zba,I3m,Lnrat,M1bit1
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
c     End statement functions
      delta12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      delta56=s(j5,j6)-s(j1,j2)-s(j3,j4)
      Delta3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2*s(j1,j2)*s(j3,j4)-2*s(j3,j4)*s(j5,j6)-2*s(j5,j6)*s(j1,j2)
      M1bit1=
     & +2/(zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4)*Delta3)
     & *(zb(j1,j3)/za(j5,j6)
     & *(za(j2,j1)*zba2(j1,j3,j4,j5)
     & -(zab2(j2,j3,j4,j1)*za(j1,j5)+zab2(j2,j3,j4,j2)*za(j2,j5)))
     & *(za(j1,j5)*t(j1,j2,j4)+zab(j1,j2,j6)*za(j6,j5))
     & -zb(j1,j6)*za(j2,j4)/(za(j3,j4)*zb(j5,j6))
     & *((s(j1,j4)-s(j2,j3))
     & *(zab(j1,j2,j6)*delta12-zab(j1,j5,j6)*delta56)
     & +zab2(j1,j3,j4,j2)*(zab(j2,j1,j6)*delta12-zab(j2,j5,j6)*delta56))
     & +2*zab(j2,j1,j3)*zab2(j5,j1,j2,j6)*(s(j1,j4)-s(j2,j3)))

      Master1=zb(j1,j3)
     & /(zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4)*Delta3)
     & *(2*za(j1,j2)*zab(j5,j2,j6)
     & *(t(j1,j2,j3)*delta12+s(j5,j6)*delta56)
     & +2*za(j1,j2)*za(j1,j5)*zab(j4,j5,j6)
     & *(zb(j1,j4)*delta56-2*zba(j1,j2,j3)*zb(j3,j4))
     & +s(j5,j6)*za(j2,j4)
     & *(za(j1,j5)*zb(j4,j6)*delta56-2*zab(j1,j2,j6)*zab(j5,j3,j4)))
     & *I3m(s(j1,j2),s(j3,j4),s(j5,j6))

     & +M1bit1*lnrat(-s(j1,j2),-s(j5,j6))
      return
      end
