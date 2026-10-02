!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFPMscT(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.21
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FFPMscT,zab2,zba,zab,zb22b
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t,Delta3,delta12,delta56
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
      zb22b(j1,j2,j3,j4,j5,j6)=
     & zb(j1,j2)*zab2(j2,j4,j5,j6)+zb(j1,j3)*zab2(j3,j4,j5,j6)
c     End statement functions

      Delta3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2*s(j1,j2)*s(j3,j4)-2*s(j3,j4)*s(j5,j6)-2*s(j5,j6)*s(j1,j2)
      delta12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      delta56=s(j5,j6)-s(j1,j2)-s(j3,j4)

      FFPMscT=
     & 3*s(j1,j2)*zb(j4,j6)*za(j5,j6)*zab2(j4,j1,j2,j3)
     &  /(zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4)*Delta3**2)
     & *(za(j1,j2)*zb(j5,j6)
     & *(zab(j5,j6,j1)*delta56-zab(j5,j2,j1)*delta12)
     & +t(j1,j2,j4)*(zab(j2,j5,j6)*delta56-zab(j2,j1,j6)*delta12))

     & +za(j1,j2)*zb(j4,j6)
     & /(zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4)*Delta3)
     & *(-zb(j1,j2)*za(j5,j6)
     & *(3*zab2(j4,j1,j2,j3)*zab2(j2,j1,j3,j6)
     & -za(j2,j4)*zb(j3,j6)*(t(j1,j2,j3)-t(j1,j2,j4)))

     & +zab2(j4,j1,j2,j3)
     & *(zab(j5,j2,j1)*delta12-zab(j5,j6,j1)*delta56
     & +t(j1,j2,j3)*zab2(j5,j2,j4,j1)-t(j1,j2,j4)*zab2(j5,j2,j3,j1))

     & +za(j4,j5)*t(j1,j2,j3)*(zb22b(j1,j3,j4,j1,j2,j3)
     & -zba(j1,j2,j4)*zb(j4,j3))
     & -zab2(j4,j1,j2,j3)*za(j3,j5)*t(j1,j2,j3)
     & *(zb22b(j1,j3,j4,j1,j2,j4)-zba(j1,j2,j3)*zb(j3,j4))
     & /zab2(j3,j1,j2,j4))

     & -zab(j2,j1,j3)*zab(j5,j4,j6)
     & /(zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4))

     & +3*zab(j1,j5,j6)*zab2(j4,j1,j2,j3)/(zab2(j1,j2,j3,j4)*Delta3**2)
     & *(zab(j2,j1,j3)*(zab(j5,j2,j1)*delta12-zab(j5,j6,j1)*delta56)
     & +zab(j1,j2,j3)*(zab(j5,j6,j2)*delta56-zab(j5,j1,j2)*delta12)
     & *zab2(j2,j3,j4,j1)/zab2(j1,j3,j4,j2))

     & +za(j1,j5)/(zab2(j1,j2,j3,j4)*Delta3)
     & *(-zb(j1,j3)
     & *(zab(j2,j1,j3)*zab(j4,j5,j6)+zab(j4,j2,j3)*zab(j2,j5,j6)
     & -zab2(j4,j1,j2,j3)*(zab(j2,j1,j6)-zab(j2,j5,j6)))

     & +zb(j2,j3)*zab2(j2,j3,j4,j1)/zab2(j1,j3,j4,j2)
     & *(3*zab(j1,j5,j6)*zab(j4,j2,j3)-zab(j1,j2,j6)*zab2(j4,j1,j2,j3)
     & +zb(j2,j3)*za(j1,j4)*zab(j1,j5,j6)*(t(j2,j3,j4)-t(j1,j3,j4))
     & /zab2(j1,j3,j4,j2)))

      return
      end
