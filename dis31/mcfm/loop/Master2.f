!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function Master2(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq.7.3
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: Master2,Master2a,zab2,zba2,zab,zba,zb12b,za12a
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lsm1_2mh,I3m,M2bit1,M2bit2,M2bit3
      real(dp):: t,delta12,delta34,delta56,Delta3
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
      zb12b(j1,j2,j3,j4,j5)=zb(j1,j2)*zab2(j2,j3,j4,j5)
      za12a(j1,j2,j3,j4,j5)=za(j1,j2)*zba2(j2,j3,j4,j5)
c     End statement functions

      delta12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      delta34=s(j3,j4)-s(j5,j6)-s(j1,j2)
      delta56=s(j5,j6)-s(j1,j2)-s(j3,j4)
      Delta3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2*s(j1,j2)*s(j3,j4)-2*s(j3,j4)*s(j5,j6)-2*s(j5,j6)*s(j1,j2)
      M2bit1=t(j1,j2,j3)*zb(j4,j6)*zab(j5,j4,j1)
     & /(zb(j2,j3)*zab2(j4,j1,j2,j3))
      M2bit2=
     & -zab(j3,j2,j1)*zab2(j5,j1,j2,j3)**2
     & /(zb(j2,j3)*za(j5,j6)*zab2(j4,j1,j2,j3)**2)
      M2bit3=
     & -zab(j5,j4,j1)*zab2(j5,j1,j2,j3)*t(j1,j2,j3)
     & /(zb(j2,j3)*za(j5,j6)*zab2(j4,j1,j2,j3)**2)

      Master2=Master2a(j1,j2,j3,j4,j5,j6,za,zb)

      Master2=Master2
     & +zab2(j4,j2,j3,j1)*zab2(j5,j1,j2,j3)**2
     & /(zb(j2,j3)*za(j5,j6)*zab2(j4,j1,j2,j3)**3)
     & *Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))

     & +zb(j1,j2)/(zb(j2,j3)*zab2(j4,j1,j2,j3))
     & *(3*za(j1,j2)*zb(j3,j4)*zab2(j3,j1,j2,j4)*zab2(j5,j1,j2,j6)
     & /Delta3**2
     & *(zab(j4,j2,j1)*delta12-zab(j4,j3,j1)*delta34)

     & -1._dp/Delta3
     & *(zab(j2,j1,j3)*zab2(j3,j1,j2,j4)*za(j4,j5)*zb(j3,j6)
     & /zab2(j4,j1,j2,j3)*(t(j1,j2,j4)-t(j1,j2,j3))
     & +za(j2,j3)*zb(j4,j6)
     & *(zab(j5,j4,j3)*delta34-zab(j5,j6,j3)*delta56)

     & +zab2(j3,j1,j2,j4)*( zab(j2,j4,j3)*zab2(j5,j1,j2,j6)
     & +3*zab(j2,j1,j3)*zab(j5,j4,j6))

     & -zab(j2,j1,j4)*za(j3,j5)*(zb12b(j6,j4,j1,j2,j3)
     & +2*zba(j6,j5,j4)*zb(j4,j3)
     & +zb(j3,j6)*delta12))
     & +za12a(j2,j3,j1,j2,j5)*zb(j4,j6)/t(j1,j2,j3))
     & *I3m(s(j1,j2),s(j3,j4),s(j5,j6))


     & +M2bit1*L1(-s(j5,j6),-t(j1,j2,j3))/t(j1,j2,j3)**2
     & +M2bit2*L0(-t(j1,j2,j3),-s(j1,j2))/s(j1,j2)
     & +M2bit3*L0(-t(j1,j2,j3),-s(j5,j6))/s(j5,j6)

     & -0.5_dp*(zab(j4,j2,j1)*delta12-zab(j4,j3,j1)*delta34)
     & /(zb(j2,j3)*za(j3,j4)*zab2(j4,j1,j2,j3)*Delta3)
     & *(za(j3,j4)*zb(j4,j6)**2/zb(j5,j6)
     &  +za(j3,j5)**2*zb(j3,j4)/za(j5,j6))

     & -za(j3,j5)*zb(j4,j6)
     & *(zb(j1,j3)*delta56-2*zba(j1,j2,j4)*zb(j4,j3))
     & /(zb(j2,j3)*zab2(j4,j1,j2,j3)*Delta3)

      return
      end
