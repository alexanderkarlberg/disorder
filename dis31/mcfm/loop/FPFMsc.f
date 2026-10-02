!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPFMsc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.14
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FPFMsc,Master2,Master3,zab2,zab
      complex(dp):: FPFMsc1,FPFMsc2,FPFMsc3,FPFMsc4,FPFMsc5
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lnrat
      real(dp):: t,s234,s23,s56
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

c      FPFMsc1=
c     & +0.5_dp*zab(j5,j1,j2)/(zb(j3,j4)*zab2(j1,j2,j3,j4))
c     & *zab2(j5,j3,j4,j2)/za(j5,j6)

c      FPFMsc2=
c     & +0.5_dp*zab(j5,j1,j2)/(zb(j3,j4)*zab2(j1,j2,j3,j4))
c     & *zb21b(j2,j3,j4,j1,j6)

c      FPFMsc3=
c     & +0.5_dp*zb(j2,j3)*za(j3,j4)*zab2(j5,j2,j3,j4)
c     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4))
c     & *zab2(j5,j3,j4,j2)/t(j2,j3,j4)

c      FPFMsc4=
c     & +0.5_dp*zb(j2,j3)*za(j3,j4)*zab2(j5,j2,j3,j4)
c     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4))
c     & *zab(j5,j4,j2)

c      FPFMsc5=
c     & -0.5_dp*zb(j1,j2)/(zb(j1,j4)*zab2(j1,j2,j3,j4)*zab2(j2,j1,j4,j3))
c     & *(zab2(j5,j3,j4,j1)*za(j1,j5)/za(j5,j6)
c     & -zb21b(j6,j1,j2,j4,j6)/zb(j5,j6))

      s234=t(j2,j3,j4)
      s23=s(j2,j3)
      s56=s(j5,j6)
      FPFMsc=za(j1,j3)/za(j2,j3)*Master2(j3,j2,j1,j4,j6,j5,zb,za)
     &      +za(j1,j4)/za(j2,j4)*Master3(j3,j2,j1,j4,j6,j5,zb,za)

     & +FPFMsc1(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s56)/s56
     & +FPFMsc2(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s56,-s234)/s234**2
     & +FPFMsc3(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s23)/s23
     & +FPFMsc4(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s234,-s23)/s23**2
     & +FPFMsc5(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s23,-s56)
c       fpfm_poly
     & +0.5_dp*(
     & -zb(j1,j2)*za(j1,j5)*za(j3,j5)
     & /(za(j2,j3)*zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4))
     & -zb(j1,j3)**2*za(j3,j5)**2
     & /(zb(j1,j4)*za(j2,j3)*zb(j3,j4)*za(j5,j6)*zab2(j2,j1,j4,j3))
     & +zab(j4,j3,j1)*za(j1,j5)*za(j3,j5)
     & /(za(j2,j3)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j2,j1,j4,j3))

     & -za(j1,j3)*za(j3,j5)*za(j4,j5)
     & /(za(j1,j2)*za(j2,j3)*za(j5,j6)*zab2(j1,j2,j3,j4))
     & -s(j3,j4)*za(j1,j5)*za(j4,j5)
     & /(za(j1,j2)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j2,j1,j4,j3))
     & -zb(j1,j2)*za(j1,j5)*zab2(j5,j3,j4,j2)
     & /(zb(j3,j4)*za(j5,j6)*s234*zab2(j1,j2,j3,j4))

     & -za(j1,j3)*zb(j1,j6)**2
     & /(zb(j1,j4)*za(j1,j2)*za(j2,j3)*zb(j3,j4)*zb(j5,j6))
     & +zab(j3,j4,j6)*zab2(j1,j2,j3,j6)
     & /(za(j1,j2)*za(j2,j3)*zb(j3,j4)*zb(j5,j6)*zab2(j1,j2,j3,j4))
     & -s(j1,j3)*zab(j3,j4,j6)*zb(j3,j6)
     & /(za(j2,j3)*zb(j3,j4)*zb(j5,j6)
     & *zab2(j1,j2,j3,j4)*zab2(j2,j1,j4,j3))

     & +za(j1,j4)*zab(j3,j4,j6)*zb(j3,j6)
     & /(za(j1,j2)*zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j2,j1,j4,j3))
     & +zab(j3,j1,j6)*zb(j1,j2)*zb(j4,j6)
     & /(zb(j1,j4)*za(j2,j3)*zb(j3,j4)*zb(j5,j6)*zab2(j1,j2,j3,j4))
     & -zab2(j4,j1,j2,j6)*zb(j2,j6)*za(j2,j4)
     & /(za(j1,j2)*zb(j5,j6)*t(j1,j2,j4)*zab2(j2,j1,j4,j3)))

      return
      end


      function FPFMsc1(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.14
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMsc1,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FPFMsc1=
     & +0.5_dp*zab(j5,j1,j2)/(zb(j3,j4)*zab2(j1,j2,j3,j4))
     & *zab2(j5,j3,j4,j2)/za(j5,j6)

      return
      end


      function FPFMsc2(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.14
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMsc2,zab2,zba2,zab,zb21b
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zb21b(j1,j2,j3,j4,j5)=zba2(j1,j2,j3,j4)*zb(j4,j5)
c     End statement functions

      FPFMsc2=
     & +0.5_dp*zab(j5,j1,j2)/(zb(j3,j4)*zab2(j1,j2,j3,j4))
     & *zb21b(j2,j3,j4,j1,j6)

      return
      end


      function FPFMsc3(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.14
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMsc3,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp)::t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FPFMsc3=
     & +0.5_dp*zb(j2,j3)*za(j3,j4)*zab2(j5,j2,j3,j4)
     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4))
     & *zab2(j5,j3,j4,j2)/t(j2,j3,j4)

      return
      end


      function FPFMsc4(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.14
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMsc4,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FPFMsc4=
     & +0.5_dp*zb(j2,j3)*za(j3,j4)*zab2(j5,j2,j3,j4)
     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4))
     & *zab(j5,j4,j2)

      return
      end


      function FPFMsc5(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.14
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMsc5,zab2,zba2,zb21b
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zb21b(j1,j2,j3,j4,j5)=zba2(j1,j2,j3,j4)*zb(j4,j5)
c     End statement functions

      FPFMsc5=
     & -0.5_dp*zb(j1,j2)/(zb(j1,j4)*zab2(j1,j2,j3,j4)*zab2(j2,j1,j4,j3))
     & *(zab2(j5,j3,j4,j1)*za(j1,j5)/za(j5,j6)
     & -zb21b(j6,j1,j2,j4,j6)/zb(j5,j6))

      return
      end



