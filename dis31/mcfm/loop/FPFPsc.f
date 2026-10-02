!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPFPsc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FPFPsc,zab2
      complex(dp):: FPFPsc1,FPFPsc2,FPFPsc3,FPFPsc4,FPFPsc5
      complex(dp):: FPFPsc6,FPFPsc7,FPFPsc8,FPFPsc9,FPFPsc10,FPFPsc11
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lsm1_2me,Lnrat,Ls1
      real(dp):: t,s123,s234,s23,s56
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

c      FPFPsc1=
c     & -0.5_dp*za(j1,j3)*za(j5,j3)**2
c     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))

c      FPFPsc2=
c     & -za(j1,j3)*za(j5,j1)**2*zb(j1,j2)
c     & /(za(j1,j2)*za(j5,j6)*za(j4,j1)**2)

c      FPFPsc3=
c     & +za(j3,j4)*za(j5,j1)**2*zb(j2,j4)
c     & /(za(j1,j2)*za(j5,j6)*za(j4,j1)**2)


c      FPFPsc4=
c     & -za(j1,j3)*za(j5,j1)*za(j5,j4)*zab2(j3,j1,j2,j4)
c     & /(za(j1,j2)*za(j2,j3)*za(j5,j6)*za(j4,j1)**2)

c      FPFPsc5=
c     & +za(j1,j3)*zb(j6,j4)*za(j5,j3)/(za(j1,j2)*za(j2,j3)*za(j4,j1))


c      FPFPsc6=
c     & +za(j1,j3)*za(j5,j1)**2*zab2(j3,j2,j4,j1)
c     & /(za(j1,j2)*za(j2,j3)*za(j5,j6)*za(j4,j1)**2)

c      FPFPsc7=
c     & -za(j1,j3)**2*za(j5,j3)*zb(j6,j1)
c     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))

c      FPFPsc8=
c     & -0.5_dp*za(j2,j3)*za(j5,j4)**2*zb(j4,j2)**2
c     & /(za(j1,j4)*za(j5,j6)*za(j2,j4))


c      FPFPsc9=
c     & -0.5_dp*za(j5,j2)**2*zb(j2,j4)**2*za(j4,j3)
c     & /(za(j1,j2)*za(j5,j6)*za(j2,j4))


c      FPFPsc10=
c     & -za(j1,j3)*za(j3,j4)*za(j5,j6)*zb(j6,j4)**2
c     & /(za(j1,j2)*za(j2,j3)*za(j4,j1))

c      FPFPsc11=
c     & -0.5_dp*za(j1,j3)**3*za(j5,j6)*zb(j6,j1)**2
c     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))

      s123=t(j1,j2,j3)
      s234=t(j2,j3,j4)
      s23=s(j2,j3)
      s56=s(j5,j6)

      FPFPsc=za(j3,j4)*za(j5,j2)**2*zb(j2,j4)**2
     & /(za(j1,j2)*za(j5,j6)*za(j2,j4))
     & *Ls1(-s23,-s234,-s(j3,j4),-s234)
     & /s234**2

     & -za(j1,j3)*za(j3,j4)*za(j5,j1)**2
     & /(za(j1,j2)*za(j2,j3)*za(j5,j6)*za(j4,j1)**3)
     & *Lsm1_2me(s123,s234,s23,s56)

     & +FPFPsc1(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s234,-s56)
     & +FPFPsc2(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s123,-s23)/s23
     & +FPFPsc3(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s23)/s23
     & +FPFPsc4(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s123,-s56)/s56
     & +FPFPsc5(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s123,-s56)/s56
     & +FPFPsc6(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s56)/s56
     & +FPFPsc7(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s56)/s56
     & +FPFPsc8(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s23,-s234)/s234**2
     & +FPFPsc9(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s(j3,j4),-s234)/s234**2
     & +FPFPsc10(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s123,-s56)/s56**2
     & +FPFPsc11(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s234,-s56)/s56**2

     & +0.5_dp*zb(j2,j4)*zab2(j1,j2,j4,j6)*zab2(j3,j2,j4,j6)
     & /(za(j1,j2)*za(j4,j1)*zb(j5,j6)*t(j1,j2,j4)*s234)

     & +0.5_dp*za(j2,j4)*zb(j2,j4)**2*za(j5,j1)*za(j5,j3)
     & /(za(j1,j2)*za(j4,j1)*za(j5,j6)*t(j1,j2,j4)*s234)

     & +0.5_dp*za(j1,j3)*za(j3,j5)**2
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))

      return
      end


      function FPFPsc1(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc1
      integer:: j1,j2,j3,j4,j5,j6

      FPFPsc1=
     & -0.5_dp*za(j1,j3)*za(j5,j3)**2
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))

      return
      end


      function FPFPsc2(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc2
      integer:: j1,j2,j3,j4,j5,j6

      FPFPsc2=
     & -za(j1,j3)*za(j5,j1)**2*zb(j1,j2)
     & /(za(j1,j2)*za(j5,j6)*za(j4,j1)**2)

      return
      end


      function FPFPsc3(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc3
      integer:: j1,j2,j3,j4,j5,j6

      FPFPsc3=
     & +za(j3,j4)*za(j5,j1)**2*zb(j2,j4)
     & /(za(j1,j2)*za(j5,j6)*za(j4,j1)**2)

      return
      end


      function FPFPsc4(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc4,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FPFPsc4=
     & -za(j1,j3)*za(j5,j1)*za(j5,j4)*zab2(j3,j1,j2,j4)
     & /(za(j1,j2)*za(j2,j3)*za(j5,j6)*za(j4,j1)**2)

      return
      end


      function FPFPsc5(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc5
      integer:: j1,j2,j3,j4,j5,j6

      FPFPsc5=
     & +za(j1,j3)*zb(j6,j4)*za(j5,j3)/(za(j1,j2)*za(j2,j3)*za(j4,j1))

      return
      end


      function FPFPsc6(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc6,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FPFPsc6=
     & +za(j1,j3)*za(j5,j1)**2*zab2(j3,j2,j4,j1)
     & /(za(j1,j2)*za(j2,j3)*za(j5,j6)*za(j4,j1)**2)

      return
      end


      function FPFPsc7(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc7
      integer:: j1,j2,j3,j4,j5,j6

      FPFPsc7=
     & -za(j1,j3)**2*za(j5,j3)*zb(j6,j1)
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))

      return
      end


      function FPFPsc8(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc8
      integer:: j1,j2,j3,j4,j5,j6

      FPFPsc8=
     & -0.5_dp*za(j2,j3)*za(j5,j4)**2*zb(j4,j2)**2
     & /(za(j1,j4)*za(j5,j6)*za(j2,j4))

      return
      end


      function FPFPsc9(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc9
      integer:: j1,j2,j3,j4,j5,j6

      FPFPsc9=
     & -0.5_dp*za(j5,j2)**2*zb(j2,j4)**2*za(j4,j3)
     & /(za(j1,j2)*za(j5,j6)*za(j2,j4))

      return
      end


      function FPFPsc10(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc10
      integer:: j1,j2,j3,j4,j5,j6

      FPFPsc10=
     & -za(j1,j3)*za(j3,j4)*za(j5,j6)*zb(j6,j4)**2
     & /(za(j1,j2)*za(j2,j3)*za(j4,j1))

      return
      end


      function FPFPsc11(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFPsc11
      integer:: j1,j2,j3,j4,j5,j6

      FPFPsc11=
     & -0.5_dp*za(j1,j3)**3*za(j5,j6)*zb(j6,j1)**2
     & /(za(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j1))

      return
      end


