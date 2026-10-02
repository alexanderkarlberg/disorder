!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFPPsc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FFPPsc,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lsm1_2me,Lnrat,Lsm1
      complex(dp):: FFPPsc1,FFPPsc2,FFPPsc3,FFPPsc4,FFPPsc5
      complex(dp):: FFPPsc6,FFPPsc7,FFPPsc8,FFPPsc9,FFPPsc10
      complex(dp):: FFPPsc11,FFPPsc12,FFPPsc13,FFPPsc14,FFPPsc15
      real(dp):: t,s12,s23,s56,s123,s234,s124
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

c      FFPPsc1=
c     & +((za(j4,j1)*za(j5,j2)-za(j2,j4)*za(j5,j1))*zab(j5,j4,j3)
c     &  /(za(j3,j4)*za(j4,j1)**2*za(j5,j6)))
c      FFPPsc2=
c     & -0.5_dp*zab(j5,j4,j3)**2*za(j3,j2)
c     & /(za(j3,j4)*za(j4,j1)*za(j5,j6))
c      FFPPsc3=
c     & -0.5_dp*za(j5,j2)**2
c     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))
c      FFPPsc4=
c     & +(za(j2,j4)*za(j5,j1)**2*zab2(j2,j3,j4,j1)
c     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)**2*za(j5,j6))
c     & -za(j5,j2)*zab(j2,j1,j6)/(za(j2,j3)*za(j3,j4)*za(j4,j1)))
c      FFPPsc5=
c     & -0.5_dp*za(j5,j6)*zab(j2,j1,j6)**2
c     & /(za(j2,j3)*za(j3,j4)*za(j4,j1))
c      FFPPsc6=
c     & +za(j2,j3)*zab2(j5,j1,j2,j3)*za(j3,j5)
c     & /(za(j3,j4)**2*za(j1,j3)*za(j5,j6))
c      FFPPsc7=
c     & - za(j1,j2)**2*zab2(j5,j2,j3,j1)*za(j1,j5)
c     & /(za(j2,j3)*za(j4,j1)**2*za(j1,j3)*za(j5,j6))
c      FFPPsc8=
c     & +(za(j2,j4)*zab(j5,j4,j6)/(za(j3,j4)**2*za(j4,j1))
c     & -(2*za(j4,j1)*za(j5,j2)+za(j5,j1)*za(j2,j4))*zab(j2,j4,j6)
c     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)**2)
c     & +za(j5,j2)*za(j5,j4)*zab2(j2,j1,j3,j4)
c     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6)))
c      FFPPsc9=
c     & -za(j5,j2)**2/(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))
c      FFPPsc10=
c     & -za(j5,j6)*zab(j2,j4,j6)**2
c     & /(za(j2,j3)*za(j3,j4)*za(j4,j1))
c      FFPPsc11=
c     & -za(j2,j4)*za(j5,j4)*zab2(j5,j1,j2,j4)
c     & /(za(j3,j4)**2*za(j4,j1)*za(j5,j6))

c      FFPPsc12=
c     & -za(j2,j4)*zab(j5,j3,j6)/(za(j3,j4)**2*za(j4,j1))

c      FFPPsc13=
c     & -za(j2,j3)*zb(j3,j6)**2*za(j5,j6)
c     & /(za(j3,j4)*za(j4,j1))

c      FFPPsc14=
c     & +za(j5,j1)*za(j5,j2)/(za(j4,j1)*za(j1,j3)*za(j5,j6))
c     & *(1._dp/za(j3,j4))
c      FFPPsc15=
c     & +za(j5,j1)*za(j5,j2)/(za(j4,j1)*za(j1,j3)*za(j5,j6))
c     & *(-za(j1,j2)/(za(j2,j3)*za(j4,j1)))

      s123=t(j1,j2,j3)
      s124=t(j1,j2,j4)
      s234=t(j2,j3,j4)
      s12=s(j1,j2)
      s23=s(j2,j3)
      s56=s(j5,j6)
      FFPPsc=-za(j2,j4)**2*za(j5,j1)**2
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)**3*za(j5,j6))
     & *Lsm1_2me(s123,s234,s23,s56)

     & -za(j2,j3)*za(j5,j1)**2
     & /(za(j3,j4)*za(j4,j1)*za(j1,j3)**2*za(j5,j6))
     & *Lsm1(-s12,-s123,-s23,-s123)

     & -za(j2,j3)*za(j5,j4)**2/(za(j3,j4)**3*za(j4,j1)*za(j5,j6))
     & *Lsm1_2me(s124,s123,s12,s56)

     & +FFPPsc1(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s23)/s23
     & +FFPPsc2(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s23,-s234)/s234**2
     & +FFPPsc3(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s234,-s56)
     & +FFPPsc4(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s56)/s56
     & +FFPPsc5(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s234,-s56)/s56**2
     & +FFPPsc6(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s123,-s12)/s12
     & +FFPPsc7(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s123,-s23)/s23
     & +FFPPsc8(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s123,-s56)/ s56
     & +FFPPsc9(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s123,-s56)
     & +FFPPsc10(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s123,-s56)/s56**2
     & +FFPPsc11(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s124,-s12)/s12
     & +FFPPsc12(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s124,-s56)/s56
     & +FFPPsc13(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s124,-s56)/s56**2
     & +FFPPsc14(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s12,-s56)
     & +FFPPsc15(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s23,-s56)

     & +0.5_dp*(
     & zb(j3,j4)*zab2(j4,j1,j3,j6)*zab2(j2,j3,j4,j6)
     & /(za(j3,j4)*za(j4,j1)*s234*t(j1,j3,j4)*zb(j5,j6))

     & +zb(j3,j4)*za(j5,j2)*zab(j5,j1,j3)
     & /(za(j4,j1)*s234*t(j1,j3,j4)*za(j5,j6))

     & +za(j5,j2)*zab(j5,j4,j3)
     & /(za(j3,j4)*za(j4,j1)*s234*za(j5,j6))

     & +za(j2,j5)**2/(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6)))

      return
      end


      function FFPPsc1(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc1,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPPsc1=
     & +((za(j4,j1)*za(j5,j2)-za(j2,j4)*za(j5,j1))*zab(j5,j4,j3)
     &  /(za(j3,j4)*za(j4,j1)**2*za(j5,j6)))

      return
      end


      function FFPPsc2(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc2,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPPsc2=
     & -0.5_dp*zab(j5,j4,j3)**2*za(j3,j2)
     & /(za(j3,j4)*za(j4,j1)*za(j5,j6))

      return
      end


      function FFPPsc3(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc3
      integer:: j1,j2,j3,j4,j5,j6

      FFPPsc3=
     & -0.5_dp*za(j5,j2)**2
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))

      return
      end


      function FFPPsc4(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc4,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPPsc4=
     & +(za(j2,j4)*za(j5,j1)**2*zab2(j2,j3,j4,j1)
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)**2*za(j5,j6))
     & -za(j5,j2)*zab(j2,j1,j6)/(za(j2,j3)*za(j3,j4)*za(j4,j1)))

      return
      end


      function FFPPsc5(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc5,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPPsc5=
     & -0.5_dp*za(j5,j6)*zab(j2,j1,j6)**2
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1))

      return
      end


      function FFPPsc6(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc6,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPPsc6=
     & +za(j2,j3)*zab2(j5,j1,j2,j3)*za(j3,j5)
     & /(za(j3,j4)**2*za(j1,j3)*za(j5,j6))

      return
      end


      function FFPPsc7(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc7,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPPsc7=
     & - za(j1,j2)**2*zab2(j5,j2,j3,j1)*za(j1,j5)
     & /(za(j2,j3)*za(j4,j1)**2*za(j1,j3)*za(j5,j6))

      return
      end


      function FFPPsc8(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc8,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPPsc8=
     & +(za(j2,j4)*zab(j5,j4,j6)/(za(j3,j4)**2*za(j4,j1))
     & -(2*za(j4,j1)*za(j5,j2)+za(j5,j1)*za(j2,j4))*zab(j2,j4,j6)
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)**2)
     & +za(j5,j2)*za(j5,j4)*zab2(j2,j1,j3,j4)
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6)))

      return
      end


      function FFPPsc9(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc9
      integer:: j1,j2,j3,j4,j5,j6

      FFPPsc9=
     & -za(j5,j2)**2/(za(j2,j3)*za(j3,j4)*za(j4,j1)*za(j5,j6))

      return
      end


      function FFPPsc10(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc10,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPPsc10=
     & -za(j5,j6)*zab(j2,j4,j6)**2
     & /(za(j2,j3)*za(j3,j4)*za(j4,j1))

      return
      end


      function FFPPsc11(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc11,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPPsc11=
     & -za(j2,j4)*za(j5,j4)*zab2(j5,j1,j2,j4)
     & /(za(j3,j4)**2*za(j4,j1)*za(j5,j6))

      return
      end


      function FFPPsc12(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc12,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FFPPsc12=
     & -za(j2,j4)*zab(j5,j3,j6)/(za(j3,j4)**2*za(j4,j1))

      return
      end


      function FFPPsc13(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc13
      integer:: j1,j2,j3,j4,j5,j6

      FFPPsc13=
     & -za(j2,j3)*zb(j3,j6)**2*za(j5,j6)
     & /(za(j3,j4)*za(j4,j1))

      return
      end


      function FFPPsc14(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc14
      integer:: j1,j2,j3,j4,j5,j6

      FFPPsc14=
     & +za(j5,j1)*za(j5,j2)/(za(j3,j4)*za(j4,j1)*za(j1,j3)*za(j5,j6))

      return
      end


      function FFPPsc15(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.6
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPPsc15
      integer:: j1,j2,j3,j4,j5,j6

      FFPPsc15=
     & +za(j5,j1)*za(j5,j2)/(za(j4,j1)*za(j1,j3)*za(j5,j6))
     & *(-za(j1,j2)/(za(j2,j3)*za(j4,j1)))

      return
      end

