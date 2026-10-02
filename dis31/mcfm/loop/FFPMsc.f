!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FFPMsc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq10.20
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp)::FFPMsc,FFPMsc_unsym
      integer:: j1,j2,j3,j4,j5,j6
      FFPMsc=FFPMsc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
     &      +FFPMsc_unsym(j2,j1,j4,j3,j6,j5,zb,za)
      return
      end

      function FFPMsc_unsym(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FFPMsc_unsym,Master2,Master2a,Master3a,zab2,zba2
      complex(dp):: zab,zba,FFPMscT,FFPMscTtilde,za12a
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,L1,Lsm1_2mh,I3m,Lnrat,Ls1,Lsm1
      complex(dp)::FFPMsc1,FFPMsc2,FFPMsc3,FFPMsc4,FFPMsc5
      complex(dp)::FFPMsc6,FFPMsc7,FFPMsc8,FFPMsc9,FFPMsc10,FFPMsc11
      real(dp):: t,delta12,delta56,deltat14,deltat56
      real(dp):: Delta3,Deltat3,s12,s23,s34,s56,s123,s234
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
      za12a(j1,j2,j3,j4,j5)=za(j1,j2)*zba2(j2,j3,j4,j5)
c     End statement functions

      deltat14=s(j1,j4)-s(j2,j3)-s(j5,j6)
      deltat56=s(j5,j6)-s(j1,j4)-s(j2,j3)
      delta12=s(j1,j2)-s(j5,j6)-s(j3,j4)
      delta56=s(j5,j6)-s(j1,j2)-s(j3,j4)
      Delta3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2*s(j1,j2)*s(j3,j4)-2*s(j3,j4)*s(j5,j6)-2*s(j5,j6)*s(j1,j2)
      Deltat3=s(j1,j4)**2+s(j2,j3)**2+s(j5,j6)**2
     & -2*s(j1,j4)*s(j2,j3)-2*s(j2,j3)*s(j5,j6)-2*s(j5,j6)*s(j1,j4)

      s123=t(j1,j2,j3)
      s234=t(j2,j3,j4)
      s12=s(j1,j2)
      s23=s(j2,j3)
      s34=s(j3,j4)
      s56=s(j5,j6)

c      FFPMsc1=
c     & +zb(j1,j3)/(zb(j5,j6)*za(j1,j3)*zab2(j3,j1,j2,j4))
c     * *(za(j1,j2)*zab2(j3,j1,j2,j6)**2/zab2(j3,j1,j2,j4))
c      FFPMsc2=
c     & +zb(j1,j3)/(zb(j5,j6)*za(j1,j3)*zab2(j3,j1,j2,j4))
c     & *(-za(j2,j3)*zab2(j1,j2,j3,j6)**2/zab2(j1,j2,j3,j4))

c      FFPMsc3=
c     & +0.5_dp*zb(j2,j3)*za(j2,j4)**2/(za(j5,j6)*zb(j2,j4)*t(j2,j3,j4))
c     & *(zb(j2,j3)*zab2(j5,j2,j3,j4)**2/zab2(j1,j2,j3,j4))

c      FFPMsc4=
c     & +0.5_dp*zb(j2,j3)*za(j2,j4)**2/(za(j5,j6)*zb(j2,j4)*t(j2,j3,j4))
c     & *(-zb(j3,j4)*zab2(j5,j3,j4,j2)**2/zab2(j1,j3,j4,j2))

c      FFPMsc5=
c     & +za(j2,j4)*(zb(j4,j6)*t(j1,j2,j3))**2
c     & /(zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)**2)

c      FFPMsc6=
c     & -zb(j2,j3)**2*za(j2,j4)*za(j1,j5)**2*t(j2,j3,j4)
c     & /(za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2)**2)

c      FFPMsc7=
c     & +0.5_dp*zb(j1,j3)*zb(j2,j3)*za(j1,j5)**2
c     & *t(j2,j3,j4)*zab2(j1,j2,j4,j3)
c     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2))
c     & *zb(j1,j3)/zb(j2,j3)

c      FFPMsc8=
c     & +0.5_dp*zb(j1,j3)*zb(j2,j3)*za(j1,j5)**2
c     & *t(j2,j3,j4)*zab2(j1,j2,j4,j3)
c     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2))
c     & *(-2._dp/zab2(j1,j3,j4,j2))

c       FFPMsc9=
c     & +0.5_dp*zb(j1,j3)*zb(j4,j6)
c     & /(zab2(j1,j2,j3,j4)*zab2(j3,j1,j4,j2)*zab2(j3,j1,j2,j4))
c     & *(za(j1,j3)*za(j4,j5)
c     & +(zab2(j3,j1,j4,j6)*t(j1,j2,j3)-zab2(j3,j1,j2,j6)*t(j1,j3,j4))
c     & /(zb(j1,j4)*zb(j5,j6)))

c      FFPMsc10=
c     & +(3*zb(j2,j3)*zab2(j2,j3,j4,j1)/(zab2(j3,j1,j2,j4)*Deltat3**2)
c     & *(za(j4,j5)*(zab(j2,j3,j6)*deltat23-zab(j2,j5,j6)*deltat56)
c     & -za(j1,j4)*zb(j1,j6)
c     & *(za(j2,j5)*deltat14 -2*zab(j2,j3,j6)*za(j6,j5)))

c     & +0.5_dp/(zab2(j3,j1,j2,j4)*Deltat3)
c     * *((-za(j2,j5)*za(j4,j5)/za(j5,j6)
c     & +zb(j1,j6)*zab2(j2,j1,j4,j6)/(zb(j1,j4)*zb(j5,j6)))
c     & *(zb(j1,j3)*deltat23-zb21b(j1,j5,j6,j2,j3))

c     & -2._dp*zb(j3,j6)*zab2(j2,j3,j4,j1)/(zb(j1,j4)*zb(j5,j6))
c     & *(zb(j1,j6)*deltat23-2*zba(j1,j4,j5)*zb(j5,j6))
c     & -4*zab(j4,j1,j3)*zb(j1,j6)*za(j2,j5)))

c      FFPMsc11=
c     & -0.5_dp*za(j3,j5)*zab2(j4,j1,j2,j3)*zab2(j5,j2,j4,j3)
c     & /(s(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4))

      FFPMsc_unsym=zb(j2,j4)*zab2(j3,j2,j4,j1)
     & /(zb(j1,j4)*zab2(j3,j1,j2,j4))
     & *Master2(j1,j4,j2,j3,j5,j6,za,zb)

     & -zb(j2,j3)**2*za(j2,j4)**2*zab2(j5,j2,j3,j4)**2
     & /(za(j5,j6)*zb(j2,j4)*s234*zab2(j1,j2,j3,j4))
     & *Ls1(-s23,-s234,-s34,-s234)
     & /s234**2

     & +(za(j2,j3)*zb(j4,j6)*s123)**2
     & /(za(j2,j3)*zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)**3)
     & *Lsm1_2mh(s34,s123,s12,s56)

     & -(zb(j2,j3)*za(j1,j5))**2*s234*zab2(j1,j2,j4,j3)
     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2)**3)
     & *Lsm1_2mh(s12,s234,s34,s56 )

     & +(za(j2,j3)*zab2(j1,j2,j3,j6))**2
     & /(za(j2,j3)*zb(j5,j6)*za(j1,j3)**2
     & *zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4))
     & *Lsm1(-s12,-s123,-s23,-s123)

     & +FFPMscT(j1,j2,j3,j4,j5,j6,za,zb)*I3m(s12,s34,s56)
     & +FFPMscTtilde(j1,j2,j3,j4,j5,j6,za,zb)*I3m(s(j1,j4),s23,s56)

     & +FFPMsc1(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s123,-s12)/s12
     & +FFPMsc2(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s123,-s23)/s23
     & +FFPMsc3(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s234,-s23)/s23**2
     & +FFPMsc4(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s234,-s34)/s34**2
     & +FFPMsc5(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s123,-s56)/s56
     & +FFPMsc6(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s34)/s34
     & +FFPMsc7(j1,j2,j3,j4,j5,j6,za,zb)*L1(-s56,-s234)/s234**2
     & +FFPMsc8(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s56,-s234)/s234
     & +FFPMsc9(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s23,-s56)
     & +FFPMsc10(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s23,-s56)
     & +FFPMsc11(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s(j1,j2),-s56)

     & -Master2a(j1,j2,j4,j3,j5,j6,za,zb)
     & -Master3a(j1,j2,j4,j3,j5,j6,za,zb)

     & +0.5_dp*zb(j4,j6)*zab2(j4,j1,j2,j3)
     & /(zb(j5,j6)*zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4)*Delta3)
     & *(zab2(j3,j1,j2,j6)*(za(j2,j4)*delta56
     & -2*zab(j2,j1,j3)*za(j3,j4))/za(j3,j4)
     & +za(j1,j2)*zb(j4,j6)
     & *(zb(j1,j3)*delta56-2*zba(j1,j2,j4)*zb(j4,j3))/zb(j3,j4))

     & +0.5_dp*za(j1,j5)*zab2(j4,j1,j2,j3)
     & /(za(j5,j6)*zab2(j1,j3,j4,j2)*zab2(j1,j2,j3,j4)*Delta3)
     & *((zab(j5,j2,j1)*delta12-zab(j5,j6,j1)*delta56)
     & *za(j1,j4)/za(j3,j4)
     & -zb(j2,j3)*(za(j2,j5)*delta56-2*za12a(j2,j1,j3,j4,j5)))

     & +0.5_dp*zab2(j2,j3,j4,j1)/(zab2(j3,j1,j2,j4)*Deltat3)
     & *(0.5_dp*deltat56*za(j2,j5)*zb(j1,j6)/(zb(j1,j4)*za(j2,j3))
     & +deltat14*za(j4,j5)*za(j2,j5)/(za(j2,j3)*za(j5,j6))
     & +za(j4,j5)*zb(j3,j6))

     & +0.5_dp/(zab2(j3,j1,j2,j4)*zab2(j1,j2,j3,j4))
     & *(-za(j1,j2)*zb(j1,j3)**2*zb(j4,j6)**2
     & /(zb(j1,j4)*zb(j3,j4)*zb(j5,j6))
     & -za(j2,j5)*za(j4,j5)*zab2(j2,j1,j3,j4)/(za(j2,j3)*za(j5,j6))
     & +za(j2,j4)*za(j4,j5)*zb(j4,j6)/za(j3,j4))

     & -0.5_dp*za(j2,j4)**2*zab2(j5,j3,j4,j2)*zab2(j5,j2,j3,j4)
     & /(za(j2,j3)*za(j3,j4)*za(j5,j6)*t(j2,j3,j4)*zb(j2,j4)
     & *zab2(j1,j2,j3,j4))

     & +0.5_dp*za(j1,j5)*za(j4,j5)*zab2(j4,j1,j2,j3)
     & /(za(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2)*zab2(j1,j2,j3,j4))

      return
      end

      function FFPMsc1(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMsc1,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc1=
     & +zb(j1,j3)/(zb(j5,j6)*za(j1,j3)*zab2(j3,j1,j2,j4))
     & *(za(j1,j2)*zab2(j3,j1,j2,j6)**2/zab2(j3,j1,j2,j4))

      return
      end


      function FFPMsc2(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMsc2,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc2=
     & +zb(j1,j3)/(zb(j5,j6)*za(j1,j3)*zab2(j3,j1,j2,j4))
     & *(-za(j2,j3)*zab2(j1,j2,j3,j6)**2/zab2(j1,j2,j3,j4))

      return
      end


      function FFPMsc3(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMsc3,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc3=
     & +0.5_dp*zb(j2,j3)*za(j2,j4)**2/(za(j5,j6)*zb(j2,j4)*t(j2,j3,j4))
     & *(zb(j2,j3)*zab2(j5,j2,j3,j4)**2/zab2(j1,j2,j3,j4))

      return
      end


      function FFPMsc4(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMsc4,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc4=
     & +0.5_dp*zb(j2,j3)*za(j2,j4)**2/(za(j5,j6)*zb(j2,j4)*t(j2,j3,j4))
     & *(-zb(j3,j4)*zab2(j5,j3,j4,j2)**2/zab2(j1,j3,j4,j2))

      return
      end


      function FFPMsc5(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMsc5,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc5=
     & +za(j2,j4)*(zb(j4,j6)*t(j1,j2,j3))**2
     & /(zb(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4)**2)

      return
      end


      function FFPMsc6(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMsc6,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc6=
     & -zb(j2,j3)**2*za(j2,j4)*za(j1,j5)**2*t(j2,j3,j4)
     & /(za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2)**2)

      return
      end


      function FFPMsc7(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMsc7,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc7=
     & +0.5_dp*zb(j1,j3)*zb(j2,j3)*za(j1,j5)**2
     & *t(j2,j3,j4)*zab2(j1,j2,j4,j3)
     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2))
     & *zb(j1,j3)/zb(j2,j3)

      return
      end


      function FFPMsc8(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMsc8,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc8=
     & +0.5_dp*zb(j1,j3)*zb(j2,j3)*za(j1,j5)**2
     & *t(j2,j3,j4)*zab2(j1,j2,j4,j3)
     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4)*zab2(j1,j3,j4,j2))
     & *(-2._dp/zab2(j1,j3,j4,j2))

      return
      end


      function FFPMsc9(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FFPMsc9,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc9=
     & +0.5_dp*zb(j1,j3)*zb(j4,j6)
     & /(zab2(j1,j2,j3,j4)*zab2(j3,j1,j4,j2)*zab2(j3,j1,j2,j4))
     & *(za(j1,j3)*za(j4,j5)
     & +(zab2(j3,j1,j4,j6)*t(j1,j2,j3)-zab2(j3,j1,j2,j6)*t(j1,j3,j4))
     & /(zb(j1,j4)*zb(j5,j6)))

      return
      end


      function FFPMsc10(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FFPMsc10,zab2,zba2
      complex(dp):: zab,zba,zb21b
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: deltat14,deltat23,deltat56
      real(dp):: Deltat3
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zba2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
      zba(j1,j2,j3)=zb(j1,j2)*za(j2,j3)
      zb21b(j1,j2,j3,j4,j5)=zba2(j1,j2,j3,j4)*zb(j4,j5)
c     End statement functions

      deltat14=s(j1,j4)-s(j2,j3)-s(j5,j6)
      deltat23=s(j2,j3)-s(j5,j6)-s(j1,j4)
      deltat56=s(j5,j6)-s(j1,j4)-s(j2,j3)
      Deltat3=s(j1,j4)**2+s(j2,j3)**2+s(j5,j6)**2
     & -2*s(j1,j4)*s(j2,j3)-2*s(j2,j3)*s(j5,j6)-2*s(j5,j6)*s(j1,j4)

      FFPMsc10=
     & +(3*zb(j2,j3)*zab2(j2,j3,j4,j1)/(zab2(j3,j1,j2,j4)*Deltat3**2)
     & *(za(j4,j5)*(zab(j2,j3,j6)*deltat23-zab(j2,j5,j6)*deltat56)
     & -za(j1,j4)*zb(j1,j6)
     & *(za(j2,j5)*deltat14 -2*zab(j2,j3,j6)*za(j6,j5)))

     & +0.5_dp/(zab2(j3,j1,j2,j4)*Deltat3)
     & *((-za(j2,j5)*za(j4,j5)/za(j5,j6)
     & +zb(j1,j6)*zab2(j2,j1,j4,j6)/(zb(j1,j4)*zb(j5,j6)))
     & *(zb(j1,j3)*deltat23-zb21b(j1,j5,j6,j2,j3))

     & -2._dp*zb(j3,j6)*zab2(j2,j3,j4,j1)/(zb(j1,j4)*zb(j5,j6))
     & *(zb(j1,j6)*deltat23-2*zba(j1,j4,j5)*zb(j5,j6))
     & -4*zab(j4,j1,j3)*zb(j1,j6)*za(j2,j5)))

      return
      end


      function FFPMsc11(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FFPMsc11,zab2
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FFPMsc11=
     & -0.5_dp*za(j3,j5)*zab2(j4,j1,j2,j3)*zab2(j5,j2,j4,j3)
     & /(s(j3,j4)*za(j5,j6)*zab2(j1,j3,j4,j2)*zab2(j3,j1,j2,j4))

      return
      end


