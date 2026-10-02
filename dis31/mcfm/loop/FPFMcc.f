!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function FPFMcc(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.10
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: FPFMcc,Master1,FPFMccT,FPFMccTtilde,zab2
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: L0,Lsm1_2mh,I3m,Lnrat,Lsm1
      complex(dp):: FPFMcc1,FPFMcc2,FPFMcc3,FPFMcc4,FPFMcc5
      real(dp):: t,s124,s234
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

c      FPFMcc1=
c     & -2*zb(j1,j2)*za(j1,j4)*zab2(j2,j1,j4,j6)*zab2(j4,j1,j2,j6)
c     & /(za(j1,j2)*zb(j5,j6)*t(j1,j2,j4)*zab2(j2,j1,j4,j3))

c      FPFMcc2=
c     &  +2*zb(j1,j6)*zab2(j1,j3,j4,j2)*zab2(j5,j3,j4,j2)
c     & /(zb(j3,j4)*t(j2,j3,j4)*zab2(j1,j2,j3,j4))

c      FPFMcc3=
c     & +2*zab(j5,j4,j2)*zab2(j5,j3,j4,j2)
c     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4))

c      FPFMcc4=
c     & +2*zab(j4,j3,j6)*zab2(j4,j1,j2,j6)
c     & /(za(j1,j2)*zb(j5,j6)*zab2(j2,j1,j4,j3))

c      FPFMcc5=
c     & -2*zab2(j5,j3,j4,j2)**2
c     & /(zb(j3,j4)*za(j5,j6)*t(j2,j3,j4)*zab2(j1,j2,j3,j4))

      s124=t(j1,j2,j4)
      s234=t(j2,j3,j4)

      FPFMcc=za(j1,j3)*zab2(j3,j1,j2,j6)**2
     & /(za(j1,j2)*za(j2,j3)*zb(j5,j6)
     & *zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4))
     & *Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))
c       fpfm_gf_b4
     & +(zab2(j5,j3,j4,j2)**2
     & /(zb(j3,j4)*za(j5,j6)*s234*zab2(j1,j2,j3,j4))
     & -(za(j3,j4)*zb(j1,j6))**2
     & /(za(j2,j3)*zb(j5,j6)*s234*zab2(j2,j3,j4,j1)))
     & *Lsm1_2mh(s(j1,j2),s234,s(j3,j4),s(j5,j6))
c       fpfm_gf_b3a
     & +(za(j1,j3)**3*(zb(j4,j6)*t(j1,j2,j3))**2
     & /(za(j1,j2)*za(j2,j3)*zb(j5,j6)
     & *zab2(j1,j2,j3,j4)**3*zab2(j3,j1,j2,j4))
     & -za(j1,j3)*zab2(j1,j2,j3,j6)**2*zab2(j3,j1,j2,j4)
     & /(za(j1,j2)*za(j2,j3)*zb(j5,j6)*zab2(j1,j2,j3,j4)**3))
     & *Lsm1_2mh( (s(j1,j4)),(t(j1,j2,j3)),(s(j2,j3)),(s(j5,j6)))
c       fpfm_gf_b1a
     & -((zab2(j2,j1,j4,j6)*zab2(j4,j1,j2,j3))**2
     & /(za(j1,j2)*zb(j5,j6)*s124*zab2(j2,j1,j4,j3)**3)
     & +(zb(j1,j2)*za(j3,j5))**2
     & /(zb(j1,j4)*za(j5,j6)*s124*zab2(j3,j1,j2,j4))
     & -(za(j2,j4)*zb(j3,j6))**2*s124
     & /(za(j1,j2)*zb(j5,j6)*zab2(j2,j1,j4,j3)**3))
     & *Lsm1_2mh(s(j2,j3),s124,s(j1,j4),s(j5,j6))
c       fpfm_gf_b5
     & +za(j1,j3)*zab2(j3,j1,j2,j6)**2
     & /(za(j1,j2)*za(j2,j3)*zb(j5,j6)
     & *zab2(j1,j2,j3,j4)*zab2(j3,j1,j2,j4))
     & *Lsm1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))
c       fpfm_gf_b1
     & +zab2(j5,j3,j4,j2)**2
     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4))
     & *Lsm1(-s(j2,j3),-s234,-s(j3,j4),-s234)
     & /s234
c       fpfm_gf_b5a
     & -(zb(j1,j2)*za(j3,j5))**2
     & /(za(j5,j6)*zb(j1,j4)*zab2(j3,j1,j2,j4))
     & *Lsm1(-s(j1,j4),-s124,-s(j1,j2),-s124)
     & /s124
c       fpfm_gf_3m
     & +FPFMccT(j1,j2,j3,j4,j5,j6,za,zb)
     & *I3m(s(j1,j2),s(j3,j4),s(j5,j6))
     & +FPFMccTtilde(j1,j2,j3,j4,j5,j6,za,zb)
     & *I3m(s(j1,j4),s(j2,j3),s(j5,j6))
c       master function terms
     & -Master1(j1,j4,j2,j3,j5,j6,za,zb)
     & -Master1(j3,j2,j4,j1,j6,j5,zb,za)
c       fpfm_gf_lzero_terms
     & +FPFMcc1(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s124,-s(j1,j4))/s(j1,j4)
     & +FPFMcc2(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j5,j6))/s(j5,j6)
     & +FPFMcc3(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s234,-s(j2,j3))/s(j2,j3)
     & +FPFMcc4(j1,j2,j3,j4,j5,j6,za,zb)*L0(-s124,-s(j5,j6))/s(j5,j6)
     & +FPFMcc5(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s(j2,j3),-s(j5,j6))


      return
      end

      function FPFMcc1(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.10
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMcc1,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FPFMcc1=
     & -2*zb(j1,j2)*za(j1,j4)*zab2(j2,j1,j4,j6)*zab2(j4,j1,j2,j6)
     & /(za(j1,j2)*zb(j5,j6)*t(j1,j2,j4)*zab2(j2,j1,j4,j3))

      return
      end


      function FPFMcc2(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.10
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMcc2,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FPFMcc2=
     &  +2*zb(j1,j6)*zab2(j1,j3,j4,j2)*zab2(j5,j3,j4,j2)
     & /(zb(j3,j4)*t(j2,j3,j4)*zab2(j1,j2,j3,j4))

      return
      end


      function FPFMcc3(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.10
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMcc3,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FPFMcc3=
     & +2*zab(j5,j4,j2)*zab2(j5,j3,j4,j2)
     & /(zb(j3,j4)*za(j5,j6)*zab2(j1,j2,j3,j4))

      return
      end


      function FPFMcc4(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.10
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMcc4,zab2,zab
      integer:: j1,j2,j3,j4,j5,j6
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      zab(j1,j2,j3)=za(j1,j2)*zb(j2,j3)
c     End statement functions

      FPFMcc4=
     & +2*zab(j4,j3,j6)*zab2(j4,j1,j2,j6)
     & /(za(j1,j2)*zb(j5,j6)*zab2(j2,j1,j4,j3))

      return
      end


      function FPFMcc5(j1,j2,j3,j4,j5,j6,za,zb)
c     Z. Bern, L. Dixon, D.A. Kosower, e-Print:hep-ph/9708239
c     Eq9.10
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: FPFMcc5,zab2
      integer:: j1,j2,j3,j4,j5,j6
      real(dp):: t
c     Statement functions
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
c     End statement functions

      FPFMcc5=
     & -2*zab2(j5,j3,j4,j2)**2
     & /(zb(j3,j4)*za(j5,j6)*t(j2,j3,j4)*zab2(j1,j2,j3,j4))

      return
      end


