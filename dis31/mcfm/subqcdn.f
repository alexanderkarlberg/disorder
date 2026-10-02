!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine subqcdn(p1,p2,p3,p4,p5,p6,nDp5,za,zb,zab,zba,
     & qcdab,qcdba)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'nwz.f'
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
c     first argument is the gluon polarization
c     2nd argument is the fermion line
c     3nd argument is the lepton line

c     1 is left handed
c     2 is right handed

      complex(dp):: qcdab(2,2,2),qcdba(2,2,2),ab36,ab36x,ab64x,
     & ab64,ab35x,ab35,ab54,ab54x,ab31,ab24,aba23
      real(dp):: s34,t15,t25,t16,t26,t56,t156,t256,nDp5
      integer:: p1,p2,p3,p4,p5,p6,plep,pglu,pfer


      do pglu=1,2
      do pfer=1,2
      do plep=1,2
      qcdab(pglu,pfer,plep)=czip
      qcdba(pglu,pfer,plep)=czip
      enddo
      enddo
      enddo



      s34=s(p3,p4)
      t15=s(p1,p5)
      t16=s(p1,p6)
      t25=s(p2,p5)
      t26=s(p2,p6)
      t56=s(p5,p6)
      t156=t15+t16+t56
      t256=t25+t26+t56

      ab24=za(p1,p2)*zb(p1,p4)+za(p2,p5)*zb(p4,p5)
      ab31=za(p2,p3)*zb(p1,p2)-za(p3,p5)*zb(p1,p5)
      ab35=za(p1,p3)*zb(p1,p5)-za(p3,p4)*zb(p4,p5)
      ab36=za(p1,p3)*zb(p1,p6)-za(p3,p4)*zb(p4,p6)
      ab54=za(p2,p5)*zb(p2,p4)+za(p3,p5)*zb(p3,p4)
      ab64=za(p2,p6)*zb(p2,p4)+za(p3,p6)*zb(p3,p4)
      aba23=zab(p2,p1)*za(p1,p3)-zab(p2,p4)*za(p3,p4)


      qcdab(1,1,1) =
     & +za(p2,p5)/zb(p1,p5)*za(p3,p4)*zb(p1,p4)*zb(p1,p4)
     & *2._dp/t256/t56*nDp5
     & -(zab(p1,p1)*zb(p1,p4)-zab(p6,p1)*zb(p4,p6))
     & *ab31*za(p2,p5)/zb(p1,p5)/t25/t16
     & +ab36*zab(p5,p1)/zb(p1,p5)*za(p2,p6)*zb(p1,p4)/t256/t56
     & -(za(p5,p6)*zb(p1,p6)/t56+za(p2,p5)*zb(p1,p2)/t25)
     &  *aba23*zb(p1,p4)/zb(p1,p5)/t256
     & +(zab(p5,p4)*za(p3,p4)-zab(p5,p1)*za(p1,p3))
     & *za(p2,p5)*zb(p1,p4)/t256/t25
     & +(zab(p3,p1)*zb(p3,p4)+zab(p2,p1)*zb(p2,p4))
     & *za(p2,p3)*za(p5,p6)*zb(p1,p6)/zb(p1,p5)/t156/t56
     & -ab64*zab(p5,p1)*za(p2,p3)*zb(p1,p6)/zb(p1,p5)/t156/t56
     & -ab54*zab(p6,p1)*za(p2,p3)*zb(p1,p6)/zb(p1,p5)/t156/t16
      qcdab(1,1,1) =2._dp*qcdab(1,1,1)/s34

c      aRLAB=[TATB]*(2._dp*sw^2*gzle*Lu)*qcdab(2,1,1)
      qcdab(2,1,1)=
     & (-2._dp*za(p2,p3)*za(p2,p3)*zb(p1,p5)*zb(p3,p4)*nDp5
     & +zba(p5,p2)*za(p2,p3)*zb(p1,p6)*ab64
     & +za(p2,p3)*za(p2,p3)*zb(p3,p4)*t56/t16
     & *(zab(p1,p1)*zb(p1,p5)-zab(p6,p1)*zb(p5,p6))
     & -za(p2,p3)*za(p2,p6)*zb(p5,p6)
     & *(zab(p2,p1)*zb(p2,p4)+zab(p3,p1)*zb(p3,p4)))/t156
     & -(ab36*zba(p5,p2)*za(p2,p6)*zb(p1,p4)
     & -za(p2,p6)*zb(p1,p4)*zb(p5,p6)*aba23)/t256
      qcdab(2,1,1)=2._dp*qcdab(2,1,1)/za(p2,p5)/s34/t56

c      aLLBA=[TBTA]*qcdba(1,1,1)*(2._dp*sw^2*gzle*Lu)
      qcdba(1,1,1)=
     & (-2._dp*za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p1,p4)*nDp5
     & +(zab(p2,p2)*za(p2,p5)-zab(p2,p6)*za(p5,p6))
     & *za(p3,p4)*zb(p1,p4)*zb(p1,p4)*t56/t26
     & -ab36*zab(p5,p1)*za(p2,p6)*zb(p1,p4)
     & +aba23*za(p5,p6)*zb(p1,p4)*zb(p1,p6))/t256
     & +(zab(p5,p1)*za(p2,p3)*zb(p1,p6)*ab64
     & -(zab(p2,p1)*zb(p2,p4)+zab(p3,p1)*zb(p3,p4))
     & *za(p2,p3)*za(p5,p6)*zb(p1,p6))/t156
      qcdba(1,1,1)=2._dp*qcdba(1,1,1)/zb(p1,p5)/t56/s34

c      aRLBA=[TBTA]*(2._dp*sw^2*gzle*Lu)*qcdba(2,1,1)
      qcdba(2,1,1)=
     & za(p1,p2)*za(p2,p3)*t56/za(p1,p5)/za(p2,p5)
     & *(zab(p2,p1)*zb(p2,p4)/t156+zab(p2,p2)*zb(p1,p4)/t26)
     & +(zab(p2,p2)*za(p2,p3)*zb(p4,p5)-zab(p2,p6)*za(p3,p6)
     & *ab24/za(p2,p5))*t56/za(p1,p5)/t26
     & +((zab(p3,p1)*za(p1,p2)*zb(p3,p4)/za(p2,p5)
     & -(zab(p2,p5)*zb(p2,p4)+zab(p3,p5)*zb(p3,p4)))
     & *(za(p2,p3)*t56/za(p1,p5))
     & +(zab(p2,p1)*zb(p2,p4)+zab(p3,p1)*zb(p3,p4))
     & *za(p2,p3)*za(p2,p6)*zb(p5,p6)/za(p2,p5)
     & +(-zba(p5,p2)*za(p2,p3)*zb(p1,p6)*ab64
     & +2._dp*nDp5*za(p2,p3)*za(p2,p3)*zb(p1,p5)*zb(p3,p4))
     & /za(p2,p5))/t156
     & +(-zab(p2,p6)*zb(p1,p4)*ab35
     & *t56/za(p2,p5)/zb(p2,p6)
     & -aba23*za(p2,p6)*zb(p1,p4)*zb(p5,p6)/za(p2,p5)
     & +zba(p5,p2)*za(p2,p6)*zb(p1,p4)
     & *ab36/za(p2,p5))/t256
      qcdba(2,1,1)=2._dp*qcdba(2,1,1)/t56/s34

      if ((nwz == +1) .or. (nwz == -1)) return

      ab36x=za(p3,p4)*zb(p4,p6)-za(p2,p3)*zb(p2,p6)
      ab64x=za(p1,p6)*zb(p1,p4)+za(p3,p6)*zb(p3,p4)
      ab35x=za(p3,p4)*zb(p4,p5)-za(p2,p3)*zb(p2,p5)
      ab54x=za(p1,p5)*zb(p1,p4)+za(p3,p5)*zb(p3,p4)



c      aLRAB = qcdab(1,2,1)*sw/cw*(-2._dp*Qu*sw^2)*gzle
      qcdab(1,2,1) =
     & +ab64x*zab(p5,p2)*za(p1,p3)*zb(p2,p6)/t256/t56
     &  +ab36x*zab(p5,p2)*za(p1,p6)*zb(p2,p4)/t156/t56
     &  +(zba(p1,p1)*za(p1,p5)-zba(p6,p1)*za(p5,p6))
     &  *za(p3,p4)*zb(p2,p4)*zb(p2,p4)/t156/t16
     &  +(zba(p4,p1)*za(p3,p4)-zba(p2,p1)*za(p2,p3))
     &  *zb(p2,p4)*zb(p2,p6)/zb(p5,p6)/t156
     &  +(zba(p2,p1)*zb(p1,p4)+zba(p2,p3)*zb(p3,p4))
     &  *za(p1,p3)*zb(p2,p6)/zb(p5,p6)/t256
     &  -za(p1,p5)*za(p3,p4)*zb(p2,p4)*zb(p2,p4)
     &  *2._dp*nDp5/t156/t56
      qcdab(1,2,1) = 2._dp/s34*qcdab(1,2,1)/zb(p2,p5)

c      aRRAB = qcdab(2,2,1)*sw/cw*(-2._dp*Qu*sw^2)*gzle
      qcdab(2,2,1) =
     & -(zba(p1,p1)*za(p1,p3)-zba(p6,p1)*za(p3,p6))
     & *za(p1,p2)*zb(p2,p4)/za(p1,p5)/za(p2,p5)/t16
     & +(zba(p1,p1)*za(p1,p3)-zba(p6,p1)*za(p3,p6))
     & *zb(p4,p5)/za(p2,p5)/t16
     & -(zba(p2,p3)*zb(p3,p4)+zba(p2,p1)*zb(p1,p4))
     & *za(p1,p2)*za(p1,p3)/za(p1,p5)/za(p2,p5)/t256
     & -(zba(p2,p1)*zb(p1,p4)+zba(p2,p3)*zb(p3,p4))
     & *za(p1,p3)*za(p1,p6)/za(p1,p5)/za(p5,p6)/t256
     & +(zba(p2,p1)*za(p2,p3)- zba(p4,p1)*za(p3,p4))
     & *za(p1,p6)*zb(p2,p4)/za(p1,p5)/za(p5,p6)/t156
     & -ab64x
     & *zba(p5,p1)*za(p1,p3)*zb(p2,p6)/za(p1,p5)/t256/t56
     & -(zba(p5,p1)*zb(p1,p4)+zba(p5,p3)*zb(p3,p4))
     & *za(p1,p3)/za(p2,p5)/t256
     & -ab36x*zba(p5,p1)*za(p1,p6)*zb(p2,p4)/za(p1,p5)/t156/t56
     & +ab35x*zba(p6,p1)*zb(p2,p4)/za(p1,p5)/zb(p1,p6)/t156
     & +za(p1,p3)*za(p1,p3)*zb(p2,p5)*zb(p3,p4)/za(p1,p5)
     & *(2._dp/t256/t56*nDp5)
      qcdab(2,2,1) = 2._dp*qcdab(2,2,1)/s34


c      aLRBA = qcdba(1,2,1)*sw/cw*(-2._dp*Qu*sw^2)*gzle
      qcdba(1,2,1) =
     & -ab64x*zab(p5,p2)*za(p1,p3)*zb(p2,p6)/zb(p2,p5)/t256/t56
     & -ab36x*zab(p5,p2)*za(p1,p6)*zb(p2,p4)/zb(p2,p5)/t156/t56
     & -(zba(p2,p1)*zb(p1,p4)+zba(p2,p3)*zb(p3,p4))
     & *za(p1,p3)*zb(p2,p6)/zb(p2,p5)/zb(p5,p6)/t256
     & +ab54x*zba(p2,p6)*za(p1,p3)/za(p2,p6)/zb(p2,p5)/t256
     & +(zba(p2,p5)*za(p2,p3)-zba(p4,p5)*za(p3,p4))
     & *zb(p2,p4)/zb(p1,p5)/t156
     & +(zba(p4,p1)*za(p3,p4)-zba(p2,p1)*za(p2,p3))
     & *zb(p1,p2)*zb(p2,p4)/zb(p1,p5)/zb(p2,p5)/t156
     & +(zba(p2,p1)*za(p2,p3)-zba(p4,p1)*za(p3,p4))
     & *zb(p2,p4)*zb(p2,p6)/zb(p2,p5)/zb(p5,p6)/t156
     & +(zba(p2,p6)*zb(p4,p6)- zba(p2,p2)*zb(p2,p4))
     & *za(p1,p3)*zb(p1,p2)/zb(p1,p5)/zb(p2,p5)/t26
     & +(zba(p2,p6)*zb(p4,p6)-zba(p2,p2)*zb(p2,p4))
     & *za(p3,p5)/zb(p1,p5)/t26
     & +za(p1,p5)*za(p3,p4)*zb(p2,p4)*zb(p2,p4)/zb(p2,p5)
     & *(2._dp/t156/t56*nDp5)
      qcdba(1,2,1) = 2._dp/s34*qcdba(1,2,1)

c      aRRBA = qcdba(2,2,1)*sw/cw*(-2._dp*Qu*sw^2)*gzle
      qcdba(2,2,1) =
     & +(zba(p2,p1)*zb(p1,p4)+zba(p2,p3)*zb(p3,p4))
     &  *za(p1,p3)*za(p1,p6)/za(p5,p6)/t256
     &  +(zba(p2,p2)*zb(p2,p5)-zba(p2,p6)*zb(p5,p6))
     &  *za(p1,p3)*za(p1,p3)*zb(p3,p4)/t256/t26
     &  +(zba(p4,p1)*za(p3,p4)-zba(p2,p1)*za(p2,p3))
     &  *za(p1,p6)*zb(p2,p4)/za(p5,p6)/t156
     &  +ab64x*zba(p5,p1)*za(p1,p3)*zb(p2,p6)/t256/t56
     &  +ab36x*zba(p5,p1)*za(p1,p6)*zb(p2,p4)/t156/t56
     &  -za(p1,p3)*za(p1,p3)*zb(p2,p5)*zb(p3,p4)
     &  *2._dp*nDp5/t256/t56
      qcdba(2,2,1) = 2._dp/s34*qcdba(2,2,1)/za(p1,p5)


      qcdab(2,2,2)=conjg(qcdab(1,1,1))
      qcdab(1,1,2)=conjg(qcdab(2,2,1))
      qcdab(1,2,2)=conjg(qcdab(2,1,1))
      qcdab(2,1,2)=conjg(qcdab(1,2,1))

      qcdba(2,2,2)=conjg(qcdba(1,1,1))
      qcdba(1,1,2)=conjg(qcdba(2,2,1))
      qcdba(1,2,2)=conjg(qcdba(2,1,1))
      qcdba(2,1,2)=conjg(qcdba(1,2,1))

      return
      end
