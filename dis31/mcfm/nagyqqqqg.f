!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nagyqqQQg(i1,i2,i3,i4,i5,i6,i7,a1,a2,a3,a4)
      implicit none
      include 'types.f'
c----  %\cite{Nagy:1998bb}
c----  \bibitem{Nagy:1998bb}
c----  Z.~Nagy and Z.~Trocsanyi,
c----  %``Next-to-leading order calculation of four-jet
c----  observables in electron  positron annihilation,''
c----  Phys.\ Rev.\ D {\bf 59}, 014020 (1999)
c----  [Erratum-ibid.\ D {\bf 62}, 099902 (2000)]
c----  [arXiv:hep-ph/9806317].
c----  %%CITATION = HEP-PH 9806317;%%
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: j,k,j1,j2,j3,j4,j5,j6,j7,i1,i2,i3,i4,i5,i6,i7
      integer:: hq,Qh,lh,hg,h2,h3,h4
      real(dp):: s167,s267,s134,s234,s345
      complex(dp)::a1(2,2,2,2),a2(2,2,2,2),a3(2,2,2,2),a4(2,2,2,2),
     & t2a,xa(mxpart,mxpart),xb(mxpart,mxpart)
c-----statement function
      t2a(j1,j2,j3,j4)=xa(j1,j2)*xb(j2,j4)+xa(j1,j3)*xb(j3,j4)

      j1=i1
      j2=i2
      j5=i5


c----hq,Qh,hg,lh


      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
         if (hq == 2) then
            do j=1,mxpart
            do k=1,mxpart
            xa(j,k)=za(j,k)
            xb(j,k)=zb(j,k)
c            h1=2
            h2=Qh
            h3=hg
            h4=lh
            enddo
            enddo
         elseif (hq == 1) then
            do j=1,mxpart
            do k=1,mxpart
            xa(j,k)=zb(k,j)
            xb(j,k)=za(k,j)
            enddo
            enddo
c            h1=3-hq
            h2=3-Qh
            h3=3-hg
            h4=3-lh
         endif
         if (h2 == 2) then
            j3=i3
            j4=i4
         elseif (h2 == 1) then
            j3=i4
            j4=i3
         endif
         if (h4 == 2) then
            j6=i6
            j7=i7
         elseif (h4 == 1) then
            j6=i7
            j7=i6
         endif
      s167=s(j1,j6)+s(j1,j7)+s(j6,j7)
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)

      s134=s(j1,j3)+s(j1,j4)+s(j3,j4)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)

      if (h3 == 2) then
c +-+-+(A50)
       a1(hq,Qh,hg,lh)=
     & -xb(j1,j5)*t2a(j4,j1,j5,j3)*t2a(j4,j2,j6,j7)*xa(j6,j2)
     &  /(xa(j4,j5)*s(j1,j5)*s(j3,j4)*s267)
     & -xb(j1,j7)*t2a(j6,j1,j7,j5)*xa(j4,j2)**2*xb(j2,j3)
     &  /(xa(j4,j5)*s(j3,j4)*s234*s167)
     & -t2a(j4,j1,j5,j7)*t2a(j6,j2,j4,j3)*xa(j4,j2)
     &  /(xa(j1,j5)*xa(j5,j4)*s(j3,j4)*s234)
     & +xb(j5,j3)*t2a(j4,j3,j5,j1)*t2a(j4,j2,j6,j7)*xa(j6,j2)
     &  /(xa(j4,j5)*s(j3,j4)*s345*s267)
     & +xb(j1,j7)
     & *(t2a(j6,j1,j7,j3)*xa(j3,j4)+t2a(j6,j1,j7,j5)*xa(j5,j4))
     & *xb(j3,j5)*xa(j4,j2)/(xa(j4,j5)*s(j3,j4)*s345*s167)

c --(A53)
      a2(hq,Qh,hg,lh)=
     & -xb(j5,j3)*t2a(j4,j3,j5,j1)*t2a(j4,j2,j6,j7)*xa(j6,j2)
     &  /(xa(j4,j5)*s(j3,j4)*s345*s267)
     & -xb(j1,j7)*t2a(j6,j1,j7,j3)*xb(j2,j5)*xa(j4,j2)**2
     &  /(xa(j4,j5)*s(j2,j5)*s(j3,j4)*s167)
     & -xb(j1,j3)*t2a(j4,j1,j3,j5)*t2a(j4,j2,j6,j7)*xa(j6,j2)
     &  /(xa(j4,j5)*s(j3,j4)*s134*s267)
     & -xb(j1,j3)*t2a(j4,j1,j3,j7)*xa(j6,j2)*xa(j4,j2)
     &  /(xa(j4,j5)*xa(j5,j2)*s(j3,j4)*s134)
     & -xb(j1,j7)
     &  *(t2a(j6,j1,j7,j3)*xa(j3,j4)+t2a(j6,j1,j7,j5)*xa(j5,j4))
     & *xb(j3,j5)*xa(j4,j2)/(xa(j4,j5)*s(j3,j4)*s345*s167)

      if (h2==2) then
c ---(A56)
      a3(hq,Qh,hg,lh)=
     & -xb(j5,j3)*t2a(j4,j3,j5,j1)*t2a(j4,j2,j6,j7)*xa(j6,j2)
     & /(xa(j4,j5)*s(j3,j5)*s345*s267)
     & -xb(j1,j7)
     & *(t2a(j6,j1,j7,j3)*xa(j3,j4)+t2a(j6,j1,j7,j5)*xa(j5,j4))
     & *xb(j3,j5)*xa(j4,j2)/(xa(j4,j5)*s(j3,j5)*s345*s167)
      a4(hq,Qh,hg,lh)=czip
      else
      a3(hq,Qh,hg,lh)=czip
c ----(A59) with 3 and 4 exchanged, to compensate for above
      a4(hq,Qh,hg,lh)=+xb(j5,j3)*t2a(j4,j3,j5,j1)*t2a(j4,j2,j6,j7)
     & *xa(j6,j2)/(xa(j4,j5)*s(j3,j5)*s345*s267)
     &  +xb(j1,j7)
     & *(t2a(j6,j1,j7,j3)*xa(j3,j4)+t2a(j6,j1,j7,j5)*xa(j5,j4))
     & *xb(j3,j5)*xa(j4,j2)/(xa(j4,j5)*s(j3,j5)*s345*s167)
      endif

      elseif (h3 == 1)  then
c ---(A51)
      a1(hq,Qh,hg,lh)=
     & +xb(j1,j3)**2*xa(j5,j1)*t2a(j4,j2,j6,j7)*xa(j6,j2)
     &  /(xb(j3,j5)*s(j1,j5)*s(j3,j4)*s267)
     & +xb(j1,j7)*t2a(j6,j1,j7,j3)*t2a(j5,j2,j4,j3)*xa(j4,j2)
     &  /(xb(j3,j5)*s(j3,j4)*s234*s167)
     & -xb(j1,j3)*xb(j1,j7)*t2a(j6,j2,j4,j3)*xa(j4,j2)
     &  /(xb(j1,j5)*xb(j5,j3)*s(j3,j4)*s234)
     & +xb(j1,j3)*xa(j5,j4)
     &  *(+xb(j3,j4)*xa(j4,j2)*xb(j2,j7)+xb(j3,j5)*xa(j5,j2)*xb(j2,j7)
     &    +xb(j3,j4)*xa(j4,j6)*xb(j6,j7)+xb(j3,j5)*xa(j5,j6)*xb(j6,j7))
     & *xa(j6,j2)/(xb(j3,j5)*s(j3,j4)*s345*s267)
     & -xb(j1,j7)*t2a(j6,j1,j7,j3)*xa(j5,j4)*t2a(j2,j4,j5,j3)
     &  /(xb(j3,j5)*s(j3,j4)*s345*s167)

c ----(A54)
      a2(hq,Qh,hg,lh)=
     & -xb(j1,j3)*xa(j5,j4)*
     & (+xb(j3,j4)*xa(j4,j2)*xb(j2,j7)+xb(j3,j5)*xa(j5,j2)*xb(j2,j7)
     &  +xb(j3,j4)*xa(j4,j6)*xb(j6,j7)+xb(j3,j5)*xa(j5,j6)*xb(j6,j7))
     & *xa(j6,j2)/(xb(j3,j5)*s(j3,j4)*s345*s267)
     & +xb(j1,j3)**2*xa(j4,j1)*t2a(j5,j2,j6,j7)*xa(j6,j2)
     &  /(xb(j3,j5)*s(j3,j4)*s134*s267)
     & +xb(j1,j7)*t2a(j6,j1,j7,j3)*xa(j5,j4)*t2a(j2,j4,j5,j3)
     &  /(xb(j3,j5)*s(j3,j4)*s345*s167)
     & -xb(j1,j3)*t2a(j4,j1,j3,j7)*t2a(j6,j2,j5,j3)
     &  /(xb(j3,j5)*xb(j5,j2)*s(j3,j4)*s134)
     & +xb(j1,j7)*t2a(j6,j1,j7,j3)*xa(j5,j2)*t2a(j4,j2,j5,j3)
     &  /(xb(j3,j5)*s(j2,j5)*s(j3,j4)*s167)


      if (h2==1) then
c ----(A57) with 3 and 4 exchanged, to compensate for above
      a3(hq,Qh,hg,lh)=-xb(j1,j3)*xa(j5,j4)*
     & (+xb(j3,j4)*xa(j4,j2)*xb(j2,j7)+xb(j3,j5)*xa(j5,j2)*xb(j2,j7)
     &  +xb(j3,j4)*xa(j4,j6)*xb(j6,j7)+xb(j3,j5)*xa(j5,j6)*xb(j6,j7))
     & *xa(j6,j2)/(xb(j3,j5)*s(j4,j5)*s345*s267)
     & +xb(j1,j7)*t2a(j6,j1,j7,j3)*xa(j5,j4)*t2a(j2,j4,j5,j3)
     & /(xb(j3,j5)*s(j4,j5)*s345*s167)
      a4(hq,Qh,hg,lh)=czip
      else
      a3(hq,Qh,hg,lh)=czip
c ---(A58)
      a4(hq,Qh,hg,lh)=+xb(j1,j3)*xa(j5,j4)*
     & (+xb(j3,j4)*xa(j4,j2)*xb(j2,j7)+xb(j3,j4)*xa(j4,j6)*xb(j6,j7)
     &  +xb(j3,j5)*xa(j5,j2)*xb(j2,j7)+xb(j3,j5)*xa(j5,j6)*xb(j6,j7))
     & *xa(j6,j2)/(xb(j3,j5)*s(j4,j5)*s345*s267)
     & -xb(j1,j7)*t2a(j6,j1,j7,j3)*xa(j5,j4)*t2a(j2,j4,j5,j3)
     &  /(xb(j3,j5)*s(j4,j5)*s345*s167)
      endif

      endif

c--- put in photon propagator
      a1(hq,Qh,hg,lh)=a1(hq,Qh,hg,lh)/s(j6,j7)
      a2(hq,Qh,hg,lh)=a2(hq,Qh,hg,lh)/s(j6,j7)
      a3(hq,Qh,hg,lh)=a3(hq,Qh,hg,lh)/s(j6,j7)
      a4(hq,Qh,hg,lh)=a4(hq,Qh,hg,lh)/s(j6,j7)
      enddo
      enddo
      enddo
      enddo

      return
      end

