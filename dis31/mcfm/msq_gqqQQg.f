!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine msq_gqqQQg(i1,i2,i3,i4,i5,i6,i7,ea,eb,nll,MN,MI)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     March, 2001.                                                     *
c     Return matrix elements squared as a function of f1 and f2        *
c     summed over helicity, using the formulae of                      *
c     Nagy and Trocsnayi, PRD59 014020 (1999)                          *
c***********************************************************************
c disorder (2026-10-02): photon-exchange version of msq_ZqqQQg. The
c couplings are the electric charges ea of the quark line i1-i2
c (flavour f1) and eb of the line i3-i4 (f3) (makemb_photon), so MN, MI are numbers (the
c former arrays over f1, f3 have one entry); identical quarks: ea = eb.
c disorder (2026-10-08): ea, eb per helicity of their line and of the
c lepton line (dis31/ew31). nll: drop the interference of the boson on the
c two lines within a pairing class (the terms A, D, F, G: they are taken
c as the sum of those with the boson on line i1-i2 only and on i3-i4
c only), keeping the direct-exchange interference of identical quarks
c (B, C, E), as dis31/me31 for photon + Z.
      include 'constants.f'
      integer:: Qh,hq,hg,lh,f1,f3,i1,i2,i3,i4,i5,i6,i7,j
      logical:: nll
      real(dp):: A(2,2,2),B(2,2,2),C(2,2,2),D(2,2,2),E(2,2,2),
     & F(2,2,2),G(2,2,2)
      real(dp):: MI,MN,ea(2,2),eb(2,2),M0(2,2,2),Mx(2,2,2),My(2,2,2),
     & Mz(2,2,2),Mxx(2,2,2),Mxy(2,2,2),e0(2,2)
      real(dp):: x,y,z
      parameter(x=xn/cf,y=half/cf,z=0.25_dp*(xn**2-two)/xn/cf**2)
      complex(dp)::
     &               mb1_1234(1,1,2,2,2,2),mb2_1234(1,1,2,2,2,2),
     &               mb1_3412(1,1,2,2,2,2),mb2_3412(1,1,2,2,2,2),
     &               mb1_3214(1,1,2,2,2,2),mb2_3214(1,1,2,2,2,2),
     &               mb1_1432(1,1,2,2,2,2),mb2_1432(1,1,2,2,2,2)


c---set everything to zero

      f1=1
      f3=1
      MI=zip
      MN=zip
      e0=zip
      do j=1,2
      A(j,f1,f3)=zip
      B(j,f1,f3)=zip
      C(j,f1,f3)=zip
      D(j,f1,f3)=zip
      E(j,f1,f3)=zip
      F(j,f1,f3)=zip
      G(j,f1,f3)=zip
      enddo

c---mb1_1234 etc, has 6 indices each with possible values 1 or 2
c---corresponding to f1,f3,hq,Qh,hg,lh
      if (nll) then
        call makemb(ea,e0)
        call same()
        call makemb(e0,eb)
        call same()
        call makemb(ea,eb)
      else
        call makemb(ea,eb)
        call same()
      endif

      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
      if (hq == Qh) then
      B(2,f1,f3)=B(2,f1,f3)-two*real(
     &+mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_1432(f1,f3,hq,Qh,hg,lh))
     &+mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3214(f3,f1,Qh,hq,hg,lh))
     &+mb1_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb1_3214(f3,f1,Qh,hq,hg,lh))
     &+mb2_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh)))

      C(2,f1,f3)=C(2,f1,f3)-two*real(
     &+mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3214(f3,f1,Qh,hq,hg,lh))
     &+mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh))
     &+mb1_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb1_1432(f1,f3,hq,Qh,hg,lh))
     &+mb2_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb2_3214(f3,f1,Qh,hq,hg,lh)))

      E(2,f1,f3)=E(2,f1,f3)-two*real(
     &+(mb1_1234(f1,f3,hq,Qh,hg,lh)+mb1_3412(f3,f1,Qh,hq,hg,lh))
     &*conjg(mb2_3214(f3,f1,Qh,hq,hg,lh)+mb2_1432(f1,f3,hq,Qh,hg,lh))
     &+(mb2_1234(f1,f3,hq,Qh,hg,lh)+mb2_3412(f3,f1,Qh,hq,hg,lh))
     &*conjg(mb1_3214(f3,f1,Qh,hq,hg,lh)+mb1_1432(f1,f3,hq,Qh,hg,lh)))
      endif
      enddo
      enddo
      enddo
      enddo

      do j=1,2
      M0(j,f1,f3)=B(j,f1,f3)+C(j,f1,f3)+E(j,f1,f3)
      Mx(j,f1,f3)=-0.5_dp*(3._dp*C(j,f1,f3)+2._dp*E(j,f1,f3)+B(j,f1,f3))
      My(j,f1,f3)=A(j,f1,f3)+D(j,f1,f3)
      Mz(j,f1,f3)=F(j,f1,f3)+G(j,f1,f3)
      Mxx(j,f1,f3)=0.25_dp*(2._dp*C(j,f1,f3)+E(j,f1,f3))
      Mxy(j,f1,f3)=-0.5_dp*(F(j,f1,f3)+D(j,f1,f3))
      enddo
      MN=CF**3*xn*(M0(1,f1,f3)+x*Mx(1,f1,f3)+y*My(1,f1,f3)
     & +z*Mz(1,f1,f3)+x**2*Mxx(1,f1,f3)+x*y*Mxy(1,f1,f3))
      MI=
     & CF**3*xn*(M0(2,f1,f3)+x*Mx(2,f1,f3)+y*My(2,f1,f3)
     & +z*Mz(2,f1,f3)+x**2*Mxx(2,f1,f3)+x*y*Mxy(2,f1,f3))

      return

      contains

c---the four pairings with couplings ca (line i1-i2), cb (line i3-i4)
      subroutine makemb(ca,cb)
      real(dp):: ca(2,2),cb(2,2)
      call makemb_photon(i1,i2,i3,i4,i5,i6,i7,ca,cb,mb1_1234,mb2_1234)
      call makemb_photon(i3,i2,i1,i4,i5,i6,i7,cb,ca,mb1_3214,mb2_3214)
      call makemb_photon(i3,i4,i1,i2,i5,i6,i7,cb,ca,mb1_3412,mb2_3412)
      call makemb_photon(i1,i4,i3,i2,i5,i6,i7,ca,cb,mb1_1432,mb2_1432)
      end subroutine makemb

c---the terms within a pairing class: A, D, F, G
      subroutine same()
      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
      A(1,f1,f3)=A(1,f1,f3)
     & +abs(mb1_1234(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb2_1234(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb1_3412(f3,f1,Qh,hq,hg,lh))**2
     & +abs(mb2_3412(f3,f1,Qh,hq,hg,lh))**2
      D(1,f1,f3)=D(1,f1,f3)+two*real(
     &+mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_1234(f1,f3,hq,Qh,hg,lh))
     &+mb1_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh)))
      F(1,f1,f3)=F(1,f1,f3)+two*real(
     & +mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3412(f3,f1,Qh,hq,hg,lh))
     & +mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh)))
      G(1,f1,f3)=G(1,f1,f3)+two*real(
     & +mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh))
     & +mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3412(f3,f1,Qh,hq,hg,lh)))

      A(2,f1,f3)=A(2,f1,f3)
     & +abs(mb1_1234(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb2_1234(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb1_3214(f3,f1,Qh,hq,hg,lh))**2
     & +abs(mb2_3214(f3,f1,Qh,hq,hg,lh))**2
     & +abs(mb1_1432(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb2_1432(f1,f3,hq,Qh,hg,lh))**2
     & +abs(mb1_3412(f3,f1,Qh,hq,hg,lh))**2
     & +abs(mb2_3412(f3,f1,Qh,hq,hg,lh))**2
      D(2,f1,f3)=D(2,f1,f3)+two*real(
     &+mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_1234(f1,f3,hq,Qh,hg,lh))
     &+mb1_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb2_3214(f3,f1,Qh,hq,hg,lh))
     &+mb1_1432(f1,f3,hq,Qh,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh))
     &+mb1_3412(f3,f1,Qh,hq,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh)))
      F(2,f1,f3)=F(2,f1,f3)+two*real(
     & +mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3412(f3,f1,Qh,hq,hg,lh))
     & +mb1_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb1_1432(f1,f3,hq,Qh,hg,lh))
     & +mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh))
     & +mb2_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh)))
      G(2,f1,f3)=G(2,f1,f3)+two*real(
     & +mb1_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb2_3412(f3,f1,Qh,hq,hg,lh))
     & +mb1_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb2_1432(f1,f3,hq,Qh,hg,lh))
     & +mb2_1234(f1,f3,hq,Qh,hg,lh)*conjg(mb1_3412(f3,f1,Qh,hq,hg,lh))
     & +mb2_3214(f3,f1,Qh,hq,hg,lh)*conjg(mb1_1432(f1,f3,hq,Qh,hg,lh)))
      enddo
      enddo
      enddo
      enddo
      end subroutine same

      end
