!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine makemb_photon(i1,i2,i3,i4,i5,i6,i7,ea,eb,mb1,mb2)
      implicit none
      include 'types.f'

c     Author: R.K. Ellis, March 2001
c     A subroutine calculating Nagy and Trocsanyi, PRD59 014020 (1999)
c     Eq. B.56 with a factor of 2*i*e^2*g^3/s removed
c disorder (2026-10-02): makemb with makem inlined, photon exchange only:
c the couplings ea of the line i1-i2 and eb of the line i3-i4 are arguments
c (makem: Q(f)*q1 + Z terms from common blocks); flavour-type indices of
c size 1. disorder (2026-10-08): ea, eb per helicity of their line and
c of the lepton line, ea(hq,lh), eb(Qh,lh) (dis31/ew31)
      integer:: hq,Qh,hg,lh,i1,i2,i3,i4,i5,i6,i7
      real(dp):: ea(2,2),eb(2,2)
      complex(dp):: mb1(1,1,2,2,2,2),mb2(1,1,2,2,2,2),
     &     a1(2,2,2,2),a2(2,2,2,2),a3(2,2,2,2),a4(2,2,2,2),
     &     b1(2,2,2,2),b2(2,2,2,2),b3(2,2,2,2),b4(2,2,2,2)

      call nagyqqQQg(i1,i2,i3,i4,i5,i6,i7,a1,a2,a3,a4)
      call nagyqqQQg(i3,i4,i1,i2,i5,i6,i7,b1,b2,b3,b4)
      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
      mb1(1,1,hq,Qh,hg,lh)=ea(hq,lh)*a1(hq,Qh,hg,lh)
     &                    +eb(Qh,lh)*b3(Qh,hq,hg,lh)
      mb2(1,1,hq,Qh,hg,lh)=ea(hq,lh)*a2(hq,Qh,hg,lh)
     &                    +eb(Qh,lh)*b4(Qh,hq,hg,lh)
      enddo
      enddo
      enddo
      enddo
      return
      end
