!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine a6routine(st,j1,j2,j3,j4,j5,j6,za,zb,a6sf,a6tp,a6uv)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     July, 1998.                                                      *
c***********************************************************************
c     a6sf is the sum of a6s and a6f
c     it is only the sum which is needed.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'masses.f'
      include 'epinv.f'
      include 'toploops.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp)::atree,virtsf,virtuv,virttp,Lnrat,a6sf,a6tp,a6uv,
     & tree
      character(len=2):: st
      logical:: msbar
      common/msbar/msbar

      if (st == 'sl') then
      write(6,*) 'error in a6routine',st
      stop
      endif
      tree=atree(st,j1,j2,j3,j4,j5,j6,za,zb)

      virttp=czip

      if (toplight) then
        if     (toploops == tapprox) then
          if (mt  /=  0._dp) then
            virttp=-2._dp/15._dp*s(j2,j3)/mt**2
          else
            stop 'mt=0 in a6routine'
          endif
        elseif (toploops == texact) then
c disorder (2026-10-03): exact top loops not ported (massless n_f = 5)
          stop 'a6routine: toploops = texact not available'
        elseif (toploops == tnone) then
          virttp=czip
        endif
      else
        virttp=czip
      endif

c Check that asymptotic limit is okay
c      do k=1,1000
c      s(j2,j3)=float(k)*1000._dp
c      write(6,*) sqrt(s(j2,j3)),-2._dp/15._dp*s(j2,j3)/mt**2/A6texact(s(j2,j3),mt**2)
c      enddo
c      stop

      virtsf=two/three*epinv
     & +two/three*Lnrat(musq,-s(j2,j3))+10._dp/9._dp
      virtuv=(epinv*(11._dp-two/xn*real(nf,dp))-one)/three
      if (msbar) virtuv=virtuv+two*CF/xn
c---virtuv is the infinite and finite renormalization
c---to get us to MS bar system
c--the term commented is associated with the number of legs
c  and in this program is taken care of in factorization procedure.
      a6sf=tree*virtsf
      a6tp=tree*virttp
      a6uv=tree*virtuv

      return
      end


