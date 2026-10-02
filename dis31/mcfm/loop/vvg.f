!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!


c ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0._dp
c ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0._dp
c ADJUSTED TO REMOVE ALL 1/eps**2 PIECES IF EPINV2 IN EPINV2.F = 0._dp

      function vvg(st,j1,j2,j3,j4,j5,j6)
      implicit none
      include 'types.f'
      complex(dp):: vvg
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     July, 1998.                                                      *
c***********************************************************************
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'epinv2.f'
      include 'heldefs.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: Lnrat,xl12,xl34,xl23,xl56,Vcc,Vsc
      integer st


      xl12=Lnrat(musq,-s(j1,j2))
      xl34=Lnrat(musq,-s(j3,j4))
      xl23=Lnrat(musq,-s(j2,j3))
      xl56=Lnrat(musq,-s(j5,j6))

      if(st==hqpgmgpqbm) then
      Vcc=-four
     & -(epinv*epinv2+xl12*epinv+half*xl12**2)
     & -(epinv*epinv2+xl23*epinv+half*xl23**2)
     & -(epinv*epinv2+xl34*epinv+half*xl34**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st==hqpgpgmqbm) then
      Vcc=-four
     &   -(epinv*epinv2+xl12*epinv+half*xl12**2)
     &   -(epinv*epinv2+xl23*epinv+half*xl23**2)
     &   -(epinv*epinv2+xl34*epinv+half*xl34**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st==hqpgpgpqbm) then
      Vcc=-four
     &  -(epinv*epinv2+xl12*epinv+half*xl12**2)
     &  -(epinv*epinv2+xl23*epinv+half*xl23**2)
     &  -(epinv*epinv2+xl34*epinv+half*xl34**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st==hqpgpqbmgm) then
      Vcc=-four
     & -(+(epinv*epinv2+xl12*epinv+half*xl12**2)
     &   +(epinv*epinv2+xl23*epinv+half*xl23**2))-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)


      elseif(st==hqpgpqbmgp) then
      Vcc=-four
     & -(epinv*epinv2+xl12*epinv+half*xl12**2)
     & -(epinv*epinv2+xl23*epinv+half*xl23**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st==hqpqbmgmgp) then
      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st==hqpqbmgpgm) then
      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)

      elseif(st==hqpqbmgpgp) then
      Vcc=-four-(epinv*epinv2+xl12*epinv+half*xl12**2)-two*(epinv+xl56)
      Vsc=half*(one+epinv+xl56)
      else
      write(6,*) 'unimplemented st',st
      stop
      endif

      vvg=Vcc+Vsc
      end
