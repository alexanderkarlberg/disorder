!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function Fsc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'heldefs.f'
      complex(dp):: Fsc
      integer:: j1,j2,j3,j4,j5,j6
      integer st
c Comments may be restored to test vs. old implementation
c      complex(dp):: Fsc1,Fsc2,Fsc3,Fsc4,Fsc5,Fsc6,Fsc7,Fsc8
      complex(dp):: FMPFsc,FPMFsc,FPPFsc,FPFMsc,FPFPsc,FFMPsc,FFPMsc,FFPPsc

      if(st==hqpgmgpqbm) then
c FMPFsc Eq.8.20
      Fsc=FMPFsc(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) Fsc/Fsc1(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st==hqpgpgmqbm) then
c FPMFsc Eq 8.14
      Fsc=FPMFsc(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) Fsc/Fsc2(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st==hqpgpgpqbm) then
c FPPFsc Eq.8.8
      Fsc=FPPFsc(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) Fsc/Fsc3(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st==hqpgpqbmgm) then
c FPFMsc Eq.9.14
      Fsc=FPFMsc(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) Fsc/Fsc4(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st==hqpgpqbmgp) then
c FPFPsc Eq 9.6
      Fsc=FPFPsc(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) Fsc/Fsc5(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st==hqpqbmgmgp) then
c FFMPsc Eq.10.12
      Fsc=FFMPsc(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) Fsc/Fsc6(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st==hqpqbmgpgm) then
c FFPMsc Eq 10.20
      Fsc=FFPMsc(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) Fsc/Fsc7(j1,j2,j3,j4,j5,j6,za,zb)
      elseif(st==hqpqbmgpgp) then
c FFPPsc Eq.10.6
      Fsc=FFPPsc(j1,j2,j3,j4,j5,j6,za,zb)
c      write(6,*) Fsc/Fsc8(j1,j2,j3,j4,j5,j6,za,zb)
      endif

      return
      end
