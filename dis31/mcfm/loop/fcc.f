!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function Fcc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'heldefs.f'
      complex(dp):: Fcc
      integer:: j1,j2,j3,j4,j5,j6
      integer st
c Comments may be restored to test vs. old implementation
c      complex(dp):: Fcc_qpgmgpqm,Fcc_qpgpgmqm,Fcc_qpgpgpqm,
c     & Fcc_qpgpqmgm,Fcc_qpgpqmgp,Fcc_qpqmgmgp,Fcc_qpqmgpgm,Fcc_qpqmgpgp
      complex(dp)::FMPFcc,FPMFcc,FPPFcc,FPFMcc,FPFPcc,FFMPcc,FFPMcc,FFPPcc

      if     (st==hqpgmgpqbm) then
        Fcc=FMPFcc(j1,j2,j3,j4,j5,j6,za,zb)
c        write(6,*) Fcc/Fcc_qpgmgpqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st==hqpgpgmqbm) then
        Fcc=FPMFcc(j1,j2,j3,j4,j5,j6,za,zb)
c        write(6,*) Fcc/Fcc_qpgpgmqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st==hqpgpgpqbm) then
        Fcc=FPPFcc(j1,j2,j3,j4,j5,j6,za,zb)
c        write(6,*) Fcc/Fcc_qpgpgpqm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st==hqpgpqbmgm) then
        Fcc=FPFMcc(j1,j2,j3,j4,j5,j6,za,zb)
c        write(6,*) Fcc/Fcc_qpgpqmgm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st==hqpgpqbmgp) then
        Fcc=FPFPcc(j1,j2,j3,j4,j5,j6,za,zb)
c        write(6,*) Fcc/Fcc_qpgpqmgp(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st==hqpqbmgmgp) then
        Fcc=FFMPcc(j1,j2,j3,j4,j5,j6,za,zb)
c        write(6,*) Fcc/Fcc_qpqmgmgp(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st==hqpqbmgpgm) then
        Fcc=FFPMcc(j1,j2,j3,j4,j5,j6,za,zb)
c        write(6,*) Fcc/Fcc_qpqmgpgm(j1,j2,j3,j4,j5,j6,za,zb)
      elseif (st==hqpqbmgpgp) then
        Fcc=FFPPcc(j1,j2,j3,j4,j5,j6,za,zb)
c        write(6,*) Fcc/Fcc_qpqmgpgp(j1,j2,j3,j4,j5,j6,za,zb)
      else
        write(6,*) 'Error in Fcc: argument st is ',st
        stop
      endif

      return
      end

