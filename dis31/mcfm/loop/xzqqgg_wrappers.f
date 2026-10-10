!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c disorder (2026-10-03): the helicity-configuration wrappers of MCFM 10.3's
c src/Zbb/xzqqgg.f (a6treeg1) and src/Zbb/xzqqgg_v.f (a61g1lc, a61g1slc,
c a61g1nf, a61gcol, a63g1; 2026-10-08: a64v, the boson on a closed quark
c loop, vector coupling), copied unchanged, for dis31/virt31.f90.

      function a6treeg1(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6treeg1

c----wrapper to a6treeg that also includes config st='q+g-g-qb-'
      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'heldefs.f'
      integer st
      complex(dp):: a6treeg

      if(st==hqpgmgmqbm) then
        a6treeg1=a6treeg(hqpgpgpqbm,j4,j3,j2,j1,j6,j5,zb,za)
      else
        a6treeg1=a6treeg(st,j1,j2,j3,j4,j5,j6,za,zb)
      endif

      return
      end

      function a61g1lc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a61g1lc

c----wrapper to a61g that also includes config st='qpgmgmqbm'
      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'heldefs.f'
      integer st
      complex(dp):: a61gcol

      if(st==hqpgmgmqbm) then
        a61g1lc=a61gcol(hqpgpgpqbm,j4,j3,j2,j1,j6,j5,zb,za,1)
      else
        a61g1lc=a61gcol(st,j1,j2,j3,j4,j5,j6,za,zb,1)
      endif

      return
      end

      function a61g1slc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a61g1slc

c----wrapper to a61g that also includes config st='qpgmgmqbm'
      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'heldefs.f'
      integer st
      complex(dp):: a61gcol

      if(st==hqpqbmgmgm) then
        a61g1slc=a61gcol(hqpqbmgpgp,j1,j2,j3,j4,j5,j6,za,zb,2)
      else
        a61g1slc=a61gcol(st,j4,j3,j2,j1,j6,j5,zb,za,2)
      endif

      return
      end

      function a61g1nf(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a61g1nf

c----wrapper to a61g that also includes config st='qpgmgmqbm'
      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'heldefs.f'
      integer st
      complex(dp):: a61gcol

      if(st==hqpgmgmqbm) then
        a61g1nf=a61gcol(hqpgpgpqbm,j4,j3,j2,j1,j6,j5,zb,za,3)
      else
        a61g1nf=a61gcol(st,j1,j2,j3,j4,j5,j6,za,zb,3)
      endif

      return
      end

      function a61gcol(st,j1,j2,j3,j4,j5,j6,za,zb,ncol)
      implicit none
      include 'types.f'
      complex(dp):: a61gcol

c---hep-ph/9708239, Eqn 2.13
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6,ncol
      integer st
      complex(dp):: a6g,a6sg,a6fg,a6tg

c--- use ncol=1 for leading colour piece, ncol=2 for subleading
      if     (ncol == 1) then
c--- comes with natural colour factor (1)
      a61gcol=
     & +a6g(st,j1,j2,j3,j4,j5,j6,za,zb)
      elseif (ncol == 2) then
c--- comes with natural colour factor (-1/xnsq)
      a61gcol=
     & +a6g(st,j1,j4,j3,j2,j5,j6,za,zb)
      elseif (ncol == 3) then
c--- comes with natural colour factor (1/xn)
      a61gcol=
     & -real(nf,dp)*(a6sg(st,j1,j2,j3,j4,j5,j6,za,zb)
     &             +a6fg(st,j1,j2,j3,j4,j5,j6,za,zb))
     & +a6tg(st,j1,j2,j3,j4,j5,j6,za,zb)
      endif
      return
      end

      function a63g1(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a63g1

c----wrapper to a63g that also includes config st='qpqbmgmgm'
      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'heldefs.f'
      integer st
      complex(dp):: a63g

      if(st==hqpqbmgmgm) then
        a63g1=a63g(hqpqbmgpgp,j4,j1,j2,j3,j6,j5,zb,za)
      else
        a63g1=a63g(st,j1,j4,j2,j3,j5,j6,za,zb)
      endif

      return
      end

      function a64v(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a64v

c----definition (2.13) of BDK, writes in terms of fvs and fvf
      integer:: j1,j2,j3,j4,j5,j6
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'heldefs.f'
      integer st
      complex(dp):: fvs,fvf

      if     (st==hqpqbmgmgm) then
        a64v=-fvs(hqpqbmgpgp,j4,j1,j3,j2,j6,j5,zb,za)
     &       -fvf(hqpqbmgpgp,j4,j1,j3,j2,j6,j5,zb,za)
      elseif (st==hqpqbmgmgp) then
        a64v=-fvs(hqpqbmgpgm,j1,j4,j3,j2,j5,j6,za,zb)
     &       -fvf(hqpqbmgpgm,j1,j4,j3,j2,j5,j6,za,zb)
      else
        a64v=-fvs(st,j1,j4,j2,j3,j5,j6,za,zb)
     &       -fvf(st,j1,j4,j2,j3,j5,j6,za,zb)
      endif

      return
      end
