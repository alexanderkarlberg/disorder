!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function a61g(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a61g
c---hep-ph/9708239, Eqn 2.13
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      integer st
      complex(dp):: a6g,a6sg,a6fg,a6tg
      a61g=
     & +a6g(st,j1,j2,j3,j4,j5,j6,za,zb)
     & -a6g(st,j1,j4,j3,j2,j5,j6,za,zb)/xnsq
     & +(-real(nf,dp)*(a6sg(st,j1,j2,j3,j4,j5,j6,za,zb)
     &               +a6fg(st,j1,j2,j3,j4,j5,j6,za,zb))
     & +a6tg(st,j1,j2,j3,j4,j5,j6,za,zb))/xn
      return
      end

      function a63g(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a63g
c---hep-ph/9708239, Eqn 2.13
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'heldefs.f'
      integer:: j1,j2,j3,j4,j5,j6
      integer st,stamp(6)
      complex(dp):: a6g

c      write(*,*) 'Called a63g with arguments'
c      write(*,*) 'st',st
c      write(*,*) 'j1,j4,j2,j3,j5,j6',j1,j4,j2,j3,j5,j6
c      pause
c--- although a63 is labelled as (q1,qb4,g2,g3), the individual
c---  amplitudes are labelled by (q ,g ,g ,qb) [for 1 and 2]
c---                         and (q, g ,qb, g) [for 3 and 4]
c---                         and (q,qb , g, g) [for 5 and 6]

      if     (st == hqpqbmgpgp) then
        stamp(1)=hqpgpgpqbm
        stamp(2)=hqpgpgpqbm
        stamp(5)=hqpqbmgpgp
        stamp(6)=hqpqbmgpgp
      elseif (st == hqpqbmgpgm) then
        stamp(1)=hqpgpgmqbm
        stamp(2)=hqpgmgpqbm
        stamp(5)=hqpqbmgpgm
        stamp(6)=hqpqbmgmgp
      elseif (st == hqpqbmgmgp) then
        stamp(1)=hqpgmgpqbm
        stamp(2)=hqpgpgmqbm
        stamp(5)=hqpqbmgmgp
        stamp(6)=hqpqbmgpgm
      endif

      a63g=a6g(stamp(1),j1,j2,j3,j4,j5,j6,za,zb)
     &    +a6g(stamp(2),j1,j3,j2,j4,j5,j6,za,zb)
     &    +a6g(stamp(5),j1,j4,j2,j3,j5,j6,za,zb)
     &    +a6g(stamp(6),j1,j4,j3,j2,j5,j6,za,zb)

      if     (st == hqpqbmgpgp) then
        stamp(3)=hqpgpqbmgp
        stamp(4)=hqpgpqbmgp
        a63g=a63g
     &    +a6g(stamp(3),j1,j2,j4,j3,j5,j6,za,zb)
     &    +a6g(stamp(4),j1,j3,j4,j2,j5,j6,za,zb)
      elseif (st == hqpqbmgpgm) then
        stamp(3)=hqpgpqbmgm
        stamp(4)=hqpgpqbmgm
        a63g=a63g
     &    +a6g(stamp(3),j1,j2,j4,j3,j5,j6,za,zb)
     &    +a6g(stamp(4),j4,j3,j1,j2,j6,j5,zb,za)
      elseif (st == hqpqbmgmgp) then
        stamp(3)=hqpgpqbmgm
        stamp(4)=hqpgpqbmgm
        a63g=a63g
     &    +a6g(stamp(3),j4,j2,j1,j3,j6,j5,zb,za)
     &    +a6g(stamp(4),j1,j3,j4,j2,j5,j6,za,zb)
      else
        write(*,*) 'Unimplemented st in a63g'
        stop
      endif

      return
      end

      function a64vg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a64vg
c---hep-ph/9708239, Eqn 2.13
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      integer st
      complex(dp):: a6vsg,a6vfg
      a64vg=-a6vsg(st,j1,j4,j2,j3,j5,j6,za,zb)
     &      -a6vfg(st,j1,j4,j2,j3,j5,j6,za,zb)
      return
      end

      function a64axg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a64axg
c---hep-ph/9708239, Eqn 2.13
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      integer st
      complex(dp):: a6axg
      a64axg=a6axg(st,j1,j4,j2,j3,j5,j6,za,zb)
      return
      end

      function a65axg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a65axg
c---hep-ph/9708239, Eqn 2.13
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6axslg
      a65axg=a6axslg(st,j1,j4,j2,j3,j5,j6,za,zb)
      return
      end

      function a6sg(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6sg
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'heldefs.f'
      integer st
      integer:: j1,j2,j3,j4,j5,j6
c--- modified 4/18/00, according to BDK hep-ph/9708239, eq. (8.2)
      a6sg=czip
      if (st == hqpgpgpqbm)
     &  a6sg=one/three/za(j2,j3)**2/s(j5,j6)*(
     &        -(zb(j6,j1)*za(j1,j3)+zb(j6,j2)*za(j2,j3))*zb(j3,j1)
     &         *za(j4,j5)/(s(j1,j2)+s(j1,j3)+s(j2,j3))
     &        +(za(j5,j4)*zb(j4,j3)+za(j5,j2)*zb(j2,j3))*za(j3,j4)
     &         *zb(j1,j6)/(s(j2,j3)+s(j2,j4)+s(j3,j4)))
      return
      end

      function a6fg(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6fg
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer st
      integer:: j1,j2,j3,j4,j5,j6
c--- checked 4/18/00, this piece vanishes identically, eq. (8.3)
      a6fg=czip
      return
      end

      function a6tg(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'toploops.f'
      include 'heldefs.f'
      complex(dp):: a6tg,a6sg
      integer st
      integer:: j1,j2,j3,j4,j5,j6

      a6tg=czip
      if (st == hqpgpgpqbm) then
        if (toplight) then
          if     (toploops == tapprox) then
c--- modified 4/18/00, according to BDK hep-ph/9708239, eq. (8.3)
            a6tg=one/20._dp*s(j2,j3)/mt**2
          elseif (toploops == texact) then
c disorder (2026-10-03): exact top loops not ported (massless n_f = 5)
            stop 'a6tg: toploops = texact not available'
          endif
        endif
      endif
      a6tg=a6tg*a6sg(st,j1,j2,j3,j4,j5,j6,za,zb)

c Check that asymptotic limit is okay
c      do k=1,1000
c      s(j2,j3)=float(k)*1000._dp
c      write(6,*) sqrt(s(j2,j3)),one/20._dp*s(j2,j3)/mt**2/Ftexact(s(j2,j3),mt**2)
c      enddo
c      stop

      return
      end

      function a6vfg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6vfg
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer st
      integer:: j1,j2,j3,j4,j5,j6
      a6vfg=czip
      return
      end

      function a6vsg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6vsg
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer st
      integer:: j1,j2,j3,j4,j5,j6
      a6vsg=czip
      return
      end

      function a6axg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6axg
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer st
      integer:: j1,j2,j3,j4,j5,j6
      a6axg=czip
      return
      end

      function a6axslg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6axslg
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer st
      integer:: j1,j2,j3,j4,j5,j6
      a6axslg=czip
      return
      end
