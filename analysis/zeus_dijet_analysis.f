!  ZEUS-like dijets (validation of the tau_2-sliced NNLO DIS 2+1 against
!  NNLOJET, nnlo21/validation/nnlojet_epLJJ_zeus2j.run): no FastJet.
!  Selection (as NNLOJET's ecuts_dis): inclusive kt jets (R = 1,
!  E-scheme) in the Breit frame; jets with E_T = sqrt(p_T^2 + m^2) > 8
!  GeV; then jets outside -1 < eta_lab < 2.5 dropped (proton direction
!  = +eta); then >= 2 jets and m12 > 20 GeV for the two leading jets in
!  Breit p_T. The Q^2 and y cuts (125 < Q^2 < 20000 GeV^2, 0.2 < y <
!  0.6) are disorder's phase-space limits (-Q2min -Q2max -ymin -ymax).
!  Histograms: sigma, Q2, ptavg_12, m12 (sigma per bin in the output).

      subroutine define_histograms
      implicit none
      include 'pwhg_bookhist-multi.h'
      double precision, parameter :: q2e(7) = (/125d0, 250d0, 500d0,
     $     1000d0, 2000d0, 5000d0, 20000d0/)
      double precision, parameter :: pte(5) = (/8d0, 15d0, 22d0, 30d0,
     $     60d0/)
      double precision, parameter :: mje(5) = (/20d0, 30d0, 45d0, 65d0,
     $     120d0/)
      call bookupeqbins('sigma',1d0,0d0,1d0)
      call bookup('Q2',6,q2e)
      call bookup('ptavg12',4,pte)
      call bookup('m12',4,mje)
      end

      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_parameters
      use mod_analysis
      implicit none
      integer n
      double precision dsig(maxscales), x, y, Q2
      double precision pb(0:3,4), jets(0:3,4), pj(0:3), pt2(4)
      double precision boostsgn, zsgn, et2, zel, pzl, pabs
      double precision etal, m12, ptavg, w(0:3)
      integer np, nj, i, m, i1, i2
      np = n - 3
      if (np .lt. 2) return
!     Breit frame partons, oriented so that the incoming parton is +z
      zsgn = sign(1d0, pbreit(3,2))
      do i = 1, np
         pb(:,i) = pbreit(:,3+i)
         pb(3,i) = zsgn*pb(3,i)
      enddo
      call zkt(pb, np, jets, nj)
!     lab frame: the same jets boosted; use the lab momenta of the
!     partons clustered identically: recluster in the lab is not
!     equivalent, so map each Breit jet to the lab via the sum of the
!     lab momenta of its constituents (zkt returns the assignment)
      boostsgn = sign(1d0, plab(3,2))
      m = 0
      do i = 1, nj
         pj = jets(:,i)
         et2 = pj(1)**2 + pj(2)**2 + max(pj(0)**2 - pj(1)**2
     $        - pj(2)**2 - pj(3)**2, 0d0)
         if (et2 .le. 64d0) cycle
         call zlabjet(i, w)
         zel = w(0)
         pzl = boostsgn*w(3)
         pabs = sqrt(w(1)**2 + w(2)**2 + w(3)**2)
         if (pabs - abs(pzl) .le. 0d0) cycle
         etal = 0.5d0*log((pabs + pzl)/(pabs - pzl))
         if (etal .le. -1d0 .or. etal .ge. 2.5d0) cycle
         m = m + 1
         jets(:,m) = pj
         pt2(m) = pj(1)**2 + pj(2)**2
      enddo
      if (m .lt. 2) return
      i1 = maxloc(pt2(1:m), 1)
      pt2(i1) = -1d0
      i2 = maxloc(pt2(1:m), 1)
      pj = jets(:,i1) + jets(:,i2)
      m12 = sqrt(max(pj(0)**2 - pj(1)**2 - pj(2)**2 - pj(3)**2, 0d0))
      if (m12 .le. 20d0) return
      ptavg = (sqrt(jets(1,i1)**2 + jets(2,i1)**2)
     $     + sqrt(jets(1,i2)**2 + jets(2,i2)**2))/2d0
      call filld('sigma', 0.5d0, dsig)
      call filld('Q2', Q2, dsig)
      call filld('ptavg12', ptavg, dsig)
      call filld('m12', m12, dsig)
      end

!     inclusive kt (R = 1, E-scheme) on np partons; records which parton
!     ends in which jet (common /zktmap/) for the lab-frame momenta
      subroutine zkt(p, np, jets, nj)
      implicit none
      integer np, nj
      double precision p(0:3,4), jets(0:3,4)
      double precision zq(0:3,4), pt2(4), yr(4), ph(4), dmin, zd, zdy, dph
      integer act(4), m, i, j, ii, jj, own(4), k
      logical beam
      double precision, parameter :: pi = 3.141592653589793d0
      integer jetof(4), njm, npm
      common /zktmap/ jetof, njm, npm
      jetof = 0
      npm = np
      zq(:,1:np) = p(:,1:np)
      do i = 1, np
         act(i) = i
         own(i) = i
      enddo
      m = np
      nj = 0
      do while (m .gt. 0)
         do i = 1, m
            k = act(i)
            pt2(i) = zq(1,k)**2 + zq(2,k)**2
            yr(i) = 0.5d0*log(max(zq(0,k) + zq(3,k), 1d-300)
     $           /max(zq(0,k) - zq(3,k), 1d-300))
            ph(i) = atan2(zq(2,k), zq(1,k))
         enddo
!        beam distances first (ties to the beam): DISENT's counter-events
!        have partons exactly soft or collinear to the beam (p_T = 0),
!        for which d_iB = 0 = d_ij
         dmin = huge(1d0)
         ii = 0
         jj = 0
         beam = .true.
         do i = 1, m
            if (pt2(i) .lt. dmin) then
               dmin = pt2(i)
               ii = i
            endif
         enddo
         do i = 1, m
            do j = i + 1, m
               if (min(pt2(i), pt2(j)) .le. 0d0) cycle
               zdy = yr(i) - yr(j)
               dph = abs(ph(i) - ph(j))
               if (dph .gt. pi) dph = 2d0*pi - dph
               zd = min(pt2(i), pt2(j))*(zdy**2 + dph**2)
               if (zd .lt. dmin) then
                  dmin = zd
                  ii = i
                  jj = j
                  beam = .false.
               endif
            enddo
         enddo
         if (beam) then
            nj = nj + 1
            jets(:,nj) = zq(:,act(ii))
            do k = 1, np
               if (own(k) .eq. act(ii)) jetof(k) = nj
            enddo
            act(ii) = act(m)
            m = m - 1
         else
            zq(:,act(ii)) = zq(:,act(ii)) + zq(:,act(jj))
            do k = 1, np
               if (own(k) .eq. act(jj)) own(k) = act(ii)
            enddo
            act(jj) = act(m)
            m = m - 1
         endif
      enddo
      njm = nj
      end

!     lab-frame four-momentum of Breit jet i: sum of the lab momenta of
!     its partons
      subroutine zlabjet(i, w)
      use mod_analysis
      implicit none
      integer i, k
      double precision w(0:3)
      integer jetof(4), njm, npm
      common /zktmap/ jetof, njm, npm
      w = 0d0
!     only this event's partons (the map is reset in zkt)
      do k = 1, npm
         if (jetof(k) .eq. i) w = w + plab(:,3+k)
      enddo
      end
