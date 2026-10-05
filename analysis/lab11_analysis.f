!  Lab-frame 1+1 jet observables for the validation of tau_1 slicing
!  against disorder -p2b (slicing/nnlo11.f90 -integrated, obs11): no
!  FastJet. Lab frame with the proton along +z; partons with E <= 0 or
!  p_T = 0 (DISENT's zero momenta and exactly beam-collinear partons of
!  counter-configurations) are dropped; anti-k_t R = 1 (E scheme) with the
!  rapidity y; jets with p_T > 5 GeV and -1 < y < 2.5.
!  Histograms: sigma, >=1 jet, leading-jet p_T, leading-jet y, >=2 jets
!  (sigma per bin in the output for the counting histograms).

      subroutine define_histograms
      implicit none
      include 'pwhg_bookhist-multi.h'
      double precision, parameter :: pte(8) = (/5d0, 8d0, 11d0, 15d0,
     $     20d0, 30d0, 50d0, 100d0/)
      double precision, parameter :: ye(7) = (/-1d0, -0.5d0, 0d0, 0.5d0,
     $     1d0, 1.5d0, 2.5d0/)
      call bookupeqbins('sigma',1d0,0d0,1d0)
      call bookupeqbins('njet1',1d0,0d0,1d0)
      call bookup('ptlead',7,pte)
      call bookup('ylead',6,ye)
      call bookupeqbins('njet2',1d0,0d0,1d0)
      end

      subroutine user_analysis(n,dsig,x,y,Q2)
      use mod_parameters
      use mod_analysis
      implicit none
      integer n
      double precision dsig(maxscales), x, y, Q2
      double precision jv(0:3,4), pt(4), yr(4), ph(4), zs, dmin, dd
      double precision ptl, yl, dph
      logical alive(4), isjet(4)
      integer np, nq, i, j, ia, ja, nj5
      double precision, parameter :: zpi = 3.141592653589793d0
      np = n - 3
      zs = sign(1d0, plab(3,2))
      nq = 0
      do i = 1, np
         if (plab(0,3+i) .le. 0d0) cycle
         if (plab(1,3+i)**2 + plab(2,3+i)**2 .le.
     $        (1d-12*plab(0,3+i))**2) cycle
         nq = nq + 1
         jv(:,nq) = plab(:,3+i)
         jv(3,nq) = zs*jv(3,nq)
      enddo
      do i = 1, 4
         alive(i) = i .le. nq
         isjet(i) = .false.
      enddo
      do while (count(alive) .gt. 0)
         call kin
         dmin = huge(1d0)
         ia = 0
         ja = 0
         do i = 1, nq
            if (.not. alive(i)) cycle
            dd = 1d0/max(pt(i), 1d-300)**2
            if (dd .lt. dmin) then
               dmin = dd
               ia = i
               ja = 0
            endif
            do j = i + 1, nq
               if (.not. alive(j)) cycle
               dph = abs(ph(i) - ph(j))
               if (dph .gt. zpi) dph = 2d0*zpi - dph
               dd = min(1d0/max(pt(i), 1d-300)**2,
     $              1d0/max(pt(j), 1d-300)**2)
     $              *((yr(i) - yr(j))**2 + dph**2)
               if (dd .lt. dmin) then
                  dmin = dd
                  ia = i
                  ja = j
               endif
            enddo
         enddo
         if (ja .eq. 0) then
            alive(ia) = .false.
            isjet(ia) = .true.
         else
            jv(:,ia) = jv(:,ia) + jv(:,ja)
            alive(ja) = .false.
         endif
      enddo
      do i = 1, 4
         alive(i) = isjet(i)
      enddo
      call kin
      nj5 = 0
      ptl = -1d0
      yl = 0d0
      do i = 1, nq
         if (.not. isjet(i)) cycle
         if (pt(i) .gt. 5d0 .and. yr(i) .gt. -1d0 .and. yr(i) .lt. 2.5d0)
     $        then
            nj5 = nj5 + 1
            if (pt(i) .gt. ptl) then
               ptl = pt(i)
               yl = yr(i)
            endif
         endif
      enddo
      call filld('sigma', 0.5d0, dsig)
      if (nj5 .ge. 1) then
         call filld('njet1', 0.5d0, dsig)
         call filld('ptlead', ptl, dsig)
         call filld('ylead', yl, dsig)
      endif
      if (nj5 .ge. 2) call filld('njet2', 0.5d0, dsig)
      contains
      subroutine kin
      integer m
      do m = 1, nq
         if (.not. alive(m)) cycle
         pt(m) = sqrt(jv(1,m)**2 + jv(2,m)**2)
         yr(m) = 0.5d0*log((jv(0,m) + jv(3,m))
     $        /max(jv(0,m) - jv(3,m), 1d-300))
         ph(m) = atan2(jv(2,m), jv(1,m))
      enddo
      end subroutine kin
      end
