!----------------------------------------------------------------------
! Observables for the disorder / NNLOJET / POWHEG-BOX-RES(DIS)
! comparisons. This file is included by the disorder analysis
! (analysis/cmp_nnlojet_powheg.f) and, as a copy, by the POWHEG-BOX
! analysis (pwhg_analysis_cmp.f), so that both codes compute exactly
! the same observables. The definitions follow NNLOJET
! (driver/core/EvalDIS.f90 in NNLOJET 1.0.2):
!   - event shapes in the Breit frame, current hemisphere = partons
!     whose longitudinal momentum points along q, normalised to the
!     sum of |p| in the current hemisphere (E-scheme):
!       tau_zE = 1 - sum|pz| / sum|p|           (dis_thrust)
!       B_zE   = sum pt / (2 sum|p|)            (dis_JB)
!       rho_E  = ((sum|p|)^2 - |sum p|^2)/(2 sum|p|)^2  (dis_JM2)
!     undefined (not filled) if the current hemisphere is empty;
!     these are not infrared safe beyond their leading order (a soft
!     gluon alone in an otherwise empty current hemisphere gives a
!     finite tau_zE, B_zE), so they are also booked with the cut
!     E_current > cmp_Ecfrac * Q on the energy in the current hemisphere
!     (names *_Ec; H1's definition; NNLOJET: dis_eventshapes = 0.1 in
!     the PROCESS block);
!   - lab-frame jets: anti-kt, R = 1, E-scheme recombination, in
!     (rapidity, azimuth), from the outgoing partons only; pseudo-
!     rapidities are given with the proton direction as positive.
!
! All momenta are (E,px,py,pz) in the lab frame, with either beam
! orientation.
!----------------------------------------------------------------------

!     Kinematics and cuts. Returns pass=.false. outside the cuts.
      subroutine cmp_kinematics(kin, kout, pin, s, x, y, Q2, pass)
      implicit none
      double precision kin(0:3), kout(0:3), pin(0:3), s, x, y, Q2
      logical pass
      double precision q(0:3), cmp_dot
      include 'cmp_obs_cuts.h'
      q = kin - kout
      Q2 = -cmp_dot(q,q)
      y = cmp_dot(pin,q) / cmp_dot(pin,kin)
      x = Q2 / (s * y)
      pass = Q2 .ge. cmp_Q2min .and. Q2 .le. cmp_Q2max .and.
     $       y .ge. cmp_ymin .and. y .le. cmp_ymax
      end

!     Event shapes in the Breit frame. valid=.false. if the current
!     hemisphere is empty.
      subroutine cmp_event_shapes(kin, kout, pin, x, npart, pout,
     $     tauzE, BzE, rhoE, valid)
      implicit none
      integer npart
      double precision kin(0:3), kout(0:3), pin(0:3), x
      double precision pout(0:3,npart), tauzE, BzE, rhoE
      logical valid
      double precision q(0:3), pb(0:3), qb(0:3), bmat(0:3,0:3)
      double precision sabsp, sabspz, spt, ptot(1:3), absp
      integer i, ncur

      q = kin - kout
      call cmp_breit_matrix(q, pin, x, bmat)
      qb = matmul(bmat, q)
      sabsp = 0d0; sabspz = 0d0; spt = 0d0; ptot = 0d0; ncur = 0
      do i = 1, npart
         pb = matmul(bmat, pout(:,i))
         if (pb(3)*qb(3) .gt. 0d0) then
            ncur = ncur + 1
            absp = sqrt(pb(1)**2 + pb(2)**2 + pb(3)**2)
            sabsp  = sabsp + absp
            sabspz = sabspz + abs(pb(3))
            spt    = spt + sqrt(pb(1)**2 + pb(2)**2)
            ptot   = ptot + pb(1:3)
         endif
      enddo
      valid = ncur .gt. 0
      if (.not. valid) return
      tauzE = 1d0 - sabspz / sabsp
      BzE   = spt / (2d0 * sabsp)
      rhoE  = (sabsp**2 - sum(ptot**2)) / (2d0 * sabsp)**2
      end

!     Energy in the current hemisphere of the Breit frame (as in
!     NNLOJET's dis_eventshapes cut, driver/core/ecuts.f).
      subroutine cmp_current_energy(kin, kout, pin, x, npart, pout,
     $     ecur)
      implicit none
      integer npart
      double precision kin(0:3), kout(0:3), pin(0:3), x
      double precision pout(0:3,npart), ecur
      double precision q(0:3), pb(0:3), qb(0:3), bmat(0:3,0:3)
      integer i
      q = kin - kout
      call cmp_breit_matrix(q, pin, x, bmat)
      qb = matmul(bmat, q)
      ecur = 0d0
      do i = 1, npart
         pb = matmul(bmat, pout(:,i))
         if (pb(3)*qb(3) .gt. 0d0) ecur = ecur + pb(0)
      enddo
      end

!     Lorentz transformation to the Breit frame: the rest frame of
!     2xP + q (P the proton momentum, along pin), rotated such that q
!     lies along the z axis. There q = (0,0,0,+Q) and the proton moves
!     along -z.
      subroutine cmp_breit_matrix(q, pin, x, bmat)
      implicit none
      double precision q(0:3), pin(0:3), x, bmat(0:3,0:3)
      double precision b(0:3), mb, g, bv(3), bb, lam(0:3,0:3)
      double precision qr(0:3), n(3), rot(0:3,0:3), cmp_dot
      double precision xP(0:3)
      integer i, j
!     xP along the incoming parton, normalised such that 2 xP.q = Q^2
      xP = pin * (-cmp_dot(q,q)) / (2d0 * cmp_dot(pin,q))
      b = 2d0 * xP + q
      mb = sqrt(cmp_dot(b,b))
!     pure boost to the rest frame of b
      g = b(0) / mb
      bv = b(1:3) / b(0)
      bb = sum(bv**2)
      lam = 0d0
      lam(0,0) = g
      lam(0,1:3) = -g * bv
      lam(1:3,0) = -g * bv
      do i = 1, 3
         do j = 1, 3
            lam(i,j) = (g - 1d0) * bv(i) * bv(j) / bb
         enddo
         lam(i,i) = lam(i,i) + 1d0
      enddo
      qr = matmul(lam, q)
!     rotation taking the direction of qr onto +z
      n = qr(1:3) / sqrt(sum(qr(1:3)**2))
      call cmp_rot_to_z(n, rot)
      bmat = matmul(rot, lam)
      end

!     Rotation matrix (acting on 4-vectors) taking unit vector n to +z
      subroutine cmp_rot_to_z(n, rot)
      implicit none
      double precision n(3), rot(0:3,0:3)
      double precision a(3), s, c, k(3), kk(3,3)
      integer i, j
      rot = 0d0
      rot(0,0) = 1d0
      c = n(3)
!     axis k = n x z, sin = |k|
      a = (/ n(2), -n(1), 0d0 /)
      s = sqrt(sum(a**2))
      if (s .lt. 1d-14) then
         do i = 1, 3
            rot(i,i) = sign(1d0, c)
         enddo
         if (c .lt. 0d0) rot(1,1) = 1d0 ! rotation by pi about x
         return
      endif
      k = a / s
      kk = 0d0
      kk(1,2) = -k(3); kk(1,3) =  k(2)
      kk(2,1) =  k(3); kk(2,3) = -k(1)
      kk(3,1) = -k(2); kk(3,2) =  k(1)
      do i = 1, 3
         do j = 1, 3
            rot(i,j) = s * kk(i,j) + (1d0 - c) * sum(kk(i,:)*kk(:,j))
         enddo
         rot(i,i) = rot(i,i) + 1d0
      enddo
      end

!     anti-kt (R, E-scheme) clustering of massless or massive
!     momenta in (rapidity, azimuth). Jets are returned ordered in
!     decreasing pt; only jets with pt > ptmin are kept.
      subroutine cmp_antikt(npart, p, R, ptmin, njets, pj)
      implicit none
      integer npart, njets
      double precision p(0:3,npart), R, ptmin, pj(0:3,npart)
      double precision w(0:3,npart), d, dmin, dij, pt2i, pt2j
      double precision cmp_rap, cmp_phi, drap, dphi, tmp(0:3)
      double precision pi
      parameter (pi = 3.141592653589793d0)
      logical alive(npart)
      integer i, j, imin, jmin, nw
      w = p
!     Partons with zero transverse momentum (DISENT passes exactly zero
!     momenta for unused slots in its O(as^2) collinear counterterms)
!     have infinite beam distance and cannot form a jet above ptmin:
!     leave them out, otherwise the loop below never terminates.
      alive = w(1,:)**2 + w(2,:)**2 .gt. 0d0
      nw = count(alive)
      njets = 0
      do while (nw .gt. 0)
         dmin = 1d300; imin = 0; jmin = 0
         do i = 1, npart
            if (.not. alive(i)) cycle
            pt2i = w(1,i)**2 + w(2,i)**2
            d = 1d0 / pt2i                ! beam distance
            if (d .lt. dmin) then
               dmin = d; imin = i; jmin = 0
            endif
            do j = i+1, npart
               if (.not. alive(j)) cycle
               pt2j = w(1,j)**2 + w(2,j)**2
               drap = cmp_rap(w(:,i)) - cmp_rap(w(:,j))
               dphi = abs(cmp_phi(w(:,i)) - cmp_phi(w(:,j)))
               if (dphi .gt. pi) dphi = 2d0*pi - dphi
               dij = min(1d0/pt2i, 1d0/pt2j) * (drap**2 + dphi**2)/R**2
               if (dij .lt. dmin) then
                  dmin = dij; imin = i; jmin = j
               endif
            enddo
         enddo
         if (jmin .eq. 0) then
            if (w(1,imin)**2 + w(2,imin)**2 .gt. ptmin**2) then
               njets = njets + 1
               pj(:,njets) = w(:,imin)
            endif
            alive(imin) = .false.
         else
            w(:,imin) = w(:,imin) + w(:,jmin)
            alive(jmin) = .false.
         endif
         nw = count(alive)
      enddo
!     order in pt
      do i = 1, njets
         do j = i+1, njets
            if (pj(1,j)**2+pj(2,j)**2 .gt. pj(1,i)**2+pj(2,i)**2) then
               tmp = pj(:,i); pj(:,i) = pj(:,j); pj(:,j) = tmp
            endif
         enddo
      enddo
      end

      double precision function cmp_dot(a, b)
      implicit none
      double precision a(0:3), b(0:3)
      cmp_dot = a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
      end

      double precision function cmp_rap(p)
      implicit none
      double precision p(0:3)
      cmp_rap = 0.5d0 * log((p(0) + p(3)) / (p(0) - p(3)))
      end

      double precision function cmp_phi(p)
      implicit none
      double precision p(0:3)
      cmp_phi = atan2(p(2), p(1))
      end

!     pseudorapidity, positive in the proton direction (protsign = +1
!     if the proton moves along +z in the lab frame, -1 otherwise)
      double precision function cmp_eta(p, protsign)
      implicit none
      double precision p(0:3), protsign, pabs, pz
      pabs = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
      pz = protsign * p(3)
      cmp_eta = 0.5d0 * log((pabs + pz) / (pabs - pz))
      end

!     Book the histograms (same names and binning in all codes)
      subroutine cmp_book
      implicit none
      include 'cmp_obs_cuts.h'
      double precision Q2bins(10), xbins(10), ptbins(10)
      data Q2bins /150d0, 200d0, 300d0, 500d0, 800d0, 1300d0, 2000d0,
     $     3500d0, 6000d0, 15000d0/
      data xbins /0.0015d0, 0.003d0, 0.006d0, 0.012d0, 0.025d0,
     $     0.05d0, 0.1d0, 0.2d0, 0.4d0, 1d0/
      data ptbins /5d0, 10d0, 15d0, 20d0, 30d0, 40d0, 50d0, 70d0,
     $     100d0, 150d0/
      call bookupeqbins('sig', 1d0, 0d0, 1d0)
      call bookup('Q2', 9, Q2bins)
      call bookup('x', 9, xbins)
      call bookupeqbins('y', 0.1d0, 0.1d0, 0.9d0)
      call bookup('ptj1_lab', 9, ptbins)
      call bookupeqbins('etaj1_lab', 0.5d0, -2d0, 3d0)
      call bookupeqbins('tauzE', 0.04d0, 0.02d0, 0.98d0)
      call bookupeqbins('BzE', 0.02d0, 0.02d0, 0.5d0)
      call bookupeqbins('rhoE', 0.01d0, 0.01d0, 0.25d0)
      call bookupeqbins('tauzE_Ec', 0.04d0, 0.02d0, 0.98d0)
      call bookupeqbins('BzE_Ec', 0.02d0, 0.02d0, 0.5d0)
      call bookupeqbins('rhoE_Ec', 0.01d0, 0.01d0, 0.25d0)
      end
