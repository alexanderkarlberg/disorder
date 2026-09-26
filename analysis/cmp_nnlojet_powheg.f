!----------------------------------------------------------------------
! disorder analysis for the comparisons with NNLOJET and POWHEG-BOX-RES
! (DIS). The observables are defined in cmp_obs_core.f, which the
! POWHEG-BOX analysis includes as well. Build with
!   cmake .. -DANALYSIS=cmp_nnlojet_powheg.f
!----------------------------------------------------------------------
      subroutine define_histograms
      implicit none
      call cmp_book
      end

      subroutine user_analysis(n,dsig,xdis,ydis,Q2dis)
      use mod_parameters
      use mod_analysis
      implicit none
      integer n
      double precision dsig(maxscales), xdis, ydis, Q2dis
      include 'cmp_obs_cuts.h'
      double precision kin(0:3), kout(0:3), pin(0:3), pout(0:3,4)
      double precision x, y, Q2, pj(0:3,4), protsign
      double precision tauzE, BzE, rhoE, cmp_eta
      integer npart, njets
      logical pass, valid

      kin  = plab(:,1)
      pin  = plab(:,2)
      kout = plab(:,3)
      npart = n - 3
      pout(:,1:npart) = plab(:,4:n)

      call cmp_kinematics(kin, kout, pin, S, x, y, Q2, pass)
      if (.not. pass) return

      call filld('sig', 0.5d0, dsig)
      call filld('Q2', Q2, dsig)
      call filld('x', x, dsig)
      call filld('y', y, dsig)

      protsign = sign(1d0, pin(3))
      call cmp_antikt(npart, pout(:,1:npart), cmp_Rjet, cmp_ptjmin,
     $     njets, pj)
      if (njets .ge. 1) then
         call filld('ptj1_lab', sqrt(pj(1,1)**2 + pj(2,1)**2), dsig)
         call filld('etaj1_lab', cmp_eta(pj(:,1), protsign), dsig)
      endif

      call cmp_event_shapes(kin, kout, pin, x, npart, pout(:,1:npart),
     $     tauzE, BzE, rhoE, valid)
      if (valid) then
         call filld('tauzE', tauzE, dsig)
         call filld('BzE', BzE, dsig)
         call filld('rhoE', rhoE, dsig)
      endif
      end

      include 'cmp_obs_core.f'
