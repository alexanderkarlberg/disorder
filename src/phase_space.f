      module phase_space
      use parameters
      implicit none
      
      PRIVATE

      PUBLIC gen_phsp_born, gen_phsp_real, p2bmomenta, mbreit2lab,
     $     mlab2breit, mbreit2labdisent
     
      contains

      subroutine gen_phsp_born(xborn, x, y, Qsq, Qvec, jacobian, plab, pbreit)
      implicit none
      double precision x, y, Qsq, jacobian, plab(0:3,2+2), pbreit(0:3,2+2)
      double precision xborn(2), jac, Qval, Qvec(0:3)
      double precision pi
      parameter (pi = 3.14159265358979323846d0)

      jac = 1d0
      pbreit = 0d0
      plab = 0d0
      Qvec = 0d0

!     Generate first momentum fraction x
      if(xmin.eq.xmax) then
         x = xmin
         jac = 1d0 * jac
!      elseif(xmin.lt.ymax) then
      elseif(xmin.lt.xmax) then
         x = xmin * exp(xborn(1)*log(xmax/xmin))
         jac = jac * log(xmax/xmin) * x
      endif

!     Check available phase space for y/x
      ymn = ymin
      ymx = ymax

      if(ymin.lt. Q2min/(x * s)) ymn = Q2min/(x * s)
      if(ymax.gt. Q2max/(x * s)) ymx = Q2max/(x * s)

      if(ymn.eq.ymx) then
         y = ymn
         jac = 1d0 * jac
      elseif(ymn.lt.ymx) then
         y = ymn * exp(xborn(2)*log(ymx/ymn))
         jac = jac * log(ymx/ymn) * y
      else
         print*, ymn,ymx,x,s,Q2min,Q2max
         stop 'No phase space available for y'
      endif
! For numerical stability
      if(Q2min.eq.Q2max) then
         Qsq = Q2min
      else
         Qsq = y * x * s
         jac = jac * x * s
      endif
      Qval = sqrt(Qsq)
      jacobian = jac 

      if(y.lt.ymin.or.y.gt.ymax) jacobian = 0d0
      if(Qsq.lt.Q2min.or.Qsq.gt.Q2max) jacobian = 0d0
      if(x.lt.xmin.or.x.gt.xmax) jacobian = 0d0

      ! Now set up the momenta in the Breit frame first

      ! Incoming lepton
      pbreit(0,1) = (2.0d0 - y) * 0.5d0 / y
      pbreit(1,1) = - sqrt(1-y)/y
      pbreit(2,1) = 0d0
      pbreit(3,1) = - 0.5d0
      pbreit(:,1) = Qval * pbreit(:,1)
      ! Incoming parton
      pbreit(0,2) = Qval * 0.5d0
      pbreit(3,2) = -Qval * 0.5d0
      ! Outgoing lepton
      pbreit(0,3) = (2.0d0 - y) * 0.5d0 / y
      pbreit(1,3) = - sqrt(1-y)/y
      pbreit(2,3) = 0d0
      pbreit(3,3) = 0.5d0
      pbreit(:,3) = Qval * pbreit(:,3)
      ! Outgoing parton
      pbreit(0,4) = Qval * 0.5d0
      pbreit(3,4) = Qval * 0.5d0 

      Qvec(0) = y * (El - x * Eh)
      Qvec(1) = Qval * sqrt(1 - y) 
      Qvec(2) = 0d0
      Qvec(3) = y * (El + x * Eh)

!     In lab frame by hand
      ! Incoming lepton
      plab(0,1) = El
      plab(3,1) = El
      ! Incoming parton
      plab(0,2) = x * Eh
      plab(3,2) = -x * Eh
      ! Outgoing lepton
      plab(:,3) = plab(:,1) - Qvec(:)
      ! Outgoing parton
      plab(:,4) = plab(:,2) + Qvec

      !print*, 'Momenta in Breit frame'
      !print*, 'Incoming lepton:', pbreit(:,1)
      !print*, 'Incoming parton:', pbreit(:,2)
      !print*, 'Outgoing lepton:', pbreit(:,3)
      !print*, 'Outgoing parton:', pbreit(:,4)
      !print*, 'Momenta in lab frame'
      !print*, 'Incoming lepton:', plab(:,1)
      !print*, 'Incoming parton:', plab(:,2)
      !print*, 'Outgoing lepton:', plab(:,3)
      !print*, 'Outgoing parton:', plab(:,4)
      !print*, 'Qvec:           ', Qvec(:)
      !print*, 'Qsq x 2:        ', Qsq, invmsq(Qvec)
      !print*, 'x,y:            ', x, y


      end subroutine gen_phsp_born

      ! It looks like something is bot right in this phase space. In
      ! particular if one plots the current Energy (defined as the sum
      ! of the energy of all partons in the direction of Q in the Breit
      ! frame) this is clearly not bounded by Q as it should be. 
      subroutine gen_phsp_real(xreal, x, y, Qsq, Qvec, csi, z, jacobian, plab, pbreit)
      implicit none
      double precision x, y, Qsq, jacobian, plab(0:3,2+3), pbreit(0:3,2+3)
      double precision xreal(2), jac, Qval, Qvec(0:3)
      ! The two variables of hep-ph/9704297
      double precision csi, z, z0, z3, z0bar, z3bar, zperp, oocsi
      DOUBLE PRECISION pbreit0(0:3,3), plab0(0:3,3)
      
      pbreit = 0d0
      plab = 0d0
      pbreit0 = 0d0
      Qval = sqrt(Qsq)
      ! We start by setting up the leptons as they are not affected by the emission.
      ! This routine takes x,y,Qsq and Qvec as inputs
      ! Incoming lepton
      pbreit(0,1) = (2.0d0 - y) * 0.5d0 / y
      pbreit(1,1) = - sqrt(1-y)/y
      pbreit(2,1) = 0d0
      pbreit(3,1) = - 0.5d0
      pbreit(:,1) = Qval * pbreit(:,1)
      ! Outgoing lepton
      pbreit(0,3) = (2.0d0 - y) * 0.5d0 / y
      pbreit(1,3) = - sqrt(1-y)/y
      pbreit(2,3) = 0d0
      pbreit(3,3) = 0.5d0
      pbreit(:,3) = Qval * pbreit(:,3)
      ! Incoming lepton
      plab(0,1) = El
      plab(3,1) = El
      ! Outgoing lepton
      plab(:,3) = plab(:,1) - Qvec(:)

      jacobian = 1d0

      ! csi is in range x < csi < 1 and 0 < z < 1

      csi = xreal(1) * (1d0 - x) + x
      jacobian = jacobian * (1d0 - x)
      z = xreal(2)

      ! eq. 3.6
      oocsi = 1d0 / csi
      z0 = 2d0 * z - 1d0 + (1d0 - z) * oocsi
      z3 = 1d0 - (1d0 - z) * oocsi
      z0bar = 1d0 - 2d0 * z + z * oocsi
      z3bar = 1d0 - z * oocsi
      zperp = 2d0 * sqrt(z * (1d0 - z)*(1d0 - csi) * oocsi)

      ! Now we construct the partons (incoming and outgoing) in the Breit frame from eq. 3.5
      ! Incoming parton
      pbreit0(0,1) = oocsi
      pbreit0(3,1) = - oocsi
      ! First outgoing parton
      pbreit0(0,2) = z0
      pbreit0(1,2) = zperp
      pbreit0(3,2) = z3
      ! Second outgoing parton
      pbreit0(0,3) = z0bar
      pbreit0(1,3) = -zperp
      pbreit0(3,3) = z3bar

      pbreit0 = Qval * 0.5d0 * pbreit0

      ! Now transfer these to the right array
      pbreit(:,2) = pbreit0(:,1)
      pbreit(:,4) = pbreit0(:,2)
      pbreit(:,5) = pbreit0(:,3)

      call mbreit2lab(3,Qvec,pbreit0,plab0,.false.)

      ! Now transfer these to the right array
      plab(:,2) = plab0(:,1)
      plab(:,4) = plab0(:,2)
      plab(:,5) = plab0(:,3)


      !print*, 'Momenta in Breit frame'
      !print*, 'Incoming lepton :', pbreit(:,1)
      !print*, 'Incoming parton :', pbreit(:,2)
      !print*, 'Outgoing lepton :', pbreit(:,3)
      !print*, 'Outgoing parton1:', pbreit(:,4)
      !print*, 'Outgoing parton2:', pbreit(:,5)
      !print*, 'Momenta in lab frame'
      !print*, 'Incoming lepton :', plab(:,1)
      !print*, 'Incoming parton :', plab(:,2)
      !print*, 'Outgoing lepton :', plab(:,3)
      !print*, 'Outgoing parton1:', plab(:,4)
      !print*, 'Outgoing parton1:', plab(:,5)
      !print*, 'Qvec:           ', Qvec(:)
      !print*, 'Qsq x 2:        ', Qsq, invmsq(Qvec)
      !print*, 'x,y:            ', x, y

      end subroutine gen_phsp_real


      double precision function invmsq(p) 
      implicit none
      double precision p(0:3)

      invmsq = p(0)**2 - p(1)**2 - p(2)**2 - p(3)**2
      end function invmsq

!     Q in the lab frame (normal DIS variable), p in the lab frame and
!     pout in the breit frame. This uses appendix seven of the DIS book
!     by A Cooper-Sakar. Page 208 Appendix 7.11. 
      subroutine mlab2breit(m,Q,p,pout,isPlus)
      implicit none
      integer i,j,m
      logical isPlus
      double precision Q(0:3), p(0:3,m), pin(0:3,m), pout(0:3,m), Qval
      double precision cosphi, sinphi, invQval, norm, q0, q1, q2, q3
      double precision z(3), Qb(0:3), modrot
      parameter (z = (/0d0,0d0,1d0/))
      double precision trans(0:3,0:3)

      Qval = DSqrt(-invmsq(Q))
      invQval = 1d0/Qval
      Qb = Q                    ! Local copy
      pin = p                   ! Local copy
!     The routine assumes "+" like incoming direction. If that is not
!     the case, we revert the z-component.
      q0 = Qb(0)
      q1 = Qb(1)
      q2 = Qb(2)
      q3 = Qb(3)
      if(.not.isPlus) then
         q3 = -Qb(3)
      endif
      norm = 1d0/(q0 - q3)


!     Compute composite Lorentz transformation
      trans(0,0) = q0*invQval + Qval*norm
      trans(0,1) = -q1*invQval
      trans(0,2) = -q2*invQval
      trans(0,3) = -q3*invQval -Qval*norm

      trans(1,0) = -q1*norm
      trans(1,1) = 1d0
      trans(1,2) = 0d0
      trans(1,3) = q1*norm

      trans(2,0) = -q2*norm
      trans(2,1) = 0d0
      trans(2,2) = 1d0
      trans(2,3) = q2*norm

      trans(3,0) = q0*invQval
      trans(3,1) = -q1*invQval
      trans(3,2) = -q2*invQval
      trans(3,3) = -q3*invQval
      
      if (.not.isPlus) then
         trans(:,3) = - trans(:,3)
         trans(3,:) = - trans(3,:)
      endif
      
!     Perform transformation
      do j=1,m
         pout(:,j) = MATMUL(trans,pin(:,j))
      enddo
      end

!     Q in the lab frame (normal DIS variable), p in the breit frame and
!     pout in the lab frame. This uses appendix seven of the DIS book by
!     A Cooper-Sakar. Page 208 Appendix 7.11. This is the inverse of the
!     transformation mlab2breit2.
      subroutine mbreit2lab(m,Q,p,pout,isPlus)
      implicit none
      integer i,j,m
      logical isPlus
      double precision Q(0:3), p(0:3,m), pin(0:3,m), pout(0:3,m), Qval
      double precision cosphi, sinphi, norm, invQval
      double precision z(3), Qb(0:3), modrot, q0, q1, q2, q3
      parameter (z = (/0d0,0d0,1d0/))
      double precision trans(0:3,0:3)

      Qval = DSqrt(-invmsq(Q))
      invQval = 1d0/Qval
      Qb = Q                    ! Local copy
      pin = p                   ! Local copy
!     The routine assumes "+" like incoming direction. If that is not
!     the case, we revert the z-component.
      q0 = Qb(0)
      q1 = Qb(1)
      q2 = Qb(2)
      q3 = Qb(3)
      if(.not.isPlus) then
         q3 = -Qb(3)
      endif
      norm = 1d0/(q0 - q3)


!     Compute composite Lorentz transformation
      trans(0,0) = q0*invQval + Qval*norm
      trans(0,1) = q1*norm
      trans(0,2) = q2*norm
      trans(0,3) = -q0*invQval

      trans(1,0) = q1*invQval
      trans(1,1) = 1d0
      trans(1,2) = 0d0
      trans(1,3) = -q1*invQval

      trans(2,0) = q2*invQval
      trans(2,1) = 0d0
      trans(2,2) = 1d0
      trans(2,3) = -q2*invQval

      trans(3,0) = Qval*norm + q3*invQval
      trans(3,1) = q1*norm
      trans(3,2) = q2*norm
      trans(3,3) = -q3*invQval
      
      if (.not.isPlus) then
         trans(:,3) = - trans(:,3)
         trans(3,:) = - trans(3,:)
      endif
      
!     Perform transformation
      do j=1,m
         pout(:,j) = MATMUL(trans,pin(:,j))
      enddo
      end
      
      ! This is not too nice. The routines above do not for disent as
      ! they give the wrong sign for the quarks due to different
      ! conventions, and the leptons have to be set by hand. For now we
      ! fix it by hand, but perhaps worth thinking of making the
      ! conventions uniform.
      subroutine mbreit2labdisent(m,Q,p,pout)
      implicit none
      integer m
      double precision Q(0:3), p(0:3,m), pin(0:3,m), pout(0:3,m),
     $     plcl(0:3,m)

      plcl = 0d0
      pin = 0d0
      pout = 0d0
      
!     First we do partons
      pin(:,1) = p(:,2)
      pin(:,2:m-2) = p(:,4:m)
      pin(3,1) = -pin(3,1)
      pin(3,2:m-2) = - pin(3,2:m-2)

      call mbreit2lab(m-2,Q,pin,plcl,.false.)
      pout(:,2) = plcl(:,1)
      pout(:,4:m) = plcl(:,2:m-2)
      
!     Then the leptons are trivial. We fix the incoming and use the Q
!     vector to get the outgoing one.
      pout(0,1) = El
      pout(3,1) = El
      pout(:,3) = pout(:,1) - Q(:)
      
      end
      
      subroutine p2bmomenta(x,y,Qsq,pbreit,plab)
      implicit none
      double precision x, y, Qsq, plab(0:3,2+2), pbreit(0:3,2+2), Qval,
     $     Qvec(0:3)
      
      pbreit = 0d0
      plab = 0d0
      Qvec = 0d0
      Qval = sqrt(Qsq)
      
!     Incoming lepton
      pbreit(0,1) = (2.0d0 - y) * 0.5d0 / y
      pbreit(1,1) = - sqrt(1-y)/y
      pbreit(2,1) = 0d0
      pbreit(3,1) = - 0.5d0
      pbreit(:,1) = Qval * pbreit(:,1)
      ! Incoming parton
      pbreit(0,2) = Qval * 0.5d0
      pbreit(3,2) = Qval * 0.5d0
      ! Outgoing lepton
      pbreit(0,3) = (2.0d0 - y) * 0.5d0 / y
      pbreit(1,3) = - sqrt(1-y)/y
      pbreit(2,3) = 0d0
      pbreit(3,3) = 0.5d0
      pbreit(:,3) = Qval * pbreit(:,3)
      ! Outgoing parton
      pbreit(0,4) = Qval * 0.5d0
      pbreit(3,4) = -Qval * 0.5d0

      Qvec(0) = y * (El - x * Eh)
      Qvec(1) = Qval * sqrt(1 - y) 
      Qvec(2) = 0d0
      Qvec(3) = y * (El + x * Eh)

!     In lab frame by hand
      ! Incoming lepton
      plab(0,1) = El
      plab(3,1) = El
      ! Incoming parton
      plab(0,2) = x * Eh
      plab(3,2) = -x * Eh
      ! Outgoing lepton
      plab(:,3) = plab(:,1) - Qvec(:)
      ! Outgoing parton
      plab(:,4) = plab(:,2) + Qvec

      end subroutine

      
      end module phase_space
