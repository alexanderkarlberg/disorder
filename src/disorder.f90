program disorder
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use sub_defs_io
  !use dummy_pdfs
  use matrix_element
  use parameters
  use phase_space
  use integration
  use types
  use mod_dsigma
  use structure_functions_gluon_only
  implicit none
  integer, parameter :: ndim = 7 ! dimension for vegas integration
  integer, parameter :: nmempdfMAX = 200 ! max number of pdfs
  integer  :: nmempdf_start, nmempdf_end, imempdf
  real(dp) :: integ, error_int, proba, tini, tfin
  real(dp) :: sigma_tot, error_tot, region(1:2*ndim)
  real(dp) :: res(0:nmempdfMAX), central, errminus, errplus, errsymm
  real(dp) :: res_scales(1:maxscales),maxscale,minscale
  integer :: iscales
  character * 100 :: analysis_name, histo_name, scale_string
  character * 3 :: pdf_string

  ! DISENT common block
  integer iev
  double precision weights(0:2)
  
  call cpu_time(tini)

  ! set up all constants and parameters from command line arguments
  call set_parameters()
  write(*,'(a)') '# Command line : '//trim(command_line())
  ! Initialise histograms
  call init_histo()
  ! Initialise PDF
  call initPDFSetByName(pdfname)

  if(run_disent.and.(.not.p2b)) then
     print*, 'Calling DISENT with'
     print*, 'NF  = ', nflav
     print*, 'NEV = ', ncall2

     call InitPDF(0)
     
     iev=0
     weights=0.0d0

     ! Careful, it seems disent runs with μR = μF
!     if(xmur.ne.xmuf.and.(order_max.ge.3)) stop 'With DISENT we must run with μR = μF'
     call DISENT(ncall2,S,nflav,user,dis_cuts,12345+iseed-1,67890+iseed-1,order_max-1,xmuf**2)
     call finalise_histograms('disent')

     weights = weights / iev

     print*, 'Born σ-2jet = ', weights(0), ' pb'
     print*, 'NLO  σ-2jet = ', weights(0) + weights(1), ' pb'
!     print*, 'Born σ-3jet = ', weights(1), ' pb'
!     print*, 'NLO  σ-3jet = ', weights(1) + weights(2), ' pb'
     !     print*, 'Born σ-4jet = ', weights(2), ' pb'
     call exit()
  elseif(run_disent.and.p2b) then
     call InitPDF(0)
     
     iev=0
     weights=0.0d0

     imempdf = nmempdf
     nmempdf_start = nmempdf
     
     write(6,*) "PDF member:",imempdf
     call InitPDF(imempdf)
     call getQ2min(0,Qmin)
     Qmin = sqrt(Qmin)
     print*, 'Qmin', Qmin
     
     if(Qmin.gt.sqrt(Q2min)) then
        print*, 'WARNING: PDF Qmin = ', Qmin
        print*, 'But running with input value of: ', sqrt(Q2min)
        stop
     endif
     
     do iscales = 1,Nscales
        call StartStrFct(sqrts*max(scales_mur(iscales),scales_muf(iscales)), order_max, nflav, xmur*scales_mur(iscales), &
             & xmuf*scales_muf(iscales), scale_choice, mz, .true., Qmin, mw, mz)
        call read_PDF()
        call InitStrFct(order_max, separate_orders = .true.)
        
        region(1:ndim)        = 0.0_dp
        region(ndim+1:2*ndim) = 1.0_dp
        sigma_tot = 0.0_dp
        error_tot = 0.0_dp
        
        
        ! vegas warmup call
        if(ncall1.gt.0.and.itmx1.gt.0) then 
           ! Skip grid generation if grid is being read from file
           writeout=.true.
           fillplots = .false.
           call vegas(region,ndim,dsigma,0,ncall1,itmx1,0,integ,error_int,proba)
           writeout=.false.
        ! set random seed to current idum value
           saveseed = idum
        else
           ! if reading in grids from first loop iteration, make sure
           ! saveseed is initialized to correct value
           saveseed = iseed
        endif
        
        ! Return if we are doing warm-up
        if(ncall2.lt.1) return
        if(itmx2.lt.1) return

        if(order_max.gt.1) then
           ! Then do the disent run
!           call DISENT(ncall2,S,nflav,user,dis_cuts,12345+iseed-1,67890+iseed-1,order_max-1,(xmuf*scales_muf(iscales))**2)
           call DISENTEXTENDED(ncall2,S,nflav,user,dis_cuts,12345&
                &+iseed-1,67890+iseed-1,order_max-1,(xmuf*scales_muf(iscales))**2,cflcl&
                &,calcl,trlcl,new_scaleuncert)
           ! Store disent result
           call pwhgaddout
        endif
        
        ! The inclusive code has much faster convergence, so we reduce
        ! the number of calls here
        ncall2 = ncall2/10
        
        ! vegas main call
        ! set random seed to saved value 
        idum     = -saveseed
        fillplots = .true.
        call vegas(region,ndim,dsigma,1,ncall2,itmx2,0,integ,error_int,proba)
        ! reset ncall2
        ncall2 = ncall2 * 10
        readin = .true.
        ! add integral to the total cross section
        sigma_tot = sigma_tot + integ
        error_tot = error_tot + error_int**2
        
        res(imempdf) = sigma_tot
        res_scales(iscales) = sigma_tot

        ! print total cross section and error into file  
        ! construct name of output file
        call pwhgaddout
        write(pdf_string,"(I3.3)") imempdf
        scale_string = scalestr(iscales)//"pdfmem_"//trim(pdf_string)
        if (order_max.eq.1) then
           histo_name="disorder_lo_seed"//seedstr//scale_string//".dat"
        else if (order_max.eq.2) then
           histo_name="disorder_nlo_seed"//seedstr//scale_string//".dat"
        else if (order_max.eq.3) then
           histo_name="disorder_nnlo_seed"//seedstr//scale_string//".dat"
        else if (order_max.eq.4) then
           histo_name="disorder_n3lo_seed"//seedstr//scale_string//".dat"
        endif
        print*, 'Outputting: ', histo_name
        call pwhgtopout(histo_name)
        call resethists
!        call finalise_histograms('disorder')
     enddo
     call exit()
  endif

  if (pdfuncert) then
     nmempdf_start = 0
     call numberPDF(nmempdf_end)
  else
     nmempdf_start = nmempdf
     nmempdf_end   = nmempdf
  endif

  if (nmempdf_end .gt. nmempdfMAX) stop "ERROR: increase nmempdfMAX"

  do imempdf = nmempdf_start,nmempdf_end
     write(6,*) "PDF member:",imempdf
     call InitPDF(imempdf)
     call getQ2min(0,Qmin)
     Qmin = sqrt(Qmin)
     print*, 'Qmin', Qmin
     
     if(Qmin.gt.sqrt(Q2min)) then
        print*, 'WARNING: PDF Qmin = ', Qmin
        print*, 'But running with input value of: ', sqrt(Q2min)
        stop
     endif
     if(new_scaleuncert) then
        call StartStrFct(sqrts, order_max, nflav, one, &
             & one, scale_choice, mz, .true., Qmin, mw, mz)
        call read_PDF()
        call InitStrFct(order_max, separate_orders = .true.)

        region(1:ndim)        = 0.0_dp
        region(ndim+1:2*ndim) = 1.0_dp
        sigma_tot = 0.0_dp
        error_tot = 0.0_dp
        
        
        ! vegas warmup call
        if(ncall1.gt.0.and.itmx1.gt.0) then 
           ! Skip grid generation if grid is being read from file
           writeout=.true.
           fillplots = .false.
           scaleuncert = .false.
           call vegas(region,ndim,dsigma,0,ncall1,itmx1,0,integ,error_int,proba)
           scaleuncert = .true.
           writeout=.false.
           ! set random seed to current idum value
           saveseed = idum
        elseif (imempdf.eq.nmempdf_start) then
           ! if reading in grids from first loop iteration, make sure
           ! saveseed is initialized to correct value
           saveseed = iseed
        endif

        if(ncall2.lt.1) return
        if(itmx2.lt.1) return
        ! vegas main call
        ! set random seed to saved value 
        idum     = -saveseed
        fillplots = .true.
        call vegas(region,ndim,dsigma,1,ncall2,itmx2,0,integ,error_int,proba)
        readin = .true.
        ! add integral to the total cross section
        sigma_tot = sigma_tot + integ
        error_tot = error_tot + error_int**2
        
        res(imempdf) = integ
        res_scales(1:nscales) = sigma_all_scales(1:nscales) 
        print*, 'sigma_tot', integ
        print*, 'res_scales', res_scales(1:nscales)
        
!        stop 'GOOD'
     else
     do iscales = 1,Nscales
        ! initialise hoppet, for now with fixed scale_choice = 1 (meaning photon virtuality)
!        call StartStrFct(sqrts, order_max, nflav, xmur*scales_mur(iscales), &
!             & xmuf*scales_muf(iscales), scale_choice, mh, .true., Qmin, mw, mz)
        call StartStrFct(sqrts, order_max, nflav, xmur*scales_mur(iscales), &
             & xmuf*scales_muf(iscales), scale_choice, mz, .true., Qmin, mw, mz)
        call read_PDF()
        call InitStrFct(order_max, separate_orders = .true.)
      if(gluon_only) then
        call StartStrFct_gluon_only(sqrts, order_max, nflav, xmur*scales_mur(iscales), &
             & xmuf*scales_muf(iscales), scale_choice, mz, .true., Qmin, mw, mz)
        call read_PDF()
        call InitStrFct_gluon_only(order_max, separate_orders = .true.)
      endif  
        
        region(1:ndim)        = 0.0_dp
        region(ndim+1:2*ndim) = 1.0_dp
        sigma_tot = 0.0_dp
        error_tot = 0.0_dp
        
        
        ! vegas warmup call
        if(ncall1.gt.0.and.itmx1.gt.0) then 
           ! Skip grid generation if grid is being read from file
           writeout=.true.
           fillplots = .false.
           call vegas(region,ndim,dsigma,0,ncall1,itmx1,0,integ,error_int,proba)
           writeout=.false.
           ! set random seed to current idum value
           saveseed = idum
        elseif (imempdf.eq.nmempdf_start) then
           ! if reading in grids from first loop iteration, make sure
           ! saveseed is initialized to correct value
           saveseed = iseed
        endif

        if(ncall2.lt.1) return
        if(itmx2.lt.1) return
        ! vegas main call
        ! set random seed to saved value 
        idum     = -saveseed
        fillplots = .true.
        call vegas(region,ndim,dsigma,1,ncall2,itmx2,0,integ,error_int,proba)
        readin = .true.
        ! add integral to the total cross section
        sigma_tot = sigma_tot + integ
        error_tot = error_tot + error_int**2
        
        res(imempdf) = sigma_tot
        res_scales(iscales) = sigma_tot
        
        ! print total cross section and error into file  
        ! construct name of output file
        call pwhgsetout
        write(pdf_string,"(I3.3)") imempdf
        scale_string = scalestr(iscales)//"pdfmem_"//trim(pdf_string)
        if (order_max.eq.1) then
           histo_name="histo_lo_seed"//seedstr//scale_string//".dat"
        else if (order_max.eq.2) then
           histo_name="histo_nlo_seed"//seedstr//scale_string//".dat"
        else if (order_max.eq.3) then
           histo_name="histo_nnlo_seed"//seedstr//scale_string//".dat"
        else if (order_max.eq.4) then
           histo_name="histo_n3lo_seed"//seedstr//scale_string//".dat"
        endif
        print*, 'Outputting: ', histo_name
        call pwhgtopout(histo_name)
        call resethists
  
     enddo
     endif
     if(imempdf.eq.nmempdf_start) then ! First PDF, this is where we compute scale uncertainties
        maxscale = maxval(res_scales(1:Nscales)) 
        minscale = minval(res_scales(1:Nscales))
        Nscales = 1
        res(imempdf) = res_scales(1) ! Copy central scale
     endif
     ! end loop over pdfs
 enddo

 analysis_name='xsct'
 if (order_max.eq.1) then
    analysis_name="xsct_lo_seed"//seedstr//".dat"
 else if (order_max.eq.2) then
    analysis_name="xsct_nlo_seed"//seedstr//".dat"
 else if (order_max.eq.3) then
    analysis_name="xsct_nnlo_seed"//seedstr//".dat"
 else if (order_max.eq.4) then
    analysis_name="xsct_n3lo_seed"//seedstr//".dat"
 endif
 !------------------------------------------------------------
  ! output results
  OPEN(UNIT=11, FILE=analysis_name, ACTION="write")
  write(11,'(a)',advance='no') '# '
  call time_stamp(11)
  write(11,'(a)') '# '//trim(command_line())
  if (order_max.eq.1) then 
     write(11,'(a,es13.6,a,es13.6,a)') '# Total LO cross-section (pb)'
  else if (order_max.eq.2) then 
     write(11,'(a,es13.6,a,es13.6,a)') '# Total NLO cross-section (pb)'
  else if (order_max.eq.3) then 
     write(11,'(a,es13.6,a,es13.6,a)') '# Total NNLO cross-section (pb)'
  else if (order_max.eq.4) then 
     write(11,'(a)') '# Total N3LO cross-section (pb)'
  endif

  if (nmempdf_start.eq.nmempdf_end.and..not.scaleuncert3.and..not.scaleuncert7.and..not.scaleuncert9) then
     write(11,'(a)') '# central     MC_error'
     write(11,'(es13.6,es13.6)') sigma_tot, sqrt(error_tot)
  elseif(nmempdf_start.eq.nmempdf_end) then
     central = res(0)
     write(11,'(a)') '# central     max          min          MC_error'
     write(11,'(7(es13.6))') central, maxscale, minscale,sqrt(error_tot)
     write(11,'(a)') ''
     write(11,'(a)') ''
     write(11,'(a)') '# Summary:'
     write(11,'(a,f16.5,a)') '# sigma =', central,' pb'
     write(11,'(a,f16.5,a,a,f9.3,a)') '# QCD scale uncertainty (+) =', maxscale-central, ' pb', &
          & ' (', ((maxscale-central)/central)*100.0_dp, ' %)'
     write(11,'(a,f16.5,a,a,f9.3,a)') '# QCD scale uncertainty (-) =', minscale-central, ' pb', &
          & ' (', ((minscale-central)/central)*100.0_dp, ' %)'
     write(11,'(a,f10.3,a)') '# MC integration uncertainty =', sqrt(error_tot)/central*100.0_dp, ' %'
  elseif(.not.scaleuncert3.and..not.scaleuncert7.and..not.scaleuncert9) then
     call getpdfuncertainty(res(nmempdf_start:nmempdf_end),central,errplus,errminus,errsymm)
     write(11,'(a)') '# central     max          min          MC_error'
     write(11,'(7(es13.6))') central, maxscale, minscale,sqrt(error_tot)
     write(11,'(a)') ''
     write(11,'(a)') ''
     write(11,'(a)') '# Summary:'
     write(11,'(a,f10.5,a)') '# sigma =', central,' pb'
     write(11,'(a,f10.3,a)') '# MC integration uncertainty =', sqrt(error_tot)/central*100.0_dp, ' %'
     write(11,'(a,f10.3,a)') '# PDF symmetric uncertainty* =', errsymm/central*100.0_dp, ' %'
     write(11,'(a)') '# (*PDF uncertainty contains alphas uncertainty if using a'
     write(11,'(a)') '#   PDF set that supports it (eg PDF4LHC15_nnlo_100_pdfas))'
  else
     call getpdfuncertainty(res(nmempdf_start:nmempdf_end),central,errplus,errminus,errsymm)
     write(11,'(a)') '# central     max          min          MC_error     PDF_err_plus PDF_err_min  PDF_err_symm'
     write(11,'(7(es13.6))') central, maxscale, minscale,sqrt(error_tot), errplus, errminus, errsymm
     write(11,'(a)') ''
     write(11,'(a)') ''
     write(11,'(a)') '# Summary:'
     write(11,'(a,f10.5,a)') '# sigma =', central,' pb'
     write(11,'(a,f9.5,a,a,f9.3,a)') '# QCD scale uncertainty (+) =', maxscale-central, ' pb', &
          & ' (', ((maxscale-central)/central)*100.0_dp, ' %)'
     write(11,'(a,f9.5,a,a,f9.3,a)') '# QCD scale uncertainty (-) =', minscale-central, ' pb', &
          & ' (', ((minscale-central)/central)*100.0_dp, ' %)'
     write(11,'(a,f10.3,a)') '# MC integration uncertainty =', sqrt(error_tot)/central*100.0_dp, ' %'
     write(11,'(a,f10.3,a)') '# PDF symmetric uncertainty* =', errsymm/central*100.0_dp, ' %'
     write(11,'(a)') '# (*PDF uncertainty contains alphas uncertainty if using a'
     write(11,'(a)') '#   PDF set that supports it (eg PDF4LHC15_nnlo_100_pdfas))'
  endif

  call cpu_time(tfin)
  write(6,'(a)')
  write(6,'(a,es9.2,a)') '==================== TOTAL TIME : ', &
       &     tfin-tini, ' s.'
  write(6,'(a)')

contains

!----------------------------------------------------------------------
  ! fill the streamlined interface PDF table (possibly using hoppet's
  ! evolution)
  subroutine read_PDF()
    use streamlined_interface
    real(dp), external :: alphasPDF
    interface
       subroutine EvolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine EvolvePDF
    end interface
    !----------------
    real(dp) :: res_lhapdf(-6:6), x, Q
    real(dp) :: res_hoppet(-6:6)
    real(dp), parameter :: mz = 91.2_dp

   ! InitRunningCoupling has to be called for the HOPPET coupling to be initialised 
   ! Default is to ask for 4 loop running and threshold corrections at quark masses.  
   call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , 4,&
        & -1000000045, quark_masses_sf(4:6), .true.)
   ! fixnf can be set to a positive number for
   ! fixed nf. -1000000045 gives variable nf
   ! and threshold corrections at quarkmasses.
   call hoppetAssign(EvolvePDF)

 end subroutine read_PDF

 subroutine dis_cuts(s,xminl,xmaxl,Q2minl,Q2maxl,yminl,ymaxl)
   real(dp), intent(in) :: s
   real(dp), intent(out) :: xminl,xmaxl,Q2minl,Q2maxl,yminl,ymaxl
   
   xminl = xmin
   xmaxl = xmax
   
   Q2minl = Q2min
   Q2maxl = Q2max
   
   yminl = 0d0
   ymaxl = 1d0
  
 end subroutine dis_cuts

 ! user-defined event analysis for disent
 subroutine user(N,NA,NT,P,S,WEIGHT)
   implicit none
   integer, intent(in) :: N, NA, NT
   real(dp), intent(in) :: s, p(4,7), weight(-6:6)
   LOGICAL SCALE_VAR
   DOUBLE PRECISION SCL_WEIGHT(3,-6:6)
   COMMON/cSCALE_VAR/SCL_WEIGHT, SCALE_VAR

   integer i

   double precision dsig, alphasPDF, totwgt
   double precision, save ::  pdfs(-6:6), eta, Q2, Qval, x, y, as2pi

   double precision, save :: etasave = -100d0
   ! For p2b
   logical, save :: recompute = .true.
   double precision, save ::  p2blab(0:3,2+2), p2bbreit(0:3,2+2), Qlab(0:3)

   if(order_max.le.2.and.NA.ge.2) return ! Disregard O(αS**2) if we are doing NLO
   if(order_max.le.1.and.NA.ge.1) return ! Disregard O(αS) if we are doing LO

   if(p2b.and.n.eq.2) return ! If we do p2b we get the Born and
                             ! virtuals from the structure functions

   if (n.eq.0) then
      iev=iev+1
      call pwhgaccumup
      recompute = .true. ! Signals that next time we have a new event cycle
      return
   endif

   ! It looks like, in a given set of calls (ie born + real + ...) eta
   ! can change, but x,y,Q2 stay the same. Eta however only changes a
   ! few times, so it is worth recomputing eta (which is cheap) and
   ! saving the pdfs (since they are more expensive to recompute).
   ETA=2*DOT(P,1,6)/S
   if(eta.ne.etasave) then
      etasave = eta
      call evolvePDF(eta,xmuf*scales_muf(iscales)*Qval,pdfs)
   endif
   
   if(recompute) then
      ! get x and Q2
      Q2=ABS(DOT(P,5,5))
      Qval=sqrt(Q2)
      X=ETA*Q2/(2*DOT(P,1,5))
      Y=DOT(P,1,5)/DOT(P,1,6)
      
      as2pi = alphasPDF(xmur*scales_mur(iscales)*Qval)/pi * 0.5d0 ! AK change this to alphasLocal at some point
      call p2bmomenta(x,y,Q2,p2bbreit,p2blab)
      Qlab(:) = p2blab(:,1) - p2blab(:,3)
      recompute = .false.
   endif

   if(order_max.eq.3.and.NA.eq.1) then ! Doing NLO in DISENT and this is the LO term. Include scale compensation.
      totwgt = dot_product(weight,pdfs)*(as2pi + two * as2pi**2 * b0 * log(xmur))
   else
      totwgt = dot_product(weight,pdfs)*as2pi**na
   endif
   dsig = totwgt * ncall2

   ! First we transfer the DISENT momenta to p2bee convention
   pbornbreit = 0 ! pborn(0:3,2+2)
   prealbreit = 0 ! preal(0:3,3+2)
   prrealbreit = 0 ! prreal(0:3,4+2)

   if(n.eq.2) then ! Born kinematics
      pbornbreit(:,1) = cshift(p(:,6),-1) ! Incoming lepton
      pbornbreit(:,2) = cshift(p(:,1),-1) ! Incoming quark
      pbornbreit(:,3) = cshift(p(:,7),-1) ! Outgoing lepton
      pbornbreit(:,4) = cshift(p(:,2),-1) ! Outgoing quark
      call mbreit2labdisent(n+2,Qlab,pbornbreit,pbornlab)
   elseif(n.eq.3) then
      prealbreit(:,1) = cshift(p(:,6),-1) ! Incoming lepton
      prealbreit(:,2) = cshift(p(:,1),-1) ! Incoming quark
      prealbreit(:,3) = cshift(p(:,7),-1) ! Outgoing lepton
      prealbreit(:,4) = cshift(p(:,2),-1) ! Outgoing quark
      prealbreit(:,5) = cshift(p(:,3),-1) ! Outgoing emission
      call mbreit2labdisent(n+2,Qlab,prealbreit,preallab)
   elseif(n.eq.4) then
      prrealbreit(:,1) = cshift(p(:,6),-1) ! Incoming lepton
      prrealbreit(:,2) = cshift(p(:,1),-1) ! Incoming quark
      prrealbreit(:,3) = cshift(p(:,7),-1) ! Outgoing lepton
      prrealbreit(:,4) = cshift(p(:,2),-1) ! Outgoing quark
      prrealbreit(:,5) = cshift(p(:,3),-1) ! Outgoing emission
      prrealbreit(:,6) = cshift(p(:,4),-1) ! Outgoing emission
      call mbreit2labdisent(n+2,Qlab,prrealbreit,prreallab)
   else
      print*, 'n = ', n
      stop 'Wrong n in user routine of DISENT'
   endif

   call user_analysis(n+2, dsig, x, y, Q2)

   if(p2b) then
      pbornbreit = p2bbreit
      pbornlab   = p2blab
      call user_analysis(2+2, -dsig, x, y, Q2)
   endif

   weights(NA) = weights(NA) + dsig ! Total weights
 END subroutine user

 subroutine finalise_histograms(prog)
   implicit none
   character (len=*) :: prog
   character(len=30) :: histname

   call pwhgsetout
   !call pwhgaddout
   if (order_max.eq.1) then
      histname=trim(prog)//'-hist'//'-LO-'//trim(seedstr)
   else if (order_max.eq.2) then
      histname=trim(prog)//'-hist'//'-NLO-'//trim(seedstr)
   else if (order_max.eq.3) then
      histname=trim(prog)//'-hist'//'-NNLO-'//trim(seedstr)
   else if (order_max.eq.4) then
      histname=trim(prog)//'-hist'//'-N3LO-'//trim(seedstr)
   endif
   call pwhgtopout(histname)

 end subroutine finalise_histograms

 ! Taken directly from DISENT
 FUNCTION DOT(P,I,J)
   IMPLICIT NONE
   !---RETURN THE DOT PRODUCT OF P(*,I) AND P(*,J)
   INTEGER I,J
   DOUBLE PRECISION DOT,P(4,7)
   DOT=P(4,I)*P(4,J)-P(3,I)*P(3,J)-P(2,I)*P(2,J)-P(1,I)*P(1,J)
 END FUNCTION DOT


end program
