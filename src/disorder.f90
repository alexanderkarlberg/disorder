program disorder
  use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use sub_defs_io
  !use dummy_pdfs
  use mod_matrix_element
  use mod_parameters
  use mod_phase_space
  use mod_analysis
  use integration
  use types
  use mod_dsigma
  use mod_pdf_setup
  use mod_disent_interface
  implicit none
  integer, parameter :: ndim = 7 ! dimension for vegas integration
  integer, parameter :: nmempdfMAX = 200 ! max number of pdfs
  integer  :: nmempdf_start, nmempdf_end, imempdf, ncall2_save
  real(dp) :: integ, error_int, proba, tini, tfin
  real(dp) :: sigma_tot, error_tot, region(1:2*ndim),xdum(4)
  real(dp) :: res(0:nmempdfMAX), central, errminus, errplus, errsymm,&
       & respdf(0:nmempdfMAX),resas(0:nmempdfMAX), central_dummy
  real(dp) :: res_scales(1:maxscales),maxscale,minscale
  real(dp) :: NC_reduced_central, NC_reduced_max, NC_reduced_min
  real(dp) :: CC_reduced_central, CC_reduced_max, CC_reduced_min
  character * 100 :: analysis_name, histo_name, scale_string
  character * 3 :: pdf_string
  logical :: scaleuncert_save

  integer vegas_ncall
  common/vegas_ncall/vegas_ncall

  call cpu_time(tini)

  ! set up all constants and parameters from command line arguments
  call set_parameters()
  ! Need to start hoppet
  if(vnf) then
    call hoppetSetPoleMassVFN(mc, mb, mt)
  else
    call hoppetSetFFN(nflav)
  endif
  call hoppetStartExtended(ymax_hoppet,dy,minQval,maxQval,dlnlnQ,nloop,&
       &         order_hoppet,factscheme_MSbar)
  call read_PDF()
  ! Initialise histograms
  call init_histo()

  call welcome_message

  call print_header(0)

  if(inclusive) then
     if (pdfuncert) then
        nmempdf_start = 0
        call numberPDF(nmempdf_end)
     else
        nmempdf_start = nmempdf
        nmempdf_end   = nmempdf
     endif
     
     if (nmempdf_end .gt. nmempdfMAX) stop "ERROR: increase nmempdfMAX"
     
     do imempdf = nmempdf_start,nmempdf_end
        call initialise_run_structure_functions
        if(do_analysis) then
           call pwhgsetout
           call finalise_histograms
           call setupmulti(1) ! Tell analysis only one scale from now on
        endif
        Nscales = 1 ! To avoid computing scale variations in the next calls
        scale_choice = 1 ! Sets up the fast tables for a central scale choice
        ncall1 = 0 ! turn off the grid computation
     enddo ! end loop over pdfs
  else
     imempdf = nmempdf
     nmempdf_start = nmempdf
     
     ! Since ncall2 sets both the calls to disent and the structure
     ! functions, and the structure functions do not need as much
     ! stats as disent, we set it locally here to a smaller value and
     ! reset before calling disent.
     ncall2_save = ncall2
     ncall2 = max(100000,ncall2)/10

     call initialise_run_structure_functions
     call pwhgaddout

     ncall2 = ncall2_save

     if(order_max.gt.1) then
        ! Then do the disent run
        call DISENTFULL(ncall2,S,nflav,user,dis_cuts,12345&
             &+iseed-1,67890+iseed-1,NPOW1,NPOW2,CUTOFF ,order_max-1&
             &,disent_muf,cflcl,calcl,trlcl,scaleuncert)
        ! Store disent result
        call pwhgaddout
     endif
     
     call finalise_histograms
  endif ! inclusive

  analysis_name='xsct'
  analysis_name=trim(prefix)//"xsct_"//trim(adjustl(order))//"_seed"//seedstr//".dat"
  
  call print_results(0, '')
  call print_results(11, analysis_name)

  
  call cpu_time(tfin)
  write(6,'(a)')
  write(6,'(a,es9.2,a)') '==================== TOTAL TIME : ', tfin&
       &-tini, ' s.'
  write(6,'(a)')

contains

 subroutine initialise_run_structure_functions
   implicit none
   
   if(toy_Q0 < 0d0) then
      write(6,*) "PDF member:",imempdf
      call InitPDF(imempdf)
      call getQ2min(0,Qmin)
      Qmin = sqrt(Qmin)
   endif

   if(Qmin.gt.sqrt(Q2min)) then
      print*, 'WARNING: PDF Qmin = ', Qmin
      print*, 'But running with input value of: ', sqrt(Q2min)
      stop
   endif

   scaleuncert_save = scaleuncert
   
   call setup_structure_functions()

   if(novegas) then ! This means Q and x fixed
      ! Need dummy random numbers
      xdum = 0.5_dp
      fillplots = .true.
      vegas_ncall = 1
      sigma_tot = dsigma(xdum, one)
      error_tot = zero
      res(imempdf) = sigma_tot
      res_scales(1:nscales) = sigma_all_scales(1:nscales) 
   else
      region(1:ndim)        = zero
      region(ndim+1:2*ndim) = one
      sigma_tot = zero
      error_tot = zero

      ! vegas warmup call
      if(ncall1.gt.0.and.itmx1.gt.0) then 
         ! Skip grid generation if grid is being read from file
         writeout=.true.
         fillplots = .false.
         scaleuncert = .false.
         call vegas(region,ndim,dsigma,0,ncall1,itmx1,0,integ,error_int,proba)
         sigma_all_scales = zero ! Reset for the production run below
         NC_reduced_sigma = zero
         CC_reduced_sigma = zero
         scaleuncert = scaleuncert_save
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
      fillplots = do_analysis !.true.
      call vegas(region,ndim,dsigma,1,ncall2,itmx2,0,integ,error_int,proba)
      readin = .true.
      ! add integral to the total cross section
      sigma_tot = sigma_tot + integ
      error_tot = error_tot + error_int**2

      res(imempdf) = integ
      res_scales(1:nscales) = sigma_all_scales(1:nscales)
   endif
   if(imempdf.eq.nmempdf_start) then ! First PDF, this is where we compute scale uncertainties
      maxscale = maxval(res_scales(1:Nscales)) 
      minscale = minval(res_scales(1:Nscales))
      res(imempdf) = res_scales(1) ! Copy central scale
      NC_reduced_central = NC_reduced_sigma(1)
      CC_reduced_central = CC_reduced_sigma(1)
      NC_reduced_max = maxval(NC_reduced_sigma(1:Nscales))
      CC_reduced_max = maxval(CC_reduced_sigma(1:Nscales))
      NC_reduced_min = minval(NC_reduced_sigma(1:Nscales))
      CC_reduced_min = minval(CC_reduced_sigma(1:Nscales))
   endif
 end subroutine initialise_run_structure_functions

 subroutine finalise_histograms
   implicit none

   ! print total cross section and error into file  
   ! construct name of output file
   write(pdf_string,"(I3.3)") imempdf
   scale_string = "_pdfmem"//trim(pdf_string)
   histo_name="disorder_"//trim(adjustl(order))//"_seed"//seedstr//scale_string//".dat"
   
   call pwhgtopout(histo_name)
   call resethists

 end subroutine finalise_histograms

 subroutine print_results(idev,filename)
   implicit none
   integer, intent(in) :: idev
   character * 100, intent(in) :: filename
   
   if(idev.gt.0) then ! This is already printed to screen earlier
      OPEN(UNIT=idev, FILE=filename, ACTION="write")
      call print_header(idev)
   endif
   write(idev,*) ''
   write(idev,*) '============================================================'
   if (order_max.eq.1) then 
      write(idev,'(a,es13.6,a,es13.6,a)') ' # Total LO cross-section'
   else if (order_max.eq.2) then 
      write(idev,'(a,es13.6,a,es13.6,a)') ' # Total NLO cross-section'
   else if (order_max.eq.3) then 
      write(idev,'(a,es13.6,a,es13.6,a)') ' # Total NNLO cross-section'
   else if (order_max.eq.4) then 
      write(idev,'(a,es13.6,a,es13.6,a)') ' # Total N3LO cross-section'
   endif

   central = res(0)

   write(idev,'(a)') ' # Summary:'
   if(NC.and.CC) then
      if(Q2min.eq.Q2max) then
         write(idev,'(a,f16.6,a)') ' # σ(NC + CC)                     =', central,' pb/GeV^2'
      else
         write(idev,'(a,f16.6,a)') ' # σ(NC + CC)                     =', central,' pb'
      endif
   elseif(NC) then
      if(Q2min.eq.Q2max) then
         write(idev,'(a,f16.6,a)') ' # σ(NC)                          =', central,' pb/GeV^2'
      else
         write(idev,'(a,f16.6,a)') ' # σ(NC)                          =', central,' pb'
      endif
   elseif(CC) then
      if(Q2min.eq.Q2max) then
         write(idev,'(a,f16.6,a)') ' # σ(CC)                          =', central,' pb/GeV^2'
      else
         write(idev,'(a,f16.6,a)') ' # σ(CC)                          =', central,' pb'
      endif
   endif
   write(idev,'(a,f14.4,a)') ' # MC integration uncertainty     =', sqrt(error_tot)/central*100.0_dp, ' %'

   if(scaleuncert) then
      write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (+)      =',&
           & ((maxscale-central)/central)*100.0_dp, ' %'
      write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (-)      =',&
           & ((minscale-central)/central)*100.0_dp, ' %'
   endif

   if(pdfuncert) then
      call getpdfuncertainty(res(nmempdf_start:nmempdf_end),central,errplus,errminus,errsymm)
      write(idev,'(a,f14.4,a)') ' # PDF symmetric uncertainty*     =', errsymm/central*100.0_dp, ' %'
      if(alphasuncert) then
        central_dummy = central
        respdf = central
        resas = central
        respdf(0:nmempdf_end-2) = res(0:nmempdf_end-2)
        resas(nmempdf_end-1:nmempdf_end) = res(nmempdf_end-1:nmempdf_end)
        call getpdfuncertainty(respdf(nmempdf_start:nmempdf_end),central_dummy,errplus,errminus,errsymm)
        write(idev,'(a,f14.4,a)') ' # Pure PDF symmetric uncertainty =', errsymm/central*100.0_dp, ' %'
        call getpdfuncertainty(resas(nmempdf_start:nmempdf_end),central_dummy,errplus,errminus,errsymm)
        write(idev,'(a,f14.4,a)') ' # Pure αS symmetric uncertainty  =', errsymm/central*100.0_dp, ' %'
     endif
      write(idev,'(a)') ' # (*PDF uncertainty contains alphas uncertainty if using a  '
      write(idev,'(a)') ' #   PDF set that supports it (eg PDF4LHC15_nnlo_100_pdfas)).'
   endif
   write(idev,*) ''
   
   if (order_max.eq.1) then 
      write(idev,'(a,es13.6,a,es13.6,a)') ' # Reduced LO cross-sections'
   else if (order_max.eq.2) then 
      write(idev,'(a,es13.6,a,es13.6,a)') ' # Reduced NLO cross-sections'
   else if (order_max.eq.3) then 
      write(idev,'(a,es13.6,a,es13.6,a)') ' # Reduced NNLO cross-sections'
   else if (order_max.eq.4) then 
      write(idev,'(a,es13.6,a,es13.6,a)') ' # Reduced N3LO cross-sections'
   endif

   if(NC.and.CC) then
      write(idev,'(a,f16.6)') ' # σ reduced (NC)                 =', NC_reduced_central
      if(scaleuncert) then
         write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (+)      =',&
              & ((NC_reduced_max-NC_reduced_central)/NC_reduced_central)&
              &*100.0_dp, ' %'
         write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (-)      =',&
              & ((NC_reduced_min-NC_reduced_central)/NC_reduced_central)&
              &*100.0_dp, ' %'
      endif
      write(idev,'(a,f16.6)') ' # σ reduced (CC)                 =', CC_reduced_central
      if(scaleuncert) then
         write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (+)      =',&
              & ((CC_reduced_max-CC_reduced_central)/CC_reduced_central)&
              &*100.0_dp, ' %'
         write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (-)      =',&
              & ((CC_reduced_min-CC_reduced_central)/CC_reduced_central)&
              &*100.0_dp, ' %'
      endif
   elseif(NC) then
      write(idev,'(a,f16.6)') ' # σ reduced (NC)                 =', NC_reduced_central
      if(scaleuncert) then
         write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (+)      =',&
              & ((NC_reduced_max-NC_reduced_central)/NC_reduced_central)&
              &*100.0_dp, ' %'
         write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (-)      =',&
              & ((NC_reduced_min-NC_reduced_central)/NC_reduced_central)&
              &*100.0_dp, ' %'
      endif
   elseif(CC) then
      write(idev,'(a,f16.6)') ' # σ reduced (CC)                 =', CC_reduced_central
      if(scaleuncert) then
         write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (+)      =',&
              & ((CC_reduced_max-CC_reduced_central)/CC_reduced_central)&
              &*100.0_dp, ' %'
         write(idev,'(a,f14.4,a)') ' # QCD scale uncertainty (-)      =',&
              & ((CC_reduced_min-CC_reduced_central)/CC_reduced_central)&
              &*100.0_dp, ' %'
      endif
   endif
      write(idev,*) '============================================================'


 end subroutine print_results
end program
