!----------------------------------------------------------------------
! Set up HOPPET: the PDF tables (toy, LHAPDF or evolved from LHAPDF)
! and the structure functions, as needed before evaluating the matrix
! element.
module mod_pdf_setup
  use hoppet, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use mod_parameters
  use types
  implicit none

  private
  public :: read_PDF, setup_structure_functions

contains

  !----------------------------------------------------------------------
  ! (Re)start HOPPET with a grid reaching the largest scale needed,
  ! fill the PDF table and initialise the structure functions.
  subroutine setup_structure_functions()
   real(dp) :: rts

   rts = sqrt(s)
   if(scaleuncert) rts = rts * maxval(scales_muf) * xmuf
   ! Need to start hoppet
   maxQval = rts 
   call hoppetStartExtended(ymax_hoppet,dy,minQval,maxQval,dlnlnQ,nloop,&
        &         order_hoppet,factscheme_MSbar)
   if(vnf) then
      call StartStrFct(order_max = order_max, scale_choice =&
           & scale_choice_hoppet, constant_mu = mz, param_coefs =&
           & .true. , wmass = mw, zmass = mz)
   else
      call StartStrFct(order_max = order_max, nflav = nflav,&
           & scale_choice = scale_choice_hoppet, constant_mu = mz,&
           & param_coefs = .true., wmass = mw, zmass = mz)      
   endif
   call read_PDF()
   call InitStrFct(order_max, separate_orders = separate_orders, xR =&
        & xmur, xF = xmuf)
  end subroutine setup_structure_functions

  !----------------------------------------------------------------------
  ! fill the streamlined interface PDF table for the structure
  ! functions only.
  subroutine read_PDF()
    use toy_pdfs
    use streamlined_interface
    real(dp), external :: alphasPDF
    interface
       subroutine EvolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine EvolvePDF
    end interface
    real(dp) :: res_lhapdf(-6:6), x, Q
    real(dp) :: res_hoppet(-6:6)
    real(dp) :: toy_pdf_at_Q0(0:grid%ny,ncompmin:ncompmax)
    real(dp) :: pdf_at_Q0(0:grid%ny,ncompmin:ncompmax)
    
    if (toy_Q0 > zero) then
       write(6,*) "WARNING: Using toy PDF"
       if(pdfname.eq.'toyHERALHC') then
          toy_pdf_at_Q0 = unpolarized_dummy_pdf(xValues(grid))
       elseif(pdfname.eq.'toyNF5') then
          toy_pdf_at_Q0 = unpolarized_toy_nf5_pdf(xValues(grid))
       else
          print*, 'Did not recognise toy PDF:', pdfname
          call exit()
       end if
       if(vnf) then
          call InitRunningCoupling(coupling, toy_alphas_Q0,&
               & toy_Q0, nloop, -1000000045, masses(4:6)&
               &, .true.)
       else
          call InitRunningCoupling(coupling, toy_alphas_Q0,&
               & toy_Q0, nloop, nflav, masses(4:6)&
               &, .true.)
       endif
       call EvolvePdfTable(tables(0), toy_Q0, toy_pdf_at_Q0, dh,&
            & coupling, nloop=nloop)
       setup_done(0)  = .true. ! This signals to HOPPET that we have set up the PDFs (since we don't use the streamlined interface)
    elseif (Q0pdf > zero) then
       write(6,*) "WARNING: Using internal HOPPET DGLAP evolution"
       call InitPDF_LHAPDF(grid, pdf_at_Q0, EvolvePDF, Q0pdf)
       
       if(vnf) then
          call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , order_max,&
               & -1000000045, masses(4:6), .true.)
       else
          call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , order_max,&
               & nflav, masses(4:6), .true.)
       end if
       !call EvolvePdfTable(tables(0), Q0pdf, pdf_at_Q0, dh, coupling, &
       !     &  muR_Q=xmuR_PDF, nloop=min(order_max,3))
       call EvolvePdfTable(tables(0), Q0pdf, pdf_at_Q0, dh, coupling, &
            &  muR_Q=xmuR_PDF, nloop=order_max)
       setup_done(0)  = .true. ! This signals to HOPPET that we have set up the PDFs (since we don't use the streamlined interface)
    else
!       if(vnf) then
!          call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , order_max,&
!               & -1000000045, masses(4:6), .true.)
!       else
!          call InitRunningCoupling(coupling, alphasPDF(MZ) , MZ , order_max,&
!               & nflav, masses(4:6), .true.)
       !       end if
       call hoppetSetCoupling(alphasPDF(MZ), MZ, order_max)
       call hoppetAssign(EvolvePDF)
    endif

 end subroutine read_PDF
end module mod_pdf_setup
