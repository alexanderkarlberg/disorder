      module mod_dsigma
      use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
      use parameters
      use types

      implicit none
     
      contains

!------------------------------------------------------------
!     dsigma function
      double precision function dsigma(xrand, vegas_weight)
      use matrix_element
      use phase_space
      implicit none
!     xrand contains a vector of random numbers in [0,1]
      real(dp) :: x, y, Qsq, Qval,jac(1:2), jacborn, jacreal, Qvec(0:3)
      real(dp) :: dsigma_real, pdf(-6:6), csi, z, dsigma_all_scales(maxscales)
      real(dp) :: xrand(4), vegas_weight
      integer vegas_ncall
      common/vegas_ncall/vegas_ncall

      dsigma = 0d0
      dsigma_all_scales = 0d0
!     generate phase space 
      call gen_phsp_born(xrand(1:2),x,y,Qsq,Qvec,jac(1),pbornlab,pbornbreit)
      jacborn = jac(1)
!     skip phase space points with vanishing jacobian or with Q < Qmin
      if (.not.(jac(1).ne.0d0).and.(Qsq.gt.(Qmin**2))) return
      if(scaleuncert) then
         dsigma_all_scales = eval_matrix_element_scale_variation(order_min,order_max, x, y, Qsq, nscales)
         dsigma_all_scales = dsigma_all_scales * gev2pb * jacborn 
         dsigma = dsigma_all_scales(1)
         dsigma_all_scales = dsigma_all_scales * vegas_weight 
         sigma_all_scales = sigma_all_scales + dsigma_all_scales
      else
         dsigma = eval_matrix_element(order_min,order_max, x, y, Qsq)
         dsigma = dsigma * gev2pb * jacborn
         dsigma_all_scales = dsigma * vegas_weight
      endif
      
!     remove the rare outliers where we get dsigma = NaN
      if (dsigma.ne.dsigma) then
         dsigma = 0d0
         dsigma_all_scales = 0d0
      endif
      if(fillplots) then
         call user_analysis(4,dsigma_all_scales*vegas_ncall,x,y,Qsq)
         call pwhgaccumup
      endif
      end function dsigma


      end module mod_dsigma
