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
      real(dp) :: dsigma_real, pdf(-6:6), csi, z
      real(dp) :: xrand(4), vegas_weight
      integer vegas_ncall
      common/vegas_ncall/vegas_ncall

      dsigma = 0d0
      dsigma_real = 0d0
!     generate phase space using phase_space module
      call gen_phsp_born(xrand(1:2),x,y,Qsq,Qvec,jac(1),pbornlab,pbornbreit)
      jacborn = jac(1)
      
!     skip phase space points with vanishing jacobian or with Q < Qmin
      if ((jac(1).ne.0d0).and.(Qsq.gt.(Qmin**2))) then
!     compute dsigma using the squared hadronic tensor
         dsigma = eval_matrix_element(order_min,order_max, x, y, Qsq)
!     convert to [pb] and add in jacobian 
         dsigma = dsigma * gev2pb * jacborn
      endif
      
!     remove the rare outliers where we get dsigma = NaN
      if (dsigma.ne.dsigma) dsigma = 0d0
      if(fillplots) call user_analysis(4,dsigma*vegas_ncall*vegas_weight,x,y,Qsq)
      !if(dsigma.gt.9d3) print*, dsigma
      if(order_max.gt.1.and.p2b.and..false.) then ! Do one emission
            call gen_phsp_real(xrand(3:4),x,y,Qsq,Qvec,csi,z,jac(2),preallab,prealbreit)
            jacreal = jac(1) * jac(2)

            dsigma_real = real_emission_matrix_element(x,y,Qsq,csi,z)
            !dsigma_real = dsigma_real * gev2nb * jacreal
            dsigma_real = dsigma_real * gev2pb * jacreal

            if(fillplots) then
                  call user_analysis(5,dsigma_real*vegas_ncall*vegas_weight,x,y,Qsq)
                  ! Projection-to-born term 
                  call user_analysis(4,-dsigma_real*vegas_ncall*vegas_weight,x,y,Qsq)
            endif
      end if

      if(fillplots) call pwhgaccumup
      end function dsigma


      end module mod_dsigma
