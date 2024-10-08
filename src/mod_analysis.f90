module mod_analysis
  use mod_parameters
  implicit none
  
  double precision, public, save :: plab(0:3,6), pbreit(0:3,6)
  public analysis, init_histo
  
contains
  
  subroutine init_histo
    implicit none
    
    call inihists

    if(scaleuncert) call setupmulti(maxscales)

    call define_histograms
    
  end subroutine init_histo
  
  subroutine analysis(n,dsigma,x,y,Qsq)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: dsigma(maxscales), x, y, Qsq
    integer i
    
    if(all(dsigma.eq.0d0)) return
    
    ! Now fill the momenta
    if(n.eq.4) then ! Born kinematics
       pbreit(:,1:n) = pbornbreit(:,1:n)
       plab(:,1:n) = pbornlab(:,1:n)
    elseif(n.eq.5) then
       pbreit(:,1:n) = prealbreit(:,1:n)
       plab(:,1:n) = preallab(:,1:n)
    elseif(n.eq.6) then
       pbreit(:,1:n) = prrealbreit(:,1:n)
       plab(:,1:n) = prreallab(:,1:n)
    else
       stop 'Wrong n in analysis'
    endif
    
    if(any(dsigma.ne.dsigma)) then
       print*, 'Encountered NaN:', dsigma
       print*, 'n, x, y, Qsq',n,x,y,Qsq
       print*, 'Breit frame momenta'
       print*, 'incoming lepton: ', pbreit(:,1)
       print*, 'incoming parton: ', pbreit(:,2)
       print*, 'outcoming lepton:', pbreit(:,3)
       print*, 'outcoming parton:', pbreit(:,4)
       do i=5,n
          print*, 'emitted parton:  ', pbreit(:,i)
       enddo
       return ! Protection against NaN which can happen when running
              ! with a very small cut off in disent
    endif

    call user_analysis(n,dsigma,x,y,Qsq)
  end subroutine analysis
  
end module mod_analysis
