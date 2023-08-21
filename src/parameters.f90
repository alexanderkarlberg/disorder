!----------------------------------------------------------------------
! A module to define all parameters of the run, which can be read from
! command line arguments
module parameters
  use sub_defs_io
  use integration
  use types
  implicit none

  private
  real(dp), parameter, public :: alpha_em = 1.0_dp/137.0_dp
  real(dp), parameter, public :: gev2pb = 389379660.0_dp
  real(dp), parameter, public :: gev2nb = 389379.66_dp
  real(dp), parameter, public :: eps    = 1.0e-14_dp
  integer, parameter, public :: maxscales = 9
  real(dp), parameter, public :: scales_mur(1:maxscales) = &
       & (/1.0_dp, 2.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 2.0_dp, 0.5_dp, 0.25_dp, 4.0_dp/)
  real(dp), parameter, public :: scales_muf(1:maxscales) = &
       & (/1.0_dp, 2.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 0.25_dp, 4.0_dp/)
  real(dp), public :: xmuf, xmur, Qmin
  real(dp), public :: sin_thw, mw, mz, w_width, z_width
  real(dp), public :: sqrts, S, Q0_cut_sq
  integer,  public :: order_min, order_max, iseed, scale_choice
  integer,  public :: nflav, ipdf, it1, itmx1, itmx2, ncall1, ncall2, nscales
  character * 4, public :: seedstr
  character * 17, public :: scalestr(maxscales)
  character(len=50), public :: pdfname
  integer, public :: nmempdf, outdev
  logical, public :: pdfuncert, scaleuncert3, scaleuncert7, scaleuncert9, fillplots&
       &, p2b, gluon_only, noZ, positron, Zonly, intonly, scaleuncert, new_scaleuncert
  real(dp), public :: Q2min, Q2max, xmin, xmax, ymin, ymax,ymn,ymx,&
       & Eh, El, sigma_all_scales(maxscales)
  character (len=4), private :: order
  real(dp), private :: Q, x, y
  real(dp), public :: pbornlab(0:3,2+2), preallab(0:3,2+3), prreallab(0:3,2+4)
  real(dp), public :: pbornbreit(0:3,2+2), prealbreit(0:3,2+3), prrealbreit(0:3,2+4)

  real(dp), public, save :: CFlcl, CAlcl, Trlcl, b0

  real(dp), public, parameter :: Qup =  2.0_dp/3.0_dp
  real(dp), public, parameter :: Qdn = -1.0_dp/3.0_dp
  real(dp), public :: sin_thw_sq, sin_2thw_sq,Ve, Ae, Ve2, Ae2, Ve2_Ae2, two_Ve_Ae ! Vector and axial couplings of the electron

  logical, public :: run_disent

    ! VEGAS common blocks
  integer, public :: ilast
  common/last_integ/ilast
  integer, public :: saveseed,idum
  COMMON /ranno/idum
  
  public :: set_parameters

contains

  ! set all parameters from input card or command line arguments
  subroutine set_parameters ()
    integer :: call, i
    real(dp) :: ran
    mw           = dble_val_opt("-mw",80.398_dp)
    mz           = dble_val_opt("-mz",91.1876_dp)
    nflav        = int_val_opt ("-nf",5)
    w_width      = dble_val_opt("-wwidth",2.141_dp)
    z_width      = dble_val_opt("-zwidth",2.4952_dp)
    CAlcl        = dble_val_opt("-CA",3.0_dp)
    CFlcl        = dble_val_opt("-CF",4.0_dp/3.0_dp)
    TRlcl        = dble_val_opt("-Tr",0.5_dp)
    
    ! Compute β0 as needed for scale compensation
    b0 = (11.0_dp * CAlcl - 4.0_dp * nflav * TRlcl) / 6.0_dp
    
    order_min    = int_val_opt ('-order-min',1)
    order_max    = int_val_opt ('-order-max',3)
    order = 'NNLO'
    sin_thw_sq = 1.0_dp - (mw/mz)**2
    sin_2thw_sq = 4.0_dp * (1.0_dp - sin_thw_sq) * sin_thw_sq
    ! if "-lo/-nlo/-nnlo/-n3lo" command is given, overwrite order_min and order_max accordingly
    if (log_val_opt("-lo")) then
       order_min = 1
       order_max = 1
       order = '  LO'
    elseif (log_val_opt("-nlo")) then
       order_min = 1
       order_max = 2
       order = ' NLO'
    elseif (log_val_opt("-nnlo")) then
       order_min = 1
       order_max = 3
       order = 'NNLO'
    elseif (log_val_opt("-n3lo")) then
       order_min = 1
       order_max = 4
       order = 'N3LO'
    endif
    
    noZ = .not.log_val_opt ("-includeZ",.false.)
    Zonly = .false.
    intonly = .false.
    if(.not.noZ) Zonly = log_val_opt ("-Zonly",.false.)
    if(.not.noZ) intonly = log_val_opt ("-intonly",.false.)
    if(Zonly.and.intonly) stop 'Cannot run with both Z and interference ONLY flags'
    positron = log_val_opt ("-positron",.false.)
    run_disent = log_val_opt("-disent",.false.)

    ! For a positron the Axial coupling flips sign
    Ae = - 0.5_dp
    if(positron) Ae = - Ae
    Ve = - 0.5_dp + 2.0_dp * sin_thw_sq 
    Ae2 = Ae**2
    Ve2 = Ve**2
    Ve2_Ae2 = Ve2 + Ae2
    two_Ve_Ae = 2.0_dp * Ve * Ae

    scale_choice = int_val_opt ('-scale-choice',1)
    readin       = log_val_opt ("-readingrid")
    xmuf         = dble_val_opt("-xmuf",1.0_dp)
    xmur         = dble_val_opt("-xmur",1.0_dp)
    pdfname      = string_val_opt("-pdf", "NNPDF30_nlo_as_0118_hera")
    ! pdfname      = string_val_opt("-pdf", "PDF4LHC15_nnlo_mc")
    nmempdf      = int_val_opt ("-nmempdf",0)
    pdfuncert    = log_val_opt ("-pdfuncert")
    scaleuncert3 = log_val_opt ("-3scaleuncert")
    scaleuncert7 = log_val_opt ("-7scaleuncert")
    scaleuncert9 = log_val_opt ("-9scaleuncert")
    new_scaleuncert = log_val_opt ("-new-scaleuncert")
    if(new_scaleuncert) scale_choice = 2
    scaleuncert = scaleuncert3.or.scaleuncert7.or.scaleuncert9
    if(scaleuncert9) then
       Nscales = 9
    elseif(scaleuncert3) then
       Nscales = 3
    elseif(scaleuncert7) then
       Nscales = 7
    else
       Nscales = 1
    endif
     
    p2b = .false.
    p2b = log_val_opt("-p2b")
    if(.not.noZ.and.p2b) stop 'Cannot do Z in p2b yet'
    if(p2b.and.pdfuncert) stop "Cannot do pdf uncertainties with p2b"
!    if(p2b.and.scaleuncert3) stop "Cannot do scale uncertainties with p2b"
!    if(p2b.and.scaleuncert7) stop "Cannot do scale uncertainties with p2b"
    gluon_only = .false.
    gluon_only = log_val_opt("-gluon_only")
    if(scaleuncert3.and.scaleuncert7) then
       write(*,*) 'WARNING: Have to do either 3 or 7 scales. Cannot do both. Doing none.'
       scaleuncert3 = .false.
       scaleuncert7 = .false.
    endif
    ncall1       = int_val_opt ("-ncall1",100000)
    ncall2       = int_val_opt ("-ncall2",100000)
    it1         = int_val_opt("-it1",1)
    itmx1        = int_val_opt ("-itmx1",3)
    itmx2        = int_val_opt ("-itmx2",1)
    iseed        = int_val_opt ("-iseed",10)

    !     Read in bounds on x,y,Q2
    Q = dble_val_opt("-Q",-1.0_dp)
    if(Q.lt.0.0_dp) then
      Q2min = (dble_val_opt("-Qmin",1.0_dp))**2
      Q2max = (dble_val_opt("-Qmax",1d10))**2
    else 
      Q2min = Q**2
      Q2max = Q2min
    endif
    x = dble_val_opt("-x",-1.0_dp)
    if(x.lt.0.0_dp) then
      xmin = dble_val_opt("-xmin",0.0_dp)
      xmax = dble_val_opt("-xmax",1.0_dp)
    else
      xmin = x
      xmax = xmin
    endif
    y = dble_val_opt("-y",-1.0_dp)
    if(y.lt.0.0_dp) then
      ymin = dble_val_opt("-ymin",0.0_dp)
      ymax = dble_val_opt("-ymax",1.0_dp)
    else
      ymin = y
      ymax = ymin
    endif

    !     Set centre of mass energy
    El = dble_val_opt("-Elep",27.5_dp)
    Eh = dble_val_opt("-Epro",820.0_dp)
    s = 4d0 * Eh * El

    !     Sanity checks
    if(Q2min.lt.0d0) stop 'Q2min negative'
    if(xmin.lt.0d0) stop 'xmin negative'
    if(ymin.lt.0d0) stop 'ymin negative'
    if(Q2min.gt.Q2max) stop 'Q2min > Q2max'
    if(xmin.gt.xmax) stop 'xmin > xmax'
    if(ymin.gt.ymax) stop 'ymin > ymax'
    if(Q2max.gt.s) Q2max = s

    if(Q2min.eq.Q2max.and.xmin.eq.xmax.and.ymin.eq.ymax) stop 'Only two variables can be constrained'

    if(s*xmax*ymax.lt.Q2min .or. s*xmin*ymin.gt.Q2max) stop 'No phase space available'

    if(Q2min.eq.Q2max.and.xmin.eq.xmax) then
       ymin = Q2min/(s*xmin)
       ymax = ymin
    else
       if(xmin*(ymax*s) .lt. Q2min) xmin = Q2min/(ymax*s)
       if(xmax*(ymin*s) .gt. Q2max) xmax = Q2max/(ymin*s)
    endif
    !         if(s*xmin*ymin .gt. Q2min) Q2min = s*xmin*ymin
    !         if(s*xmax*ymax .gt. Q2max) Q2max = s*xmax*ymax
    sqrts = sqrt(s * xmax)  
    print*, '****************************************'
    print*, 'Doing DIS @ ', trim(order)
    if(run_disent.and..not.p2b) then
       print*, 'Using DISENT'
    elseif(run_disent.and.p2b) then
       print*, 'Using DISENT with projection-to-Born'
    endif
    if(noZ)  print*, 'With γ only'
    if(.not.noZ.and.(.not.Zonly).and.(.not.intonly))  print*, 'With γ/Z'
    if(Zonly)  print*, 'With Z only'
    if(intonly)  print*, 'With γ/Z interference only'
    print*, 'xmin, xmax', xmin, xmax
    print*, 'ymin, ymax', ymin, ymax
    print*, 'Q2min, Q2max', Q2min, Q2max, 'GeV^2'
    print*, 'Electron energy: ', El, 'GeV'
    print*, 'Proton energy: ', Eh, 'GeV'
    print*, 'COM energy: ', S, 'GeV^2'
    print*, '****************************************'    

    ! various setup
    write(seedstr,"(I4.4)") iseed

    idum = -iseed ! Initial seed
    if(readin.or.it1.gt.1) then
       readin = .true.
       ! Fastforward random seed
       do i = 1,it1
          do call = 1,ncall1
             ran = ran2(idum)
          enddo
       enddo
       idum = -idum
    endif
    ingridfile ='grids-'//seedstr//'.dat'
    outgridfile='grids-'//seedstr//'.dat'
    outgridtopfile='grids-'//seedstr//'.top'
    ilast=0
    scalestr(1) = '_μR_1.0_μF_1.0_'
    scalestr(2) = '_μR_2.0_μF_2.0_'
    scalestr(3) = '_μR_0.5_μF_0.5_'
    scalestr(4) = '_μR_1.0_μF_2.0_'
    scalestr(5) = '_μR_1.0_μF_0.5_'
    scalestr(6) = '_μR_2.0_μF_1.0_'
    scalestr(7) = '_μR_0.5_μF_1.0_'
    scalestr(8) = '_μR_.25_μF_.25_'
    scalestr(9) = '_μR_4.0_μF_4.0_'

    Qmin = 1.0_dp

    sigma_all_scales = 0.0_dp

    !outdev = idev_open_opt("-out")

    if (.not.CheckAllArgsUsed(0)) call exit()
  end subroutine set_parameters

end module parameters
