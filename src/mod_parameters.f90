!----------------------------------------------------------------------
! A module to define all parameters of the run, which can be read from
! command line arguments
module mod_parameters
  use sub_defs_io
  use integration
  use types
  use streamlined_interface, CouplingValue => Value
  use dummy_pdfs
  implicit none

  private

  public print_header, welcome_message, alphasLocal

  !  real(dp), parameter, public :: gev2pb = 389379660.0_dp
  real(dp), parameter, public :: gev2pb = 389379372.1_dp
  !  real(dp), parameter, public :: gev2nb = 389379.66_dp
  real(dp), parameter, public :: gev2nb = 389379.3721_dp
  real(dp), parameter, public :: eps    = 1.0e-14_dp
  integer, parameter, public :: maxscales = 7
  real(dp), parameter, public :: scales_mur(1:maxscales) = &
       & (/1.0_dp, 2.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 2.0_dp, 0.5_dp/)
  real(dp), parameter, public :: scales_muf(1:maxscales) = &
       & (/1.0_dp, 2.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, 1.0_dp, 1.0_dp/)
  real(dp), public :: xmuf, xmur, Qmin
  real(dp), public :: mw, mz, w_width, z_width, GF, alpha_em,&
       & sin_thw_sq, sin_2thw_sq
  real(dp), public :: sqrts, S, Q0_cut_sq
  integer,  public :: order_min, order_max, iseed, scale_choice, scale_choice_hoppet
  integer,  public :: nflav, ipdf, it1, itmx1, itmx2, ncall1, ncall2, nscales
  character * 4, public :: seedstr
  character * 17, public :: scalestr(maxscales)
  character(len=50), public :: pdfname, outname, prefix
  integer, public :: nmempdf, outdev
  logical, public, save :: pdfuncert, alphasuncert, fillplots, p2b, noZ, positron,&
       & Zonly, intonly, scaleuncert, inclusive, novegas, NC, CC, vnf, help, do_analysis, separate_orders
  real(dp), public, save :: Q2min, Q2max, xmin, xmax, ymin, ymax,ymn,ymx,&
       & Eh, El, sigma_all_scales(maxscales),&
       & NC_reduced_dsigma(maxscales), CC_reduced_dsigma(maxscales),&
       & NC_reduced_sigma(maxscales), CC_reduced_sigma(maxscales)
  real(dp), public, save :: toy_Q0, Q0pdf, xmuR_PDF, Q2minPDF ! For HOPPET PDF evolution
  real(dp), public :: dy, dlnlnQ, minQval, maxQval, ymax_hoppet
  integer, public :: nloop, order_hoppet
  character (len=4), private :: order
  real(dp), private :: Q, x, y
  real(dp), public :: xlmin, xlmax
  real(dp), public :: pbornlab(0:3,2+2), preallab(0:3,2+3), prreallab(0:3,2+4)
  real(dp), public :: pbornbreit(0:3,2+2), prealbreit(0:3,2+3), prrealbreit(0:3,2+4)

  real(dp), public, save :: CFlcl, CAlcl, Trlcl, b0, NPOW1, NPOW2, CUTOFF

  real(dp), public :: Ve, Ae, Ve2, Ae2, Ve2_Ae2, two_Ve_Ae ! Vector and axial couplings of the electron

  ! VEGAS common blocks
  integer, public :: ilast
  common/last_integ/ilast
  integer, public :: saveseed,idum
  COMMON /ranno/idum

  public :: set_parameters

contains

  ! set all parameters from input card or command line arguments
  subroutine set_parameters 
    integer :: call, i
    real(dp) :: ran

    help = log_val_opt("-help",.false.)
    if(help) then
       call help_message
       call exit()
    endif

    ! Some physical parameters and constants
    mw           = dble_val_opt("-mw",80.398_dp)
    mz           = dble_val_opt("-mz",91.1876_dp)
    nflav        = int_val_opt ("-nf",5)
    vnf          = log_val_opt("-vnf",.false.)
    w_width      = dble_val_opt("-wwidth",2.141_dp)
    z_width      = dble_val_opt("-zwidth",2.4952_dp)
    alpha_em     = 1.0_dp/dble_val_opt("-one-over-alpha",137.0_dp)
    CAlcl        = dble_val_opt("-CA",3.0_dp)
    CFlcl        = dble_val_opt("-CF",4.0_dp/3.0_dp)
    TRlcl        = dble_val_opt("-Tr",0.5_dp)
    ! Compute β0 as needed for scale compensation in disent
    b0 = (11.0_dp * CAlcl - 4.0_dp * nflav * TRlcl) / 6.0_dp
    sin_thw_sq = 1.0_dp - (mw/mz)**2
    sin_2thw_sq = 4.0_dp * (1.0_dp - sin_thw_sq) * sin_thw_sq
    !    GF           = 1.1663787D-5 
    !    GF           =  3.14159265359_dp * alpha_em / sqrt(2.0_dp) / mw**2/sin_thw_sq
    GF           =  4.0_dp * atan(1.0_dp) * alpha_em / sqrt(2.0_dp) / mw**2/sin_thw_sq
    positron = log_val_opt ("-positron",.false.)
    ! For a positron the Axial coupling flips sign
    Ae = - 0.5_dp
    if(positron) Ae = - Ae
    Ve = - 0.5_dp + 2.0_dp * sin_thw_sq 
    Ae2 = Ae**2
    Ve2 = Ve**2
    Ve2_Ae2 = Ve2 + Ae2
    two_Ve_Ae = 2.0_dp * Ve * Ae

    ! The order at which we are running is read here along with
    ! whether or not we are doing inclsuvie/p2b and NC/CC.
    order_min    = int_val_opt ('-order-min',1)
    order_max    = int_val_opt ('-order-max',3)
    order = 'NNLO'
    separate_orders = log_val_opt("-separate-orders",.false.)
    ! if "-lo/-nlo/-nnlo/-n3lo" command is given, overwrite order_min and order_max accordingly
    if (log_val_opt("-lo")) then
       order_min = 1
       order_max = 1
       order = '  lo'
    elseif (log_val_opt("-nlo")) then
       order_min = 1
       order_max = 2
       order = ' nlo'
    elseif (log_val_opt("-nnlo")) then
       order_min = 1
       order_max = 3
       order = 'nnlo'
    elseif (log_val_opt("-n3lo")) then
       order_min = 1
       order_max = 4
       order = 'n3lo'
    else
       if(order_max.eq.1) order = '  lo'
       if(order_max.eq.2) order = ' nlo'
       if(order_max.eq.3) order = 'nnlo'
       if(order_max.eq.4) order = 'n3lo'
    endif
    if(order_min.ne.1) separate_orders = .true. ! In this case we always need to separate the orders
    NC = log_val_opt("-NC",.true.)
    CC = log_val_opt("-CC",.false.)
    noZ = .not.log_val_opt ("-includeZ",.false.)
    Zonly = .false.
    intonly = .false.
    if(.not.noZ) Zonly = log_val_opt ("-Zonly",.false.)
    if(.not.noZ) intonly = log_val_opt ("-intonly",.false.)
    if(Zonly.and.intonly) stop 'Cannot run with both Z and interference ONLY flags'
    if(.not.NC.and..not.CC) stop 'Need to run with either or/both of NC and CC'
    p2b = log_val_opt("-p2b",.false.)
    inclusive = .not.p2b 
    if(.not.noZ.and.p2b) stop 'Cannot do Z in p2b yet'
    if(CC.and.p2b) stop 'Cannot do CC in p2b yet'
    if(order_max.ge.4.and.p2b) stop 'Cannot run p2b at N3LO yet'
    if(vnf.and.p2b) stop 'Cannot run p2b with variable flavour'
    outname      = string_val_opt("-out", "") ! Overwite the prefix of the file name
    prefix       = string_val_opt("-prefix", "") ! Overwite the prefix of the file name
    if(prefix.eq.'') prefix       = string_val_opt("-out", "") ! Overwite the prefix of the file name
    do_analysis = .not.log_val_opt("-no-analysis",.false.)
    if(p2b.and..not.do_analysis) stop 'Should really be doing an analysis with p2b'

    ! Parameters dealing with scale variations
    scale_choice = int_val_opt ('-scale-choice',1) ! 1: Use Q. 0: Use MZ. For now fixed.
    xmuf         = dble_val_opt("-xmuf",1.0_dp)
    xmur         = dble_val_opt("-xmur",1.0_dp)
    pdfname      = string_val_opt("-pdf", "")

    toy_Q0       = dble_val_opt("-toyQ0",-1d0)
    Q0pdf        = dble_val_opt("-Q0pdf",-1d0)
    xmuR_PDF     = dble_val_opt("-xmuRPDF",1d0)
    if(pdfname.eq.'') then
       if(toy_Q0.lt.0d0) then
          call help_message
          print*, '-pdf must be specified!'
          call exit()
       endif
    else
       ! Initialise PDF
       call initPDFSetByName(pdfname)
       call getQ2min(0,Q2minPDF)
    endif

    nmempdf      = int_val_opt ("-nmempdf",0)
    pdfuncert    = log_val_opt ("-pdfuncert")
    alphasuncert = log_val_opt ("-alphasuncert")
    if(alphasuncert.and..not.pdfuncert) stop "Need -pdfuncert to run with -alphasuncert"
    scaleuncert  = log_val_opt ("-scaleuncert")
    if(scaleuncert.and.scale_choice.lt.2) scale_choice = 2
    if(scale_choice.ge.2) separate_orders = .true.
    if(scaleuncert.and.vnf) then
       print*, 'Cannot currently do automatic scale uncertainties and&
            & vnf. Please run the scale choices individually like&
            & this:'
       print*, './disorder -xmur 1.0 -xmuf 1.0 -vnf '
       stop
    endif
    Nscales =1
    if(scaleuncert) Nscales = 7
    scalestr(1) = '_μR_1.0_μF_1.0'
    scalestr(2) = '_μR_2.0_μF_2.0'
    scalestr(3) = '_μR_0.5_μF_0.5'
    scalestr(4) = '_μR_1.0_μF_2.0'
    scalestr(5) = '_μR_1.0_μF_0.5'
    scalestr(6) = '_μR_2.0_μF_1.0'
    scalestr(7) = '_μR_0.5_μF_1.0'


    if(p2b.and.pdfuncert) stop "Cannot do pdf uncertainties with p2b"

    ! Some parameters setting up the VEGAS run
    readin       = log_val_opt ("-readingrid")
    ncall1       = int_val_opt ("-ncall1",10000)
    ncall2       = int_val_opt ("-ncall2",100000)
    it1          = int_val_opt("-it1",1)
    itmx1        = int_val_opt ("-itmx1",3)
    itmx2        = 1!int_val_opt ("-itmx2",1)
    if((int_val_opt("-iseed",1).gt.1) .and. &
         & (int_val_opt("-irseq",1).gt.1)) stop 'Cannot have both rseq and  &
         &        iseed on commandline'
    iseed        = int_val_opt ("-iseed",1)
    if(iseed.eq.1) then
       iseed        = int_val_opt ("-rseq",1)
    endif
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
    ingridfile    =trim(prefix)//'grids-'//trim(adjustl(order))//'-'//seedstr//'.dat'
    outgridfile   =ingridfile
    outgridtopfile=trim(prefix)//'grids-'//trim(adjustl(order))//'-'//seedstr//'.top'
    ilast=0

    ! Setup the phase space
    !     Read in bounds on x,y,Q2
    Q = dble_val_opt("-Q",-1.0_dp)
    if(Q.lt.0.0_dp) then
       if(dble_val_opt("-Qmin",-1d0) .gt. 0d0) then
          Q2min = (dble_val_opt("-Qmin",sqrt(Q2minPDF)))**2
       else
          Q2min = dble_val_opt("-Q2min",Q2minPDF)
       endif
       if(dble_val_opt("-Qmax",-1d0) .gt. 0d0) then
          Q2max = (dble_val_opt("-Qmax",1d100))**2
       else
          Q2max = dble_val_opt("-Q2max",1d200)
       endif
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
    Eh = dble_val_opt("-Ehad",820.0_dp)
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

    novegas = .false.
    if(Q2min.eq.Q2max.and.xmin.eq.xmax) then
       ymin = Q2min/(s*xmin)
       ymax = ymin
       novegas = .true.
    elseif(Q2min.eq.Q2max.and.ymin.eq.ymax) then
       xmin = Q2min/(s*ymin)
       xmax = xmin
       novegas = .true.
    elseif(xmin.eq.xmax.and.ymin.eq.ymax) then
       Q2min = s*xmin*ymin
       Q2max = Q2min
       novegas = .true.
    else if(xmin.eq.xmax) then
       if(xmin*ymin*s .gt. Q2min) Q2min = xmin*ymin*s
       if(xmax*ymax*s .lt. Q2max) Q2max = xmax*ymax*s
    else
       if(xmin*ymax*s .lt. Q2min) xmin = Q2min/(ymax*s)
       if(xmax*ymin*s .gt. Q2max) xmax = Q2max/(ymin*s)
    endif
    !         if(s*xmin*ymin .gt. Q2min) Q2min = s*xmin*ymin
    !         if(s*xmax*ymax .gt. Q2max) Q2max = s*xmax*ymax
    sqrts = sqrt(s * xmax)  

    ! Some parameters for DISENT
    NPOW1  = dble_val_opt("-npow1",2.0_dp)
    NPOW2  = dble_val_opt("-npow2",4.0_dp)
    CUTOFF = dble_val_opt("-cutoff",1d-8)

    ! Sensible initial Qmin for the PDF
    Qmin = 1.0_dp

    sigma_all_scales = 0.0_dp ! Initialise
    NC_reduced_sigma = 0.0_dp
    CC_reduced_sigma = 0.0_dp
    NC_reduced_dsigma = 0.0_dp
    CC_reduced_dsigma = 0.0_dp

    ! For hoppetStartExtended. Could think of putting on commandline...
    ! Streamlined initialization
    ! including  parameters for x-grid
    order_hoppet = -6 
    ymax_hoppet  = 16.0_dp
    dy    = 0.05_dp  ! dble_val_opt("-dy",0.1_dp)
    dlnlnQ = dy/4.0_dp
    nloop = order_max
    minQval = min(xmuF*Qmin, Qmin)
    maxQval = max(xmuF*sqrts, sqrts)
    if(scale_choice.eq.4) maxQval = max(xmuF*sqrt(s), sqrt(s))
    scale_choice_hoppet = min(2,scale_choice)

    if (.not.CheckAllArgsUsed(0)) then
       call help_message
       call exit()
    endif

  end subroutine set_parameters

  !----------------------------------------------------------------------
  real(dp) function alphasLocal(muR)
    real(dp), intent(in) :: muR
    real(dp) :: muR_lcl, alphasPDF
    muR_lcl = max(muR,Qmin)
    if(toy_Q0 < 0d0) then
       ! we use alphas from the LHAPDF PDF
       alphasLocal = alphasPDF(muR_lcl)
    else
       ! we use alphas from HOPPET
       alphasLocal = CouplingValue(coupling, muR_lcl)
    endif
  end function alphasLocal


  subroutine print_header(idev)
    implicit none
    integer, intent(in) :: idev
    real(dp) :: Qmn, Qmx
    write(idev,'(a)',advance='no') ' # '
    call time_stamp(idev)
    write(idev,*) '#'//trim(command_line())
    write(idev,*) '# ----------------------------------------------------------'
    write(idev,'(a,a)',advance='no') ' # Doing DIS @ ', trim(order)
    if(.not.positron) write(idev,*) ' (e^- + p)'
    if(positron) write(idev,*) ' (e^+ + p)'
    if(.not.p2b) then
       write(idev,*) '# Inclusively in radiation'
    elseif(p2b) then
       write(idev,*) '# Using DISENT with projection-to-Born'
    endif
    if(CC) write(idev,*) '# Including charged current'
    if(NC) then
       write(idev,*) '# Including neutral current'
       if(noZ)  write(idev,*) '# With γ only'
       if(.not.noZ.and.(.not.Zonly).and.(.not.intonly))  write(idev,*) '# With γ/Z'
       if(Zonly)  write(idev,*) '# With Z only'
       if(intonly)  write(idev,*) '# With γ/Z interference only'
    else
       write(idev,*) '# And no neutral current'
    endif
    Qmn = dsqrt(Q2min)
    Qmx = dsqrt(Q2max)
    write(idev,'(a,F14.7,F14.7)') ' # xmin, xmax:      ', xmin, xmax
    write(idev,'(a,F14.7,F14.7)') ' # ymin, ymax:      ', ymin, ymax
    write(idev,'(a,F14.7,F14.7,a)') ' # Qmin, Qmax:      ', Qmn, Qmx, ' GeV'
    write(idev,'(a,F14.7,a)') ' # Electron energy: ', El, ' GeV'
    write(idev,'(a,F14.7,a)') ' # Proton energy:   ', Eh, ' GeV'
    write(idev,'(a,F14.7,a)') ' # COM energy:      ', S,  ' GeV^2'
    if(toy_Q0 < zero) write(idev,'(a,a)') ' # PDF:             ', trim(adjustl(pdfname))
    if(toy_Q0 > zero) write(idev,*) ' # PDF:             ', 'LHA toy PDF initialised at', toy_Q0, 'GeV'
    write(idev,'(a,F14.7,a)') ' # MZ:              ', MZ, ' GeV'
    write(idev,'(a,F14.7,a)') ' # MW:              ', MW, ' GeV'
    if(.not.vnf)   write(idev,'(a,I14)') ' # nf:              ', nflav
    if(vnf)   write(idev,'(a,a)') ' # nf:              ', 'variable'
    write(idev,'(a,F14.7)') ' # CA:              ', CAlcl
    write(idev,'(a,F14.7)') ' # CF:              ', CFlcl
    write(idev,'(a,F14.7)') ' # TR:              ', TRlcl
    write(idev,'(a,F14.7)') ' # αS(MZ):          ', alphasLocal(MZ)
    write(idev,'(a,F14.7)') ' # 1/αEM:           ', 1.0_dp/alpha_em
    write(idev,'(a,E14.7,a)') ' # GF:                  ', GF,  ' GeV^-2'
    write(idev,'(a,F14.7)') ' # sin(θ_W)^2:      ', sin_thw_sq

    write(idev,*) '# ----------------------------------------------------------'

  end subroutine print_header

  subroutine welcome_message
    write(0,'(a)') '-----------------------------------------------------------'
    write(0,'(a)') '               Welcome to disorder v. 1.0.0                '
    write(0,'(a)') '        Written by Alexander Karlberg (2023-2024)          '
    write(0,'(a)') '                                                           '
    write(0,'(a)') ' It is made available under the GNU public license,        '
    write(0,'(a)') ' with the additional request that if you use it or any     '
    write(0,'(a)') ' derivative of it in scientific work then you should cite: '
    write(0,'(a)') ' A. Karlberg (arXiv:2401.16964).                           '
    write(0,'(a)') '                                                           '
    write(0,'(a)') ' You are also encouraged to cite HOPPET, the original      '
    write(0,'(a)') ' references,for LO, NLO and NNLO splitting functions, the  '
    write(0,'(a)') ' QCD 1, 2, 3 and 4 loop beta functions and the coupling and'
    write(0,'(a)') ' PDF mass threshold matching functions. You are furthermore' 
    write(0,'(a)') ' encouraged to cite the LO, NLO, NNLO, and N3LO coefficient'
    write(0,'(a)') ' functions and the disent references.                      '
    write(0,'(a)') '-----------------------------------------------------------'
  end subroutine welcome_message

  subroutine help_message
    call welcome_message
    write(0,'(a)') ' Some common flags to use are (default values in [] and ()   '
    write(0,'(a)') ' means that the flag takes an input (dble/int/string).       '
    write(0,'(a)') ' Values are in GeV typically). -pdf is mandatory :           '
    write(0,'(a)') '                                                             '
    write(0,'(a)') ' # DIS setup:                                                '
    write(0,'(a)') ' -Q (dble) : Specify fixed Q OR                              '
    write(0,'(a)') ' -Qmin (dble) [QminPDF] : Specify minimum Q                  '
    write(0,'(a)') ' -Qmax (dble) : Specify maximum Q                            '
    write(0,'(a)') ' -Q2min (dble) [Q2minPDF] : Specify minimum Q2               '
    write(0,'(a)') ' -Q2max (dble) : Specify maximum Q2                          '
    write(0,'(a)') ' -x (dble) : Specify fixed x OR                              '
    write(0,'(a)') ' -xmin (dble) : Specify minimum x                            '
    write(0,'(a)') ' -xmax (dble) : Specify maximum x                            '
    write(0,'(a)') ' -y (dble) : Specify fixed y OR                              '
    write(0,'(a)') ' -ymin (dble) : Specify minimum y                            '
    write(0,'(a)') ' -ymax (dble) : Specify maximum y                            '
    write(0,'(a)') ' -CC [false] : Include charged current processes             '
    write(0,'(a)') ' -NC [true] : Include neutral current processes              '
    write(0,'(a)') ' -includeZ [false] : Include Z and interferences             '
    write(0,'(a)') ' -positron [false] : Incoming positron                       '
    write(0,'(a)') ' -Elep (dble) [27.5] : Energy of lepton in lab frame         '
    write(0,'(a)') ' -Ehad (dble) [820.0] : Energy of hadron in lab frame        '
    write(0,'(a)') '                                                             '
    write(0,'(a)') ' # QCD setup:                                                '
    write(0,'(a)') ' -pdf : LHAPDF name (e.g. NNPDF30_nnlo_as_0118_hera)         '
    write(0,'(a)') ' -lo/-nlo/-nnlo/-n3lo [-nnlo]: Run at LO/NLO/NNLO/N3LO       '
    write(0,'(a)') ' -scale-choice (int) [1]: Central scale: 0: MZ, 1: Q,        '
    write(0,'(a)') '                          3: Q*sqrt(1-y), 4: Q*(1-x)/x       '
    write(0,'(a)') ' -scaleuncert [false]: Do 7-point scale variation around Q   '
    write(0,'(a)') ' -pdfuncert [false] : Compute pdf uncertainties              '
    write(0,'(a)') ' -p2b [false] : Turn on disent and projection-to-Born        '
    write(0,'(a)') '                                                             '
    write(0,'(a)') ' # EW parameters:                                            '
    write(0,'(a)') ' -one-over-alpha (dble) [137.] : 1/αEM                       '
    write(0,'(a)') ' -mz (dble) [91.1876] : mass of Z boson                      '
    write(0,'(a)') ' -mw (dble) [80.398]  : mass of W boson                      '
    write(0,'(a)') '                                                             '
    write(0,'(a)') ' # Run parameters:                                           '
    write(0,'(a)') ' -prefix (string) : Add a prefix to all output               '
    write(0,'(a)') ' -ncall1 (int) [100000] : Number of calls to VEGAS warmup    '
    write(0,'(a)') ' -itmx1 (int) [3]  : Number of iterations in VEGAS warmup    '
    write(0,'(a)') ' -ncall2 (int) [100000] : Number of calls to VEGAS production'
    write(0,'(a)') ' -iseed/-rseq (int) [1] : The seed                           '
    write(0,'(a)') ' -no-analysis [false] : Turn off the analysis                ' 
    write(0,'(a)') ' -help : Print this help message                             '
    write(0,'(a)') '                                                             '
    write(0,'(a)') ' # Toy PDF parameters:                                       '
    write(0,'(a)') ' -toy_Q0 (dble) [-1.0]: If > 0 then uses the HERALHC initial '
    write(0,'(a)') ' condition at that scale and evolves using Hoppet            '
    write(0,'(a)') ' -Q0pdf (dble) [-1.0]: If > 0 then read in LHAPDF at that    '
    write(0,'(a)') ' scale and evolve using Hoppet                               '
    write(0,'(a)') '                                                             '
    write(0,'(a)') ' More specialised flags and helpful information can be       '
    write(0,'(a)') ' found in src/mod_parameters.f90. Logical flags are set to   '
    write(0,'(a)') ' true if present on the command line (-p2b turns on p2b).    '
    write(0,'(a)') ' They can be set to false by prefixing with "no", for        '
    write(0,'(a)') ' instance like this : -noNC.                                 '
    write(0,'(a)') '                                                             '
    write(0,'(a)') ' Example command line:                                       '
    write(0,'(a)') '                                                             '
    write(0,'(a)') ' ./disorder -pdf MSHT20nnlo_as118 -Q 20.0 -noNC -CC          '
  end subroutine help_message

end module mod_parameters
