set terminal pdf enhanced font "Latin Modern Roman,36" size 30cm,21cm
set datafile fortran

set output 'sigma-ratios.pdf'
set macros

#linetype 1,  linecolor rgb "dark-violet"  linewidth 1.000 dashtype solid pointtype 1 pointsize default
#linetype 2,  linecolor rgb "#009e73"  linewidth 1.000 dashtype solid pointtype 2 pointsize default
#linetype 3,  linecolor rgb "#56b4e9"  linewidth 1.000 dashtype solid pointtype 3 pointsize default
#linetype 4,  linecolor rgb "#e69f00"  linewidth 1.000 dashtype solid pointtype 4 pointsize default
#linetype 5,  linecolor rgb "#f0e442"  linewidth 1.000 dashtype solid pointtype 5 pointsize default
#linetype 6,  linecolor rgb "#0072b2"  linewidth 1.000 dashtype solid pointtype 6 pointsize default
#linetype 7,  linecolor rgb "#e51e10"  linewidth 1.000 dashtype solid pointtype 7 pointsize default
#linetype 8,  linecolor rgb "black" 

NLOline   = "lines lc rgb '#009e73' lw 2.0"
NLOlinethin   = "lines lc rgb '#009e73' lw 0.5"
NLOfill   = "filledcurves lc rgb '#009e73' lw 2.0 fs transparent solid 0.2"

NNLOline  = "lines lc rgb 'dark-violet' lw 2.0"
NNLOlinethin  = "lines lc rgb 'dark-violet' lw 0.5"
NNLOfill  = "filledcurves lc rgb 'dark-violet' lw 2.0 fs transparent solid 0.2"

N3LOline = "lines lc rgb '#009e73' lw 2.0"
N3LOlinethin = "lines lc rgb '#009e73' lw 0.5"
N3LOfill = "filledcurves lc rgb '#009e73' lw 0.2 fs transparent solid 0.2"

aN3LOline = "lines lc rgb '#e51e10' lw 2.0"
aN3LOlinethin = "lines lc rgb '#e51e10' lw 0.5"
aN3LOfill  = "filledcurves lc rgb '#e51e10' lw 2.0 fs transparent solid 0.2"


nlo='<paste NC_Q10_disorder_nlo_central.dat NC_Q10_disorder_nlo_max.dat NC_Q10_disorder_nlo_min.dat NC_Q10_disorder_nlo_pdfuncert.dat NC_Q10_disorder_nnlo_central.dat'
nnlo='<paste NC_Q10_disorder_nnlo_central.dat NC_Q10_disorder_nnlo_max.dat NC_Q10_disorder_nnlo_min.dat NC_Q10_disorder_nnlo_pdfuncert.dat NC_Q10_disorder_nnlo_central.dat'
n3lo='<paste NC_Q10_disorder_n3lo_central.dat NC_Q10_disorder_n3lo_max.dat NC_Q10_disorder_n3lo_min.dat NC_Q10_disorder_n3lo_pdfuncert.dat NC_Q10_disorder_nnlo_central.dat'
an3lo='<paste NC_an3lo_Q10_disorder_n3lo_central.dat NC_an3lo_Q10_disorder_n3lo_max.dat NC_an3lo_Q10_disorder_n3lo_min.dat NC_an3lo_Q10_disorder_n3lo_pdfuncert.dat NC_Q10_disorder_nnlo_central.dat'


reset

set mxtics
set mytics
set grid
#set log y
#set format y "10^{%T}"
set xrange [*:log(1.)]
set yrange [0.95:1.10]
set key at  -5.2,1.095

set title 'σ(Reduced NC) Q = 10 GeV'

#set label 2 '4 GeV^2 < Q^2 < 5 GeV^2' font "Latin Modern Roman,28" at -8.9,1.7
#set label 28 'e^+ (27.6 GeV) p (920 GeV) → e^+ + X' font "Latin Modern Roman,28" at -1, 0.4 right
#set label 29 'NNPDF30\_nnlo\_as\_0118\_hera' font "Latin Modern Roman,28" at -1, 0.32 right
#set label 30 'μ_R = μ_F = Q' font "Latin Modern Roman,28" at -1, 0.24 right
#set label 31 '7-point variation' font "Latin Modern Roman,28" at -1, 0.16 right

set ylabel 'Ratio to NNLO'
set xlabel 'log[x]'

ii=0
emax=8

plot 1 w @NNLOline notitle,\
     nnlo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @NNLOfill title 'NNLO',\
     nnlo i ii u (($1+$2)/2.):(($7+$16)/$19) w @NNLOlinethin not,\
     nnlo i ii u (($1+$2)/2.):(($11-$16)/$19) w @NNLOlinethin not,\
     n3lo i ii u (($1+$2)/2.):($3/$19):4 w @N3LOline not,\
     n3lo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @N3LOfill title 'N^3LO',\
     n3lo i ii u (($1+$2)/2.):(($7+$16)/$19) w @N3LOlinethin not,\
     n3lo i ii u (($1+$2)/2.):(($11-$16)/$19) w @N3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):($3/$19):4 w @aN3LOline not,\
     an3lo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @aN3LOfill title 'aN^3LO',\
     an3lo i ii u (($1+$2)/2.):(($7+$16)/$19) w @aN3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):(($11-$16)/$19) w @aN3LOlinethin not

#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @LHEFfill title 'pythia8 dipole',\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @LHEFlinethin not,\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @LHEFlinethin not,\
#     lhef every ::emax i ii u (($1+$2)/2.):3:4 w @LHEFline not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @PY8fill title 'pythia8 default',\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @PY8linethin not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @PY8linethin not,\
#     pythia8bad every ::emax i ii u (($1+$2)/2.):3:4 w @PY8line not
#     nlo i ii u (($1+$2)/2.):($3/$19):4 w @NLOline not,\
#     nlo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @NLOfill title 'NLO',\
#     nlo i ii u (($1+$2)/2.):(($7+$16)/$19) w @NLOlinethin not,\
#     nlo i ii u (($1+$2)/2.):(($11-$16)/$19) w @NLOlinethin not,\

reset

set mxtics
set mytics
set grid
#set log y
#set format y "10^{%T}"
set xrange [*:*]
set yrange [0.9:1.2]
set key at  -5.3,1.19

set title 'σ(Reduced NC) Q = 10 GeV'

#set label 2 '4 GeV^2 < Q^2 < 5 GeV^2' font "Latin Modern Roman,28" at -8.9,1.7
#set label 28 'e^+ (27.6 GeV) p (920 GeV) → e^+ + X' font "Latin Modern Roman,28" at -1, 0.4 right
#set label 29 'NNPDF30\_nnlo\_as\_0118\_hera' font "Latin Modern Roman,28" at -1, 0.32 right
#set label 30 'μ_R = μ_F = Q' font "Latin Modern Roman,28" at -1, 0.24 right
#set label 31 '7-point variation' font "Latin Modern Roman,28" at -1, 0.16 right

set ylabel 'Ratio to NNLO'
set xlabel 'x'

ii=1
emax=8

plot 1 w @NNLOline notitle,\
     nnlo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @NNLOfill title 'NNLO',\
     nnlo i ii u (($1+$2)/2.):(($7+$16)/$19) w @NNLOlinethin not,\
     nnlo i ii u (($1+$2)/2.):(($11-$16)/$19) w @NNLOlinethin not,\
     n3lo i ii u (($1+$2)/2.):($3/$19):4 w @N3LOline not,\
     n3lo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @N3LOfill title 'N3LO',\
     n3lo i ii u (($1+$2)/2.):(($7+$16)/$19) w @N3LOlinethin not,\
     n3lo i ii u (($1+$2)/2.):(($11-$16)/$19) w @N3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):($3/$19):4 w @aN3LOline not,\
     an3lo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @aN3LOfill title 'aN3LO',\
     an3lo i ii u (($1+$2)/2.):(($7+$16)/$19) w @aN3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):(($11-$16)/$19) w @aN3LOlinethin not

#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @LHEFfill title 'pythia8 dipole',\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @LHEFlinethin not,\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @LHEFlinethin not,\
#     lhef every ::emax i ii u (($1+$2)/2.):3:4 w @LHEFline not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @PY8fill title 'pythia8 default',\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @PY8linethin not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @PY8linethin not,\
#     pythia8bad every ::emax i ii u (($1+$2)/2.):3:4 w @PY8line not
#     nlo i ii u (($1+$2)/2.):($3/$19):4 w @NLOline not,\
#     nlo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @NLOfill title 'NLO',\
#     nlo i ii u (($1+$2)/2.):(($7+$16)/$19) w @NLOlinethin not,\
#     nlo i ii u (($1+$2)/2.):(($11-$16)/$19) w @NLOlinethin not,\

nlo='<paste NC_x0.01_disorder_nlo_central.dat NC_x0.01_disorder_nlo_max.dat NC_x0.01_disorder_nlo_min.dat NC_x0.01_disorder_nlo_pdfuncert.dat NC_x0.01_disorder_nnlo_central.dat'
nnlo='<paste NC_x0.01_disorder_nnlo_central.dat NC_x0.01_disorder_nnlo_max.dat NC_x0.01_disorder_nnlo_min.dat NC_x0.01_disorder_nnlo_pdfuncert.dat NC_x0.01_disorder_nnlo_central.dat'
n3lo='<paste NC_x0.01_disorder_n3lo_central.dat NC_x0.01_disorder_n3lo_max.dat NC_x0.01_disorder_n3lo_min.dat NC_x0.01_disorder_n3lo_pdfuncert.dat NC_x0.01_disorder_nnlo_central.dat'
an3lo='<paste NC_an3lo_x0.01_disorder_n3lo_central.dat NC_an3lo_x0.01_disorder_n3lo_max.dat NC_an3lo_x0.01_disorder_n3lo_min.dat NC_an3lo_x0.01_disorder_n3lo_pdfuncert.dat NC_x0.01_disorder_nnlo_central.dat'

reset

set mxtics
set mytics
set grid
#set log y
#set format y "10^{%T}"
set xrange [0.5:*]
set yrange [0.8:1.2]
set key at  3.4,1.19

set title 'σ(Reduced NC) x = 0.01'

#set label 2 '4 GeV^2 < Q^2 < 5 GeV^2' font "Latin Modern Roman,28" at -8.9,1.7
#set label 28 'e^+ (27.6 GeV) p (920 GeV) → e^+ + X' font "Latin Modern Roman,28" at -1, 0.4 right
#set label 29 'NNPDF30\_nnlo\_as\_0118\_hera' font "Latin Modern Roman,28" at -1, 0.32 right
#set label 30 'μ_R = μ_F = Q' font "Latin Modern Roman,28" at -1, 0.24 right
#set label 31 '7-point variation' font "Latin Modern Roman,28" at -1, 0.16 right

set ylabel 'Ratio to NNLO'
set xlabel 'log[Q/GeV]'

ii=2
emax=8

plot 1 w @NNLOline notitle,\
     nnlo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @NNLOfill title 'NNLO',\
     nnlo i ii u (($1+$2)/2.):(($7+$16)/$19) w @NNLOlinethin not,\
     nnlo i ii u (($1+$2)/2.):(($11-$16)/$19) w @NNLOlinethin not,\
     n3lo i ii u (($1+$2)/2.):($3/$19):4 w @N3LOline not,\
     n3lo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @N3LOfill title 'N^3LO',\
     n3lo i ii u (($1+$2)/2.):(($7+$16)/$19) w @N3LOlinethin not,\
     n3lo i ii u (($1+$2)/2.):(($11-$16)/$19) w @N3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):($3/$19):4 w @aN3LOline not,\
     an3lo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @aN3LOfill title 'aN^3LO',\
     an3lo i ii u (($1+$2)/2.):(($7+$16)/$19) w @aN3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):(($11-$16)/$19) w @aN3LOlinethin not

#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @LHEFfill title 'pythia8 dipole',\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @LHEFlinethin not,\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @LHEFlinethin not,\
#     lhef every ::emax i ii u (($1+$2)/2.):3:4 w @LHEFline not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @PY8fill title 'pythia8 default',\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @PY8linethin not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @PY8linethin not,\
#     pythia8bad every ::emax i ii u (($1+$2)/2.):3:4 w @PY8line not
#     nlo i ii u (($1+$2)/2.):($3/$19):4 w @NLOline not,\
#     nlo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @NLOfill title 'NLO',\
#     nlo i ii u (($1+$2)/2.):(($7+$16)/$19) w @NLOlinethin not,\
#     nlo i ii u (($1+$2)/2.):(($11-$16)/$19) w @NLOlinethin not,\

reset

set mxtics
set mytics
set grid
#set log y
#set format y "10^{%T}"
set xrange [*:*]
set yrange [0.9:1.2]
set key at  -5.3,1.19

set title 'σ(Reduced NC) x = 0.01'

#set label 2 '4 GeV^2 < Q^2 < 5 GeV^2' font "Latin Modern Roman,28" at -8.9,1.7
#set label 28 'e^+ (27.6 GeV) p (920 GeV) → e^+ + X' font "Latin Modern Roman,28" at -1, 0.4 right
#set label 29 'NNPDF30\_nnlo\_as\_0118\_hera' font "Latin Modern Roman,28" at -1, 0.32 right
#set label 30 'μ_R = μ_F = Q' font "Latin Modern Roman,28" at -1, 0.24 right
#set label 31 '7-point variation' font "Latin Modern Roman,28" at -1, 0.16 right

set ylabel 'Ratio to NNLO'
set xlabel 'Q [GeV]'

ii=3
emax=8

plot 1 w @NNLOline notitle,\
     nnlo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @NNLOfill title 'NNLO',\
     nnlo i ii u (($1+$2)/2.):(($7+$16)/$19) w @NNLOlinethin not,\
     nnlo i ii u (($1+$2)/2.):(($11-$16)/$19) w @NNLOlinethin not,\
     n3lo i ii u (($1+$2)/2.):($3/$19):4 w @N3LOline not,\
     n3lo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @N3LOfill title 'N3LO',\
     n3lo i ii u (($1+$2)/2.):(($7+$16)/$19) w @N3LOlinethin not,\
     n3lo i ii u (($1+$2)/2.):(($11-$16)/$19) w @N3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):($3/$19):4 w @aN3LOline not,\
     an3lo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @aN3LOfill title 'aN3LO',\
     an3lo i ii u (($1+$2)/2.):(($7+$16)/$19) w @aN3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):(($11-$16)/$19) w @aN3LOlinethin not

#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @LHEFfill title 'pythia8 dipole',\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @LHEFlinethin not,\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @LHEFlinethin not,\
#     lhef every ::emax i ii u (($1+$2)/2.):3:4 w @LHEFline not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @PY8fill title 'pythia8 default',\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @PY8linethin not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @PY8linethin not,\
#     pythia8bad every ::emax i ii u (($1+$2)/2.):3:4 w @PY8line not
#     nlo i ii u (($1+$2)/2.):($3/$19):4 w @NLOline not,\
#     nlo i ii u (($1+$2)/2.):(($7+$16)/$19):(($11-$16)/$19) w @NLOfill title 'NLO',\
#     nlo i ii u (($1+$2)/2.):(($7+$16)/$19) w @NLOlinethin not,\
#     nlo i ii u (($1+$2)/2.):(($11-$16)/$19) w @NLOlinethin not,\


nlo='<paste NC_Q10_disorder_nlo_central.dat NC_Q10_disorder_nlo_max.dat NC_Q10_disorder_nlo_min.dat NC_Q10_disorder_nlo_pdfuncert.dat NC_Q10_disorder_nnlo_central.dat'
nnlo='<paste NC_Q10_disorder_nnlo_central.dat NC_Q10_disorder_nnlo_max.dat NC_Q10_disorder_nnlo_min.dat NC_Q10_disorder_nnlo_pdfuncert.dat NC_Q10_disorder_nnlo_central.dat'
n3lo='<paste NC_Q10_disorder_n3lo_central.dat NC_Q10_disorder_n3lo_max.dat NC_Q10_disorder_n3lo_min.dat NC_Q10_disorder_n3lo_pdfuncert.dat NC_Q10_disorder_nnlo_central.dat'
an3lo='<paste NC_an3lo_Q10_disorder_n3lo_central.dat NC_an3lo_Q10_disorder_n3lo_max.dat NC_an3lo_Q10_disorder_n3lo_min.dat NC_an3lo_Q10_disorder_n3lo_pdfuncert.dat NC_Q10_disorder_nnlo_central.dat'


reset

set mxtics
set mytics
set grid
#set log y
#set format y "10^{%T}"
set xrange [*:log(1.)]
set yrange [0.9:1.2]
set key at  -5.2,1.19

set title 'σ(Reduced NC) Q = 10 GeV'

#set label 2 '4 GeV^2 < Q^2 < 5 GeV^2' font "Latin Modern Roman,28" at -8.9,1.7
#set label 28 'e^+ (27.6 GeV) p (920 GeV) → e^+ + X' font "Latin Modern Roman,28" at -1, 0.4 right
#set label 29 'NNPDF30\_nnlo\_as\_0118\_hera' font "Latin Modern Roman,28" at -1, 0.32 right
#set label 30 'μ_R = μ_F = Q' font "Latin Modern Roman,28" at -1, 0.24 right
#set label 31 '7-point variation' font "Latin Modern Roman,28" at -1, 0.16 right

set ylabel 'Ratio to NNLO'
set xlabel 'log[x]'

ii=0
emax=8

plot 1 w @NNLOline notitle,\
     nnlo i ii u (($1+$2)/2.):(($7+0.*$16)/$19):(($11-0.*$16)/$19) w @NNLOfill title 'NNLO',\
     nnlo i ii u (($1+$2)/2.):(($7+0.*$16)/$19) w @NNLOlinethin not,\
     nnlo i ii u (($1+$2)/2.):(($11-0.*$16)/$19) w @NNLOlinethin not,\
     n3lo i ii u (($1+$2)/2.):($3/$19):4 w @N3LOline not,\
     n3lo i ii u (($1+$2)/2.):(($7+0.*$16)/$19):(($11-0.*$16)/$19) w @N3LOfill title 'N^3LO',\
     n3lo i ii u (($1+$2)/2.):(($7+0.*$16)/$19) w @N3LOlinethin not,\
     n3lo i ii u (($1+$2)/2.):(($11-0.*$16)/$19) w @N3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):($3/$19):4 w @aN3LOline not,\
     an3lo i ii u (($1+$2)/2.):(($7+0.*$16)/$19):(($11-0.*$16)/$19) w @aN3LOfill title 'aN^3LO',\
     an3lo i ii u (($1+$2)/2.):(($7+0.*$16)/$19) w @aN3LOlinethin not,\
     an3lo i ii u (($1+$2)/2.):(($11-0.*$16)/$19) w @aN3LOlinethin not

#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @LHEFfill title 'pythia8 dipole',\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @LHEFlinethin not,\
#     lhefmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @LHEFlinethin not,\
#     lhef every ::emax i ii u (($1+$2)/2.):3:4 w @LHEFline not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7:3 w @PY8fill title 'pythia8 default',\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):7 w   @PY8linethin not,\
#     pythia8badmaxmin every ::emax i ii u (($1+$2)/2.):3 w   @PY8linethin not,\
#     pythia8bad every ::emax i ii u (($1+$2)/2.):3:4 w @PY8line not
#     nlo i ii u (($1+$2)/2.):($3/$19):4 w @NLOline not,\
#     nlo i ii u (($1+$2)/2.):(($7+0.*$16)/$19):(($11-0.*$16)/$19) w @NLOfill title 'NLO',\
#     nlo i ii u (($1+$2)/2.):(($7+0.*$16)/$19) w @NLOlinethin not,\
#     nlo i ii u (($1+$2)/2.):(($11-0.*$16)/$19) w @NLOlinethin not,\

set output