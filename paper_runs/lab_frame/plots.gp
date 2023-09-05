set terminal pdf enhanced font "Latin Modern Roman,36" size 30cm,21cm
set datafile fortran

set output 'lab-frame.pdf'
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

LOline  = "xerrorbars lc rgb 'dark-violet' lw 2.0"
LOlinethin  = "lines lc rgb 'dark-violet' lw 0.5"
LOfill  = "boxxyerrorbars lc rgb 'dark-violet' lw 2.0 fs transparent solid 0.2"

NLOline = "xerrorbars lc rgb '#009e73' lw 2.0"
NLOlinethin = "lines lc rgb '#009e73' lw 0.5"
NLOfill = "boxxyerrorbars lc rgb '#009e73' lw 0.2 fs transparent solid 0.2"

NNLOline = "xerrorbars lc rgb '#e51e10' lw 2.0"
NNLOlinethin = "lines lc rgb '#e51e10' lw 0.5"
NNLOfill  = "boxxyerrorbars lc rgb '#e51e10' lw 2.0 fs transparent solid 0.2"


lo='<paste NC_EIC_disorder_lo_central.dat NC_EIC_disorder_lo_max.dat NC_EIC_disorder_lo_min.dat '
nlo='<paste NC_EIC_disorder_nlo_central.dat NC_EIC_disorder_nlo_max.dat NC_EIC_disorder_nlo_min.dat '
nnlo='<paste NC_EIC_disorder_nnlo_central.dat NC_EIC_disorder_nnlo_max.dat NC_EIC_disorder_nnlo_min.dat '


reset

set mxtics
set mytics
set grid
set log y
set format y "10^{%T}"
set xrange [*:30]
set yrange [*:*]
#set key at  -5.2,1.095

set title '25 GeV^2 < Q^2 < 1000 GeV^2'

#set label 2 '4 GeV^2 < Q^2 < 5 GeV^2' font "Latin Modern Roman,28" at -8.9,1.7
#set label 28 'e^+ (27.6 GeV) p (920 GeV) → e^+ + X' font "Latin Modern Roman,28" at -1, 0.4 right
#set label 29 'NNPDF30\_nnlo\_as\_0118\_hera' font "Latin Modern Roman,28" at -1, 0.32 right
#set label 30 'μ_R = μ_F = Q' font "Latin Modern Roman,28" at -1, 0.24 right
#set label 31 '7-point variation' font "Latin Modern Roman,28" at -1, 0.16 right

set ylabel 'dσ/dp_{t,j_1} [pb/GeV]'
set xlabel 'p_{t,j_1} [GeV]'

ii=1

plot  lo i ii u (($1+$2)/2.):3:1:2:7:11 w @LOfill title 'LO',\
      lo i ii u (($1+$2)/2.):3:1:2 w @LOline not,\
      nlo i ii u (($1+$2)/2.):3:1:2:7:11 w @NLOfill title 'NLO',\
      nlo i ii u (($1+$2)/2.):3:1:2 w @NLOline not,\
      nnlo i ii u (($1+$2)/2.):3:1:2:7:11 w @NNLOfill title 'NNLO',\
      nnlo i ii u (($1+$2)/2.):3:1:2 w @NNLOline not
      
reset

set mxtics
set mytics
set grid
#set log y
#set format y "10^{%T}"
set xrange [*:*]
set yrange [*:*]
#set key at  -5.2,1.095

set title '25 GeV^2 < Q^2 < 1000 GeV^2'

#set label 2 '4 GeV^2 < Q^2 < 5 GeV^2' font "Latin Modern Roman,28" at -8.9,1.7
#set label 28 'e^+ (27.6 GeV) p (920 GeV) → e^+ + X' font "Latin Modern Roman,28" at -1, 0.4 right
#set label 29 'NNPDF30\_nnlo\_as\_0118\_hera' font "Latin Modern Roman,28" at -1, 0.32 right
#set label 30 'μ_R = μ_F = Q' font "Latin Modern Roman,28" at -1, 0.24 right
#set label 31 '7-point variation' font "Latin Modern Roman,28" at -1, 0.16 right

set ylabel 'dσ/dη_{t,j_1} [pb]'
set xlabel 'η_{t,j_1}'

ii=2

plot  lo i ii u (($1+$2)/2.):3:1:2:7:11 w @LOfill title 'LO',\
      lo i ii u (($1+$2)/2.):3:1:2 w @LOline not,\
      nlo i ii u (($1+$2)/2.):3:1:2:7:11 w @NLOfill title 'NLO',\
      nlo i ii u (($1+$2)/2.):3:1:2 w @NLOline not,\
      nnlo i ii u (($1+$2)/2.):3:1:2:7:11 w @NNLOfill title 'NNLO',\
      nnlo i ii u (($1+$2)/2.):3:1:2 w @NNLOline not
      
set output