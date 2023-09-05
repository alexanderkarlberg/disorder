set terminal pdf enhanced font "Latin Modern Roman,36" size 30cm,40cm
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
NNLOfill= "boxxyerrorbars lc rgb '#e51e10' lw 2.0 fs transparent solid 0.2"


lo='<paste NC_EIC_disorder_lo_central.dat NC_EIC_disorder_lo_max.dat NC_EIC_disorder_lo_min.dat         NC_EIC_disorder_nlo_central.dat '
nlo='<paste NC_EIC_disorder_nlo_central.dat NC_EIC_disorder_nlo_max.dat NC_EIC_disorder_nlo_min.dat     NC_EIC_disorder_nlo_central.dat '
nnlo='<paste NC_EIC_disorder_nnlo_central.dat NC_EIC_disorder_nnlo_max.dat NC_EIC_disorder_nnlo_min.dat NC_EIC_disorder_nlo_central.dat '


reset

set mxtics
set mytics
set grid
set xtics format ""
set log y
set format y "10^{%T}"
set xrange [*:30]
set yrange [*:*]
#set key at  -5.2,1.095

set title '25 GeV^2 < Q^2 < 1000 GeV^2,  0.04 < y < 0.95'

#set label 2 '4 GeV^2 < Q^2 < 5 GeV^2' font "Latin Modern Roman,28" at -8.9,1.7
set label 28 'e^- (18 GeV) p (275 GeV) → e^- + X' font "Latin Modern Roman,28" at 6, 1.8/1000. left
set label 29 'MSHT20nnlo\_as118' font "Latin Modern Roman,28" at 6, 1.0/1000. left
set label 30 'μ_R = μ_F = Q' font "Latin Modern Roman,28" at 6, 0.5555555/1000. left
set label 31 '7-point variation' font "Latin Modern Roman,28" at 6, .30864197530864197530/1000. left
set label 38 '|η_j| < 3, anti-k_t, R=0.8' font "Latin Modern Roman,28" at 6, .17146776406035665294/1000. left

set ylabel 'dσ/dp_{t,j} [mb/GeV]'
set xlabel ''
set multiplot
set origin 0.0,0.3
set size 1.0,1.0

set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.92
set bmargin at screen 0.47


ii=1

plot  lo i ii u   (($1+$2)/2.):($3/1000.):1:2:($7/1000.):($11/1000.) w @LOfill title 'LO',\
      lo i ii u   (($1+$2)/2.):($3/1000.):1:2 w @LOline not,\
      nlo i ii u  (($1+$2)/2.):($3/1000.):1:2:($7/1000.):($11/1000.) w @NLOfill title 'NLO',\
      nlo i ii u  (($1+$2)/2.):($3/1000.):1:2 w @NLOline not,\
      nnlo i ii u (($1+$2)/2.):($3/1000.):1:2 w @NNLOline not,\
      nnlo i ii u (($1+$2)/2.):($3/1000.):1:2:($7/1000.):($11/1000.) w @NNLOfill title 'NNLO',\

      
set tmargin at screen 0.47
set bmargin at screen 0.10
set nologscale y
set title ''
set format x
set format y
set ytics 0.0,0.5,1.99
set mytics 5
unset label
set xlabel 'p_{t,j [GeV]'
set ylabel 'Ratio to NLO'
set yrange [0.0:2.0]

plot  lo i ii u   (($1+$2)/2.):($3/$15):1:2:($7/$15):($11/$15) w @LOfill not,\
      lo i ii u   (($1+$2)/2.):($3/$15):1:2 w @LOline not,\
      nlo i ii u  (($1+$2)/2.):($3/$15):1:2:($7/$15):($11/$15) w @NLOfill not,\
      1 lc rgb '#009e73'  not,\
      nnlo i ii u (($1+$2)/2.):($3/$15):1:2 w @NNLOline not,\
      nnlo i ii u (($1+$2)/2.):($3/$15):1:2:($7/$15):($11/$15) w @NNLOfill not',\

unset multiplot
reset

set mxtics
set mytics
set grid
set xtics format ""
#set log y
#set format y "10^{%T}"
set xrange [*:*]
set yrange [0:*]
#set key at  -5.2,1.095

set title '25 GeV^2 < Q^2 < 1000 GeV^2,  0.04 < y < 0.95'

set label 28 'e^- (18 GeV) p (275 GeV) → e^- + X' font "Latin Modern Roman,28" at -0.6, 1.6 left
set label 29 'MSHT20nnlo\_as118' font "Latin Modern Roman,28" at -0.6, 1.3 left
set label 30 'μ_R = μ_F = Q' font "Latin Modern Roman,28" at -0.6, 1.0 left
set label 31 '7-point variation' font "Latin Modern Roman,28" at -0.6, 0.7 left
set label 38 '|η_j| < 3, anti-k_t, R=0.8' font "Latin Modern Roman,28" at -0.6, 0.4 left


set ylabel 'dσ/dη_{j} [mb]'

set multiplot
set origin 0.0,0.3
set size 1.0,1.0

set lmargin at screen 0.15
set rmargin at screen 0.95
set tmargin at screen 0.92
set bmargin at screen 0.47

ii=2

plot  lo i ii u   (($1+$2)/2.):($3/1000.):1:2:($7/1000.):($11/1000.) w @LOfill title 'LO',\
      lo i ii u   (($1+$2)/2.):($3/1000.):1:2 w @LOline not,\
      nlo i ii u  (($1+$2)/2.):($3/1000.):1:2:($7/1000.):($11/1000.) w @NLOfill title 'NLO',\
      nlo i ii u  (($1+$2)/2.):($3/1000.):1:2 w @NLOline not,\
      nnlo i ii u (($1+$2)/2.):($3/1000.):1:2 w @NNLOline not,\
      nnlo i ii u (($1+$2)/2.):($3/1000.):1:2:($7/1000.):($11/1000.) w @NNLOfill title 'NNLO',\

set tmargin at screen 0.47
set bmargin at screen 0.10
set nologscale y
set title ''
set format x
set format y
set ytics 0.0,0.5,1.99
set mytics 5
unset label
set xlabel 'η_{j}'
set ylabel 'Ratio to NLO'
set yrange [0.0:2.0]
plot  lo i ii u   (($1+$2)/2.):($3/$15):1:2:($7/$15):($11/$15) w @LOfill not,\
      lo i ii u   (($1+$2)/2.):($3/$15):1:2 w @LOline not,\
      nlo i ii u  (($1+$2)/2.):($3/$15):1:2:($7/$15):($11/$15) w @NLOfill not,\
      1 lc rgb '#009e73'  not,\
      nnlo i ii u (($1+$2)/2.):($3/$15):1:2 w @NNLOline not,\
      nnlo i ii u (($1+$2)/2.):($3/$15):1:2:($7/$15):($11/$15) w @NNLOfill not',\

unset multiplot 

set output