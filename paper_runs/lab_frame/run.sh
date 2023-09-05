#!/bin/bash

FLAGS="-p2b -pdf MSHT20nnlo_as118 -Qmin 5 -Qmax 31.62277660168379331998 -ymin 0.04 -ymax 0.95 -scaleuncert -Ehad 275 -Elep 18 -mz 91.1876 -mw 80.398 -prefix NC_EIC_"
../../disorder $FLAGS -lo   -ncall1 100000 -ncall2 100000000 &
../../disorder $FLAGS -nlo  -ncall1 100000 -ncall2 100000000 &
../../disorder $FLAGS -nnlo -ncall1 100000 -ncall2 10000000 &

