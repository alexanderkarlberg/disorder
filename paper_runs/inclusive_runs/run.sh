#!/bin/bash

FLAGS="-pdf MSHT20nnlo_as118 -Q 10 -x 0.01 -pdfuncert -scaleuncert -NC -noCC -no-analysis -Ehad 920 -Elep 27.6 -mz 91.1876 -mw 80.398 -includeZ -prefix NC_"
../../disorder $FLAGS -n3lo &
../../disorder $FLAGS -nnlo &
../../disorder $FLAGS -nlo &
../../disorder $FLAGS -lo &

FLAGS="-pdf MSHT20nnlo_as118 -Q 10 -x 0.01 -pdfuncert -scaleuncert -noNC -CC -no-analysis -Ehad 920 -Elep 27.6 -mz 91.1876 -mw 80.398 -prefix CC_"
../../disorder $FLAGS -n3lo &
../../disorder $FLAGS -nnlo &
../../disorder $FLAGS -nlo &
../../disorder $FLAGS -lo &


../../disorder -n3lo -pdf MSHT20an3lo_as118 -Q 10 -x 0.01 -pdfuncert -scaleuncert -NC -noCC -no-analysis -Ehad 920 -Elep 27.6 -mz 91.1876 -mw 80.398 -includeZ -prefix NC_an3lo_ &
../../disorder -n3lo -pdf MSHT20an3lo_as118 -Q 10 -x 0.01 -pdfuncert -scaleuncert -noNC -CC -no-analysis -Ehad 920 -Elep 27.6 -mz 91.1876 -mw 80.398 -prefix CC_an3lo_ &

