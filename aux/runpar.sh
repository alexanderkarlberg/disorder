#!/bin/bash

# Executable location
EXEC=../disorder
MERGE=../aux/mergegrids
FLAGS="-Elep 27.6 -Epro 920 -positron -noincludeZ -Qmin 50.0 -pdf NNPDF30_nnlo_as_0118_hera -lo"
# Number of cores to run on
NCORES=8
FIRSTSEED=1 # The first seed
# Number of calls for each iteration setting up grids
NCALL1=1000000
#NCALL1=0
# Number of iterations used for grid generation. 3 should work. If
# more is needed increase ncall1 instead.
ITMX=3
# Number of calls for the production run
NCALL2=1000000
> Timings.txt

if test -f  "$EXEC"
then
    echo "Found $EXEC"
else
    echo "Did not find $EXEC"
    exit 0;
fi
if test -f  "$MERGE" 
then
    echo "Found $MERGE"
else
    echo "Did not find $MERGE"
    exit 0;
fi
EXEC="$EXEC $FLAGS"
# Function to take seed and convert into string
function char {
	case $1 in
	    [1-9])  echo 000$1 ;;
	    [1-9][0-9])  echo 00$1 ;;
	    [1-9][0-9][0-9])  echo 0$1 ;;
	esac
}

if [[ "$NCALL1" -gt 1  ]]
then
    echo 'Removing old grids and log files since ncall1 > 1'
    rm *grids* *log
    # Prepare grids using ITMX iterations
    for iteration in $(seq 1 $ITMX) 
    do
	(echo -n st1 xg$igrid ' ' ; date ) >> Timings.txt
	for core in $(seq $FIRSTSEED $NCORES)
	do
	    ch=`char $core`
	    $EXEC -ncall1 $NCALL1 -itmx1 1 -it1 $iteration -iseed $core -ncall2 0 > run-xg${iteration}-${ch}.log &
	done
	wait
	# Merge grids at this stage
	$MERGE grids-*.dat > xg${iteration}.log 
	# Replace grids
	for core in $(seq $FIRSTSEED $NCORES)
	do
	    ch=`char $core`
	    cp grids-${ch}.dat xg${iteration}-grids-${ch}.dat -v >> xg${iteration}.log 
	    cp grids-${ch}.top xg${iteration}-grids-${ch}.top -v >> xg${iteration}.log 
	    cp merged-grids.top grids-${ch}.top -v >> xg${iteration}.log 
	    cp merged-grids.dat grids-${ch}.dat -v >> xg${iteration}.log 
	done
	rm merged-grids.*
	wait
    done
    wait
fi

let ITMX+=1
# Production run
if [[ "$NCALL2" -gt 1  ]]
then
    (echo -n st2 ' ' ; date ) >> Timings.txt
    for core in $(seq $FIRSTSEED $NCORES)
    do
	ch=`char $core`
	rm hist*${ch}*top -v > run-prod-${ch}.log
	$EXEC -ncall1 0 -itmx1 0 -it1 $ITMX -iseed $core -ncall2 $NCALL2 -readingrid >> run-prod-${ch}.log &
    done
    wait
fi
(echo -n end ' ' ; date ) >> Timings.txt
exit 0;	 
