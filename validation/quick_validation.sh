#!/bin/bash
# To validate the program run
# ./validate_or_generate.sh validate
#
# To generate new validation runs
# ./validate_or_generate.sh generate
#
# The script is not very flexible as it needs cmake to work with
# standard paths to hoppet, lhapdf and fastjet. If that is not the
# case, then the user should modify the below
CMAKEFLAGS=" -DNEEDS_FASTJET=ON -DANALYSIS=exclusive_lab_frame_analysis.f "

# Some colours for printout
RED='\033[0;31m'
GREEN='\033[0;32m'
PURPLE='\033[1;35m'
NC='\033[0m' # No Color

# Clean-up any semaphores and old builds
rm -rf ~/.parallel/semaphores/* build

prefix="
inclusive_nc_Q_10_x_0.01_ 
inclusive_nc_includeZ_Q_10_x_0.01_ 
inclusive_cc_Q_10_x_0.01_ 
inclusive_nc_Q_10_ 
inclusive_nc_includeZ_Q_10_ 
inclusive_cc_Q_10_ 
inclusive_nc_Qmin_1_x_0.01_ 
inclusive_nc_includeZ_Qmin_1_x_0.01_ 
inclusive_cc_Qmin_1_x_0.01_ 
inclusive_nc_Q_10_y_0.01_ 
inclusive_nc_includeZ_Q_10_y_0.01_ 
inclusive_cc_Q_10_y_0.01_ 
inclusive_nc_Q_10_x_0.01_neutrino_ 
inclusive_nc_Q_10_x_0.01_neutrino_positron_ 
inclusive_cc_Q_10_x_0.01_neutrino_ 
inclusive_cc_Q_10_x_0.01_neutrino_positron_ 
p2b_nc_Q_10_x_0.01_
p2b_nc_Q_10_x_0.01_MSHT20an3lo_as118_
"
prefixarray=($prefix)

cmdline=(
    # Some inclusive runs
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -toyQ0\ 2.0\ -Q\ 10.0\ -x\ 0.01\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -includeZ\ -positron\ -toyQ0\ 2.0\ -Q\ 10.0\ -x\ 0.01\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -CC\ -toyQ0\ 2.0\ -Q\ 10.0\ -x\ 0.01\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -toyQ0\ 2.0\ -Q\ 10.0\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -includeZ\ -toyQ0\ 2.0\ -Q\ 10.0\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -CC\ -toyQ0\ 2.0\ -Q\ 10.0\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -toyQ0\ 2.0\ -Qmin\ 1.0\ -x\ 0.01\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -includeZ\ -toyQ0\ 2.0\ -Qmin\ 1.0\ -x\ 0.01\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -CC\ -toyQ0\ 2.0\ -Qmin\ 1.0\ -x\ 0.01\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -toyQ0\ 2.0\ -Q\ 10.0\ -y\ 0.01\ -scaleuncert\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -includeZ\ -toyQ0\ 2.0\ -Q\ 10.0\ -y\ 0.01\ -scaleuncert\  
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -CC\ -toyQ0\ 2.0\ -Q\ 10.0\ -y\ 0.01\ -scaleuncert\
    # Some neutrino runs
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -toyQ0\ 2.0\ -Q\ 10.0\ -x\ 0.01\ -scaleuncert\ -neutrino\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -toyQ0\ 2.0\ -Q\ 10.0\ -x\ 0.01\ -scaleuncert\ -neutrino\ -positron\
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -noNC\ -CC\ -toyQ0\ 2.0\ -Q\ 10.0\ -x\ 0.01\ -scaleuncert\ -neutrino\ 
    -n3lo\ -yorder\ 5\ -lnlnQorder\ 4\ -noNC\ -CC\ -toyQ0\ 2.0\ -Q\ 10.0\ -x\ 0.01\ -scaleuncert\ -neutrino\ -positron\
    #Some p2b runs
    -nnlo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -toyQ0\ 2.0\ -Q\ 10.0\ -x\ 0.01\ -scaleuncert\ -p2b\ 
    -nnlo\ -yorder\ 5\ -lnlnQorder\ 4\ -NC\ -pdf\ MSHT20an3lo_as118\ -Q\ 10.0\ -x\ 0.01\ -scaleuncert\ -p2b\ 
)

if [ "${#prefixarray[@]}" -ne "${#cmdline[@]}" ]; then
    echo "Arrays are not fo the same size " ${#prefixarray[@]}  ${#cmdline[@]}
    exit 1
fi

numJobs=${#prefixarray[@]}

dir="test_runs"

rm -rf $dir
mkdir $dir

echo -e You have invoked the script to ${PURPLE}quickly validate${NC} the code

# Create build directory and compile
echo -e Building project in ${PURPLE}build${NC}
mkdir build 
cd build 
cmake ../.. $CMAKEFLAGS #> build.log
make -j #>> build.log
# Uncomment for CI debug
#ldd disorder
#ls /usr/local/lib/

# Move to directory containing reference results
cd ../$dir

echo -e "Starting ${PURPLE}$numJobs${NC} jobs (in parallel if possible)"
iJob=1
for i in $(seq 0 $((numJobs-1)))
do
    echo -e Running job number ${iJob}: ${PURPLE}../build/disorder ${cmdline[$i]} -prefix ${prefixarray[$i]}${NC}
#    sem -j 50% ../build/disorder ${cmdline[$i]} -prefix ${prefixarray[$i]} &> ${prefixarray[$i]%_}.log
    sem --use-cores-instead-of-threads -j +0 ../build/disorder ${cmdline[$i]} -prefix ${prefixarray[$i]} &> ${prefixarray[$i]%_}.log
    # Uncomment for CI debug
#    sem -j 50% ../build/disorder ${cmdline[$i]} -prefix ${prefixarray[$i]} 2>&1 | tee ${prefixarray[$i]%_}.log
    ((iJob++))
done
sem --wait 

echo -e ${PURPLE}DONE${NC} generating results

# If we are generating then nothing more to do. If we are validating then now is the time!
for file_w_path in ../ref_runs_quick/*
do
    file=${file_w_path#../ref_runs_quick/}
    echo -e Comparing output of ${PURPLE}$file${NC}
    # First remove some useless lines
    grep -v "TOTAL TIME" $file_w_path | grep -v "Stamped by" | grep -v "FastJet" > ${file}.ref
    grep -v "TOTAL TIME" $file | grep -v "Stamped by" | grep -v "FastJet" > ${file}.new
    diff  ${file}.ref ${file}.new > ${file}.diff
    checkwc=`cat ${file}.diff| wc -l `
    if [ $checkwc == "0" ]; then
	echo -e "Comparison                                                                                           ${GREEN}PASSED${NC}"
    else
	echo -e "Comparison                                                                                           ${RED}FAILED${NC}"
	cat ${file}.diff
	failed="true"
    fi
done
if [ -z $failed ]; then
    echo -e All tests ${GREEN}PASSED${NC}
else
    echo -e ERROR: At least one test ${RED}FAILED${NC}
    exit 1
fi

# Clean up
echo -e Cleaning up

rm *grids* 
cd ..
rm -rf build test_runs

exit 0
