#!/bin/bash
# Build disorder in validation/build, the way the reference results were
# made, and run the tests through ctest:
#
#   ./validate_or_generate.sh validate   # unit tests + full validation matrix
#   ./validate_or_generate.sh quick      # the same without the slow runs
#   ./validate_or_generate.sh generate   # rerun the matrix and overwrite ref_runs/
#
# The validation matrix is defined in configurations.txt and each entry
# is compared with ref_runs/ by run_validation.py. After building, single
# tests can be rerun with e.g.
#   ctest --test-dir build -R validation_inclusive_cc_Q_10 --output-on-failure
#
# cmake needs to find hoppet, lhapdf and fastjet. If they are not in
# standard paths, pass the extra flags in EXTRA_CMAKEFLAGS, e.g.
#   EXTRA_CMAKEFLAGS="-DHOPPET_CONFIG=/path/to/hoppet-config" ./validate_or_generate.sh validate
CMAKEFLAGS="-DNEEDS_FASTJET=ON -DANALYSIS=exclusive_lab_frame_analysis.f ${EXTRA_CMAKEFLAGS}"

# Some colours for printout
RED='\033[0;31m'
GREEN='\033[0;32m'
PURPLE='\033[1;35m'
NC='\033[0m' # No Color

mode=$1
case $mode in
    validate) ctest_selection="" ;;
    quick)    ctest_selection="-LE slow" ;;
    generate) ;;
    *)
        echo "Need to specify validate, quick, or generate, like this"
        echo "./validate_or_generate.sh validate"
        exit 1
        ;;
esac

cd "$(dirname "$0")"
njobs=$(nproc 2>/dev/null || echo 4)

echo -e You have invoked the script to ${PURPLE}$mode${NC} the code
echo -e Building project in ${PURPLE}build${NC}
rm -rf build
cmake -S .. -B build $CMAKEFLAGS || exit 1
cmake --build build -j $njobs || exit 1

if [ "$mode" = "generate" ]; then
    prefixes=$(python3 run_validation.py --list)
    echo -e "Regenerating ${PURPLE}$(echo $prefixes | wc -w)${NC} reference runs in ref_runs/"
    echo "$prefixes" | xargs -P $njobs -I{} \
        python3 run_validation.py --generate --disorder build/disorder \
        --prefix {} --workdir build/generate/{} || exit 1
    echo -e ${PURPLE}DONE${NC} generating results
else
    if ctest --test-dir build -j $njobs --output-on-failure $ctest_selection; then
        echo -e All tests ${GREEN}PASSED${NC}
    else
        echo -e ERROR: At least one test ${RED}FAILED${NC}
        exit 1
    fi
fi

# Clean up
echo -e Cleaning up
rm -rf build
exit 0
