#!/bin/bash

# this script should be executed from the repository's root directory

# these tests should be integrated with the rest of the tests.

serialTests=("run_shock_tube_test" "run_MHD_linear_wave_test"
             "run_HD_linear_wave_test" "run_passive_advect_sound_test"
             "run_dual_energy_cloud_test")
# The parallelTests typically just take too long to run serially
parallelTests=("run_dual_energy_shock_tube_test")
testPrefix="./input/vlct"

if [ "$CELLO_PREC" == "double" ]; then
    echo "Beginning to run vlct_tests"
else
    echo "vlct_tests should only be run when compiled with double Precision"
    exit 1
fi



echo ""
echo "Determining number of processes to use for parallel tests"
if [ -z "$VL_PARALLEL_TEST_IP_CHARM" ]; then
    N_PE_PARALLEL=4
else
    N_PE_PARALLEL=$VL_PARALLEL_TEST_IP_CHARM
fi

if [[ "$N_PE_PARALLEL" -lt 1 ]]; then
    echo "VL_PARALLEL_TEST_IP_CHARM is less than 1. Parallel tests will be skipped."
else
    echo "Using $N_PE_PARALLEL processes"
fi



if [[ ! -z "$PYTHON_VL_TEST_PREP" ]]; then
    echo ""
    echo "Preparing python environment for test"
    eval "$PYTHON_VL_TEST_PREP"
fi

echo "Checking that the python installation roughly matches test requirements"
# One last set of crude checks to make sure that environment is not obviously
# unfit for running the tests
USING_PY3=$(python -c "exec(\"import sys\nprint(int(sys.version_info.major==3))\")")
if [[ "$USING_PY3" != "1" ]]; then
    echo "the 'python' command does not launch a version of python 3."
    echo "this is required for running VLCT tests."
    exit 1
fi

YT_INSTALLED=$(pip list | grep "^yt" | wc -l )
if [[ "$YT_INSTALLED" != "1" ]]; then
    echo "yt is not installed. It is required for testing VLCT tests."
    exit 1
fi


# Now to run the actual tests:

SERIAL_FAILED=0

echo ""
echo "Executing serial tests"
for i in ${!serialTests[@]}; do
    test_path="${testPrefix}/${serialTests[$i]}.py"
    fname_test="test/${serialTests[$i]}.test-log"
    echo ""
    echo "Beginning to run ${test_path}"
    python ${test_path} > $fname_test 2>&1
    ERR_CODE=$?
    if [ $ERR_CODE -gt 0 ]; then
        ((SERIAL_FAILED++))
        echo "FAILED"
    else
        echo "PASSED"
    fi
done
N_SERIAL_TESTS="${#serialTests[@]}"
echo ""
echo "Finished running serial tests"
echo "${SERIAL_FAILED} set(s) of serial tests failed of ${N_SERIAL_TESTS}"
echo ""

PARALLEL_FAILED=0
if [[ "$N_PE_PARALLEL" -lt 1 ]]; then
    echo "SKIPPING PARALLEL TESTS"
    N_PARALLEL_TESTS=0
else
    echo ""
    echo "Starting parallel tests"
    for i in ${!parallelTests[@]}; do
        test_path="${testPrefix}/${parallelTests[$i]}.py"
        fname_test="test/${parallelTests[$i]}.test-log"
        echo ""
        echo "Beginning to run ${test_path}"
        python $test_path $N_PE_PARALLEL > $fname_test 2>&1
        ERR_CODE=$?
        if [ $ERR_CODE -gt 0 ]; then
            ((PARALLEL_FAILED++))
            echo "FAILED"
        else
            echo "PASSED"
        fi
    done
    N_PARALLEL_TESTS="${#parallelTests[@]}"
    echo ""
    echo "Finished running parallel tests"
    echo "${PARALLEL_FAILED} set(s) of parallel tests failed of $N_PARALLEL_TESTS"
    echo ""
fi

SERIAL_PASSED="$(($N_SERIAL_TESTS-$SERIAL_FAILED))"
PARALLEL_PASSED="$(($N_PARALLEL_TESTS-$PARALLEL_FAILED))"

echo ""
echo "VLCT test summary:"
echo "$SERIAL_PASSED set(s) of serial tests PASSED of $N_SERIAL_TESTS"
if [[ N_PARALLEL_TESTS -eq 0 ]]; then
    echo "No parallel tests were run"
else
    echo "$PARALLEL_PASSED set(s) of parallel tests PASSED of $N_PARALLEL_TESTS"
fi

if [[ "$SERIAL_FAILED" -eq 0 && "$PARALLEL_FAILED" -eq 0 ]]; then
    echo "All executed sets of VLCT tests have passed"
    exit 0
else
    exit 1
fi
