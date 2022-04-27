#!/bin/bash

# this script should be executed from the repository's root directory

# these tests should be integrated with the rest of the tests.

testPrefix="./input/merge_sinks"

echo "Beginning to run merge_sinks tests"

echo "Checking that the python installation roughly matches test requirements"
# One last set of crude checks to make sure that environment is not obviously
# unfit for running the tests
USING_PY3=$(python -c "exec(\"import sys\nprint(int(sys.version_info.major==3))\")")
if [[ "$USING_PY3" != "1" ]]; then
    echo "the 'python' command does not launch a version of python 3."
    echo "this is required for running merge_sinks tests."
    exit 1
fi

YT_INSTALLED=$(pip list | grep "^yt" | wc -l )
if [[ "$YT_INSTALLED" == "0" ]]; then
    echo "yt is not installed. It is required for testing merge_sinks tests."
    exit 1
fi

if [ ! -d ./test ]; then
    if [ -e ./test ]; then
        echo "A FILE called test exists that is not a directory"
        exit 1
    fi
    echo "There is no directory called ./test - Creating it now"
    mkdir test
fi



# Now to run the actual tests:

SERIAL_FAILED=0

echo ""
echo "Executing serial tests"
echo ""
echo "Beginning to run the serial stationary merge sinks test"
python $testPrefix/run_merge_sinks_test.py --launch_cmd ./bin/enzo-e --prec $CELLO_PREC --ics_type stationary > test/merge_sinks_stationary_serial.test-log 2>&1
ERR_CODE=$?				   
if [ $ERR_CODE -gt 0 ]; then
    ((SERIAL_FAILED++))
    echo "FAILED"
else
    echo "PASSED"
fi

echo "Beginning to run the serial drifting merge sinks test"
python $testPrefix/run_merge_sinks_test.py --launch_cmd ./bin/enzo-e --prec $CELLO_PREC --ics_type drift > test/merge_sinks_drift_serial.test-log 2>&1
ERR_CODE=$?				   
if [ $ERR_CODE -gt 0 ]; then
    ((SERIAL_FAILED++))
    echo "FAILED"
else
    echo "PASSED"
fi

N_SERIAL_TESTS=2
echo ""
echo "Finished running serial tests"
echo "${SERIAL_FAILED} set(s) of serial tests failed of ${N_SERIAL_TESTS}"
echo ""


PARALLEL_FAILED=0

echo ""
echo "Executing parallel tests"
echo ""
echo "Beginning to run the parallel stationary merge sinks test"
python $testPrefix/run_merge_sinks_test.py --launch_cmd "$CHARM_HOME/bin/charmrun +p4 ++local ./bin/enzo-e" --prec $CELLO_PREC --ics_type stationary > test/merge_sinks_stationary_parallel.test-log 2>&1
ERR_CODE=$?
if [ $ERR_CODE -gt 0 ]; then
    ((PARALLEL_FAILED++))
    echo "FAILED"
else
    echo "PASSED"
fi

echo "Beginning to run the parallel drifting merge sinks test"
python $testPrefix/run_merge_sinks_test.py --launch_cmd "$CHARM_HOME/bin/charmrun +p4 ++local ./bin/enzo-e" --prec $CELLO_PREC --ics_type drift > test/merge_sinks_drift_parallel.test-log 2>&1
ERR_CODE=$?
if [ $ERR_CODE -gt 0 ]; then
    ((PARALLEL_FAILED++))
    echo "FAILED"
else
    echo "PASSED"
fi

N_PARALLEL_TESTS=2
echo ""
echo "Finished running parallel tests"
echo "${PARALLEL_FAILED} set(s) of parallel tests failed of ${N_PARALLEL_TESTS}"
echo ""

SERIAL_PASSED="$(($N_SERIAL_TESTS-$SERIAL_FAILED))"
PARALLEL_PASSED="$(($N_PARALLEL_TESTS-$PARALLEL_FAILED))"

echo ""
echo "merge_sinks test summary:"
echo "$SERIAL_PASSED set(s) of serial tests PASSED of $N_SERIAL_TESTS"
if [[ N_PARALLEL_TESTS -eq 0 ]]; then
    echo "No parallel tests were run"
else
    echo "$PARALLEL_PASSED set(s) of parallel tests PASSED of $N_PARALLEL_TESTS"
fi

if [[ "$SERIAL_FAILED" -eq 0 && "$PARALLEL_FAILED" -eq 0 ]]; then
    echo "All executed sets of merge_sinks tests have passed"
    exit 0
else
    exit 1
fi
