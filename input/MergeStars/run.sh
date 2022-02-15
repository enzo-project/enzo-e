#!/bin/bash

# this script should be executed from the repository's root directory

# these tests should be integrated with the rest of the tests.

testPrefix="./input/MergeStars"

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
if [[ "$YT_INSTALLED" == "0" ]]; then
    echo "yt is not installed. It is required for testing VLCT tests."
    exit 1
fi


# Now to run the actual tests:
echo "Generating ICs for MergeStars Test"

echo $testPrefix
python $testPrefix/ics.py -r 4.0e16 -m 2.0e36 -c 1.5e17 1.5e17 1.5e17 -d 2.0e8 2.0e8 2.0e8 -i 1.0e8 -n 1000

# Run serial test
echo "Running Serial MergeStars Test"
bin/enzo-e $testPrefix/merge_stars_test_serial.in > test/MethodMergeStars/merge_stars_test_serial.out 2>&1

mpirun -np $N_PE_PARALLEL python $testPrefix/images.py -i Dir_Merge_Stars_Serial -o image_serial > images_serial.out 2>&1

mkdir -p test/MethodMergeStars/run_serial
mv image_serial* test/MethodMergeStars/run_serial/

mpirun -np $N_PE_PARALLEL python $testPrefix/mass_momentum_conservation.py -i Dir_Merge_Stars_Serial -o mmc_serial.png > mmc_serial.out 2>&1

ERR_CODE=$?
if [ $ERR_CODE -gt 0 ]; then
    SERIAL_FAILED=1
    echo "SERIAL MERGE STARS TEST FAILED"
else
    SERIAL_FAILED=0
    echo "SERIAL MERGE STARS TEST PASSED"
fi

mv mmc_serial.png test/MethodMergeStars/


# Run parallel test
echo "Running Parallel MergeStars Test"
charmrun +p$N_PE_PARALLEL ++local bin/enzo-e $testPrefix/merge_stars_test_parallel.in > test/MethodMergeStars/merge_stars_test_parallel.out 2>&1

mpirun -np $N_PE_PARALLEL python $testPrefix/images.py -i Dir_Merge_Stars_Parallel -o image_parallel > images_parallel.out 2>&1

mkdir -p test/MethodMergeStars/run_parallel
mv image_parallel* test/MethodMergeStars/run_parallel/

mpirun -np $N_PE_PARALLEL python $testPrefix/mass_momentum_conservation.py -i Dir_Merge_Stars_Parallel -o mmc_parallel.png > mmc_parallel.out 2>&1

ERR_CODE=$?
if [ $ERR_CODE -gt 0 ]; then
    PARALLEL_FAILED=1
    echo "PARALLEL MERGE STARS TEST FAILED"
else
    PARALLEL_FAILED=0
    echo "PARALLEL MERGE STARS TEST PASSED"
fi

mv mmc_parallel.png test/MethodMergeStars/

if [[ "$SERIAL_FAILED" -eq 0 && "$PARALLEL_FAILED" -eq 0 ]]; then
    echo "Both MergeStars tests have passed"
    exit 0
else
    exit 1
fi