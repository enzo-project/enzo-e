#!/bin/bash

# this script should be executed from the repository's root directory

# these tests should be integrated with the rest of the tests.

testPrefix="./input/MergeStars"

mkdir -p test/MethodMergeStars

echo "Checking that the python installation roughly matches test requirements"
# One last set of crude checks to make sure that environment is not obviously
# unfit for running the tests
USING_PY3=$(python -c "exec(\"import sys\nprint(int(sys.version_info.major==3))\")")
if [[ "$USING_PY3" != "1" ]]; then
    echo "the 'python' command does not launch a version of python 3."
    echo "this is required for running MergeStars tests."
    exit 1
fi

YT_INSTALLED=$(pip list | grep "^yt" | wc -l )
if [[ "$YT_INSTALLED" == "0" ]]; then
    echo "yt is not installed. It is required for running MergeStars tests"
    exit 1
fi

MPI4PY_INSTALLED=$(pip list | grep "^mpi4py" | wc -l )
if [[ "$MPI4PY_INSTALLED" == "0" ]]; then
    echo "mpi4py is not installed. It is required for running MergeStars tests"
    exit 1
fi


# Now to run the actual tests:
echo "Generating ICs for MergeStars Test"

python $testPrefix/ics.py -r 4.0e16 -m 2.0e36 -c 1.5e17 1.5e17 1.5e17 -d 2.0e8 2.0e8 2.0e8 -i 1.0e8 -n 1000

# Run serial test
echo "Running Serial MergeStars Test"
bin/enzo-e $testPrefix/merge_stars_test_serial.in > test/MethodMergeStars/merge_stars_test_serial.out 2>&1

cat test/MethodMergeStars/merge_stars_test_serial.out

mpirun -np 4 python $testPrefix/images.py -i Dir_Merge_Stars_Serial -o image_serial > test/MethodMergeStars/images_serial.out 2>&1

cat test/MethodMergeStars/images_serial.out

mkdir -p test/MethodMergeStars/run_serial
mv image_serial* test/MethodMergeStars/run_serial/

mpirun -np 4 python $testPrefix/mass_momentum_conservation.py -i Dir_Merge_Stars_Serial -o mmc_serial.png > test/MethodMergeStars/mmc_serial.out 2>&1

cat test/MethodMergeStars/mmc_serial.out

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
charmrun +p4 ++local bin/enzo-e $testPrefix/merge_stars_test_parallel.in > test/MethodMergeStars/merge_stars_test_parallel.out 2>&1

mpirun -np 4 python $testPrefix/images.py -i Dir_Merge_Stars_Parallel -o image_parallel > test/MethodMergeStars/images_parallel.out 2>&1

mkdir -p test/MethodMergeStars/run_parallel
mv image_parallel* test/MethodMergeStars/run_parallel/

mpirun -np 4 python $testPrefix/mass_momentum_conservation.py -i Dir_Merge_Stars_Parallel -o mmc_parallel.png > test/MethodMergeStars/mmc_parallel.out 2>&1

ERR_CODE=$?
if [ $ERR_CODE -gt 0 ]; then
    PARALLEL_FAILED=1
    echo "PARALLEL MERGE STARS TEST FAILED"
else
    PARALLEL_FAILED=0
    echo "PARALLEL MERGE STARS TEST PASSED"
fi

mv mmc_parallel.png test/MethodMergeStars/

rm -r Dir*

if [[ "$SERIAL_FAILED" -eq 0 && "$PARALLEL_FAILED" -eq 0 ]]; then
    echo "Both MergeStars tests have passed"
    exit 0
else
    exit 1
fi
