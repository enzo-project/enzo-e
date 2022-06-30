.. _new-test:

------------------------------
How to create a new ctest test
------------------------------

Introduction
============

This document provides instruction for creating tests for use with :ref:`ctest`. To make
an answer test, where results of a simulation can be run and compared from two versions
of ``enzo-e``, it is recommended to use :ref:`pytest`.

When adding a new function to the enzo-e project there is a possibility that it may affect the rest of the code, for this reason new code must be able to pass all existing tests as well as tests built for the new code. This means that integration tests are required. These tests are designed to stress the system. When making a new test these steps must be followed:

* Create input parameters
* Create a new test subdirectory
* Integrate the test into test infrastructure
* Running your new test

.. _new-test-simulation:

Create Input Parameters for the New Test.
=========================================


Create a new subdirectory for the new tests in the input directory. Create new ``.in`` file that contains input parameters for testing the new feature. This should be a simple problem such as a 2D implosion problem (See enzo-e/input/Hydro/test_implosion.in). More details on how to set input parameters see :ref:`parameter-file` . As an example of a parameter file see load-balance-4.in in the enzo-e/input/Balance directory. This file will be used later on to set the parameters of the test in the running the test section

Integrating the Test into Test Infrastructure
=============================================

To add a test to the ``ctest`` infrastructure, you have to add it to the
``test/CMakeLists.txt`` file.
Several different options exit depending on the type of the test and how it is being run.

Direct Enzo-E tests with input file
-----------------------------------

For tests that only execute the ``enzo-e`` binary with a given parameter file and
determine pass/fail by the exit code, you should use the ``setup_test_serial`` or
``setup_test_parallel`` function (in case the test should be run on more than one
compute unit).
Note that the number of compute units to run the parallel tests is set when the build
is configured via the ``PARALLEL_LAUNCHER_NPROC`` option.
Both functions take three arguments: (1) a short test name, (2) the output directory
(i.e., the directory within the ``build/test`` folder where potential output from the
test case is being stored), and (3) the input file.
A sample test case may look like ``setup_test_parallel(Bound-Periodic-2D BoundaryConditions/Periodic-2D  input/Boundary/boundary_periodic-2d.in)``.

Enzo-E test with output being analyzed/processed by Python script
-----------------------------------------------------------------

For tests that execute ``enzo-e`` but determine pass/fail by Python postprocessing,
e.g., using ``yt``.
To add this kind of test use either the ``setup_test_serial_python`` or
``setup_test_parallel_python`` function.
Again these functions have three mandatory arguments: (1) a short test name, (2) the
output directory, and (3) the Python script to be called.
Note that the Python script now includes the input file(s) that should be called
and the exit code of the script is used to determine pass/fail of the entire test.
Also note that these two function take additional arguments that are being passed to the
Python script.
A sample case is ``setup_test_serial_python(merge_stars_serial merge_stars/serial "input/merge_stars/run_merge_stars_test.py" "--prec=${PREC_STRING}")``.
We strongly encourage you to follow the overall structure used in the existing Python
scripts when setting up a new test case.

Unit tests
----------

Unit tests are self-contained test that are build as separate binaries.
In order to add these kind of test use the ``setup_test_unit`` function.
Again this function takes three arguments: (1) a short test name, (2) the output
directory, and (3) the binary name to be called.
A sample case is ``setup_test_unit(Array ArrayComponent/Array test_cello_array)``.

Note, several units test are currently commented out as there are linking issues.
For more information, see GitHub issue `#176 <https://github.com/enzo-project/enzo-e/issues/176>`_.

Running Your New Test
=====================

In order to run your new test you can call it directly, e.g., through the name
``ctest -R my-short-test-name``.

Note that in order for the test to show up in the build directory (if your build
has already been configured), you need to reconfigure the build so that the new
test case is being picked up.

To run all the tests use ``ctest`` in the build directory.
