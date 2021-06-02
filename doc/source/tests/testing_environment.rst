-----------------------
The Testing Environment
-----------------------

The testing infrastructure in Enzo-e is comprised of unit tests (which test individual functions and functionality) and integration tests which test a more holistic portion of the codebase. In this documentation, we detail the layout of the testing infrastructure. All tests are located in the input directory.

How to Run Tests
================

To run all the tests currently in the testing infrastructure run ``make test``. In order to run a test separately from the main script write ``bin/enzo-e`` followed be the location of the test to be run in the main enzo directory. For example, ``bin/enzo-e input/Cosmology/method_cosmology-1.in`` runs the method_cosmology-1 test. Some tests are intended to run on multiple processors, for these tests the number of processors must be specified like so ``charmrun +p4 bin/enzo-e input/Cosmology/method_cosmology-8.in``, this runs the method_cosmology-8 test on four processors, a test designed to be run on multiple processors in parallel. If charmrun command is not in your path charmrun's path must also be included:  ``~/Charm/bin/charmrun +p4 bin/enzo-e input/Cosmology/method_cosmology-8.in``.

How to Analyse the Test Results
===============================

``make test`` calls the build.sh script which outputs first the unit test results in a count of passes and fails, and then a count of number of integration tests, the number of tests that start successfully and the number of tests that complete.

The unit test results count all times a unit test returns a fail, including when this happens during an integration test. These results are stored in files in the test directory called pass.$CELLO_ARCH-$CELLO_PREC, incomplete.$CELLO_ARCH-$CELLO_PREC and fail.$CELLO_ARCH-$CELLO_PREC, where $CELLO_ARCH and $CELLO_PREC are the operating system and precision that Cello was set to during initial compiling of enzo code. ``build.sh`` automatically counts the number of passes incompletes and fails by counting the lines of these files.


If an integration test fails in start-up this implies that there is something wrong with how the test is set up either as an error in the SConscript file that calls the test or that the test has an undefined parameter. If a test fails while running this indicates that the feature being tested does not work.

In order to see what happened during the test, the test is recorded in the tests .unit file. Each test's .unit file is stored in a subdirectory of test that is relevant to the test, for example ``method_cosmology-1.in`` the test results are stored in ``test/MethodCosmology/test_method_cosmology-1.unit``. All outputs of tests (images and data files) are stored in the within the subdirectory in their own directory. 


What Tests are Currently Included
=================================

Currently the Enzo-e testing infrastructure tests:

.. toctree::
   adapt
   balance
   boundary
   checkpoint
   collapse
   cosmology
   gravity
   heat
   helloworld
   hierarchy
   hydro
   initialmusic
   methods
   output
   particle
   performance
   ppm
   ppml
   sedov
  
How to Add Your Own Test
========================


:ref:`new-test`.
