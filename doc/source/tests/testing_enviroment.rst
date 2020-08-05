-----------------------
The Testing Environment
-----------------------

The testing infrastructure in Enzo-e is comprised of unit tests (which test individual functions and functionality) and integration tests which test a more holistic portion of the codebase. In this documentationwe detail the layout of the testing infrastructure.

How to Run Tests
================

To run all the tests currently in the testing infrastructure run ``make test``. In order to run the a test separately write ``charmrun +p4 bin/enzo-p input/TestDirectory/test.in`` in the main enzo directory. If charmrun command is not in your path charmrun's path must also be included:  ``~/Charm/bin/charmrun +p4 bin/enzo-p input/TestDirectory/test.in``.

How Analyse Test Results
========================

``make test`` outputs first the unit test results in a count of passes and fails, and then a count of number of integration tests, the number of tests that start sucessfully and the number of tests that complete.

The unit test results count all times a unit test returns a fail, including when this happens during an integration test. 

If an integration test fails in start up this implies that there is something wrong with how the test is set up either as an error in the SConscript file that calls the test or that the test has an undefined parameter. In order to see what happened during the test, the test is recorded in the tests .unit file.

If a test fails while running this indicates that the feature being tested does not work.


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

See create new test
Create New Test:ref:`new_test`
