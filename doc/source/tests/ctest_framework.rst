.. _ctest:

-------------------
The ctest Framework
-------------------

The ``ctest`` infrastructure in Enzo-e is comprised of unit tests (which test individual functions and functionality) and integration tests which test a more holistic portion of the codebase. In this documentation, we detail the layout of the ctest infrastructure. All input files are located in the input directory and the tests themselves are defined in ``test/CMakeLists.txt``

How to Run Tests
================

To run all the tests currently in the testing infrastructure run ``ctest`` in the build directory.

Basic parameters to control ``ctest`` are:

* ``ctest -N`` to see all available tests
* ``ctest -V`` to run all tests with verbose output
* ``ctest -R someregex`` to run all tests matching ``someregex``, e.g., ``ctest -R heat``
* ``ctest -L somelabel`` to run all test with label matching ``somelabel``. For ``Enzo-E`` we currently only set ``serial`` and ``parallel`` for serial and parallel (i.e., using multiple PEs) tests, respectively

*Note*, you may need to adjust the parallel launch configuration to your environment for the parallel test.
By default, ``charmrun`` with 4 PEs will be used to launch parallel tests.
To use ``srun`` instead (e.g., for a pure Charm++ MPI build) with 8 ranks you have to set
``-DPARALLEL_LAUNCHER=srun -DPARALLEL_LAUNCHER_NPROC_ARG="-n" -DPARALLEL_LAUNCHER_NPROC=8``
when configuring your build.

In order to run a test separately from the main ``ctest`` infrastructure write ``bin/enzo-e`` followed be the location of the test to be run in the main enzo directory. For example, ``bin/enzo-e input/Cosmology/method_cosmology-1.in`` runs the method_cosmology-1 test. Some tests are intended to run on multiple processors, for these tests the number of processors must be specified like so ``charmrun +p4 bin/enzo-e input/Cosmology/method_cosmology-8.in``, this runs the method_cosmology-8 test on four processors, a test designed to be run on multiple processors in parallel. If charmrun command is not in your path charmrun's path must also be included:  ``~/Charm/bin/charmrun +p4 bin/enzo-e input/Cosmology/method_cosmology-8.in``.

In order to exclude the ``"shu_collapse"`` and ``"bb_test"`` tests (as is done when running
the tests on CircleCI, execute the following command: ``ctest -E "(shu_collapse)|(bb_test)"``.


How to Analyse the Test Results
===============================

By default all tests will be run and the output is stored in the `test` directory of the build.

``ctest`` will automatically tell you which tests passed and which tests failed. For a more verbose output of the test, you can call ``ctest`` with ``--output-on-failure``.

If an integration test fails in start-up this implies that there is something wrong with how the test is set up either as an error in the ``test/CMakeLists.txt`` file that calls the test or that the test has an undefined parameter. If a test fails while running this indicates that the feature being tested does not work.

In order to see what happened during the test, you can look at the output directory of the test, which is located in a subdirectory (named after the test) in the ``test`` directory of the build directory, for example ``method_cosmology-1.in`` the test results are stored in ``test/MethodCosmology/Cosmology-8``. This directory also contains all outputs (image and data files).

