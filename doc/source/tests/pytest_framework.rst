.. _pytest:

--------------------
The pytest Framework
--------------------

The ``pytest`` framework is primarily used for running answer tests, where a
simulation is run with two versions of ``enzo-e`` and their results are compared.
This is useful for testing problems with no analytical solution or generally
verifying that results from commonly run simulations don't drift.

It is also useful in for testing problems that do have analytic solutions (the answer test might quantify how close a simulation result is to the analytic expected solution).
While such tests do exist in the ctest-framework, they often involve more boiler-plate code.

`pytest <https://docs.pytest.org/>`__ is a Python-based framework for detecting
and running a series of tests within a source code repository. When running
``pytest``, the user can provide a directory in which ``pytest`` will look for
files named `test_*.py` and run all functions within those files whose names start
with "test". ``pytest`` will run all tests and present a summary of which ones
passed and failed. All functions that run without producing an error will be marked
as passed.

Installation
============

``pytest`` can be installed with ``pip`` or ``conda``.

.. code-block:: bash

   $ pip install pytest

.. code-block:: bash

   $ conda install pytest

Answer Testing
==============

Within ``enzo-e``, we make use of the `TestCase
<https://docs.python.org/3/library/unittest.html#unittest.TestCase>`_ class to
define a general ``EnzoETest`` class that will run a given simulation within a
temporary directory and delete that directory once finished. This class and
other useful answer testing functionality are located in the source in
`test/answer_tests/answer_testing.py`. All answer tests are located in the
other files within the `test/answer_tests` directory.

Some other functionality, that may be reused in other unrelated scripts provided in the Enzo-E repository, are provided in the ``test_utils`` subdirectory.

Running the Answer Test Suite
-----------------------------

The answer test suite is run in two stages. First, test answers must be generated
from a version of the code known to function correctly. A git tag associated with
the main repository marks a changeset for which the code is believed to produce
good results. This tag is named ``gold-standard-#``. To pull tags from the main
repository and see which tags exist, do the following:

.. code-block:: bash

   $ git fetch origin --tags
   $ git tag

To generate test answers, use the highest numbered gold standard tag.

Configuring the Answer Test Suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The behavior of the test can be configured by passing command line arguments to ``pytest`` or by setting environment variables (or by mixing both).

When invoking ``pytest``, the command line flags discussed here should be passed **after** the path to the `test/answer_tests` directory has been provided.
For the sake of example (the meaning of flags are explained below), one might invoke:

.. code-block:: bash

   $ pytest test/answer_tests  \
            --build-dir ./build \
            --answer-store

The following table lists command line flags, and where applicable, the environment variables that they are interchangable with.
In cases where both are set, the command line argument is given precedence.

.. list-table:: Configuring pytest behavior
   :widths: 10 10 30
   :header-rows: 1

   * - flag
     - env var
     - description
   * - ``--build-dir``
     - N/A
     - points to the build-directory where the target enzo-e binary was built (that binary has the path: BUILD_DIR/bin/enzo-e).
       The path to the charmrun launcher will be inferred from the `BUILD_DIR/CMakeCache.txt` file, but can be overwritten by the ``--charm`` flag or the ``CHARM_PATH`` environment variable.
       This precedence was chosen in case a user causes a change to relevant cached build-variables, but have not rebuilt Enzo-E (i.e. `CMakeCache.txt` may not be valid for the binary).
       When this flag isn't specified, the test infrastructure searches for the enzo-e binary at ENZOE_ROOT/build/bin/enzo-e, but doesn't try to infer charmrun's location from `CMakeCache.txt`.
   * - ``--local-dir``
     - ``TEST_RESULTS_DIR``
     - points to a directory in which answers will be stored/loaded
   * - ``--charm``
     - ``CHARM_PATH``
     - points to the directory in which ``charmrun`` is located
   * - ``--answer-store``
     - ``GENERATE_TEST_RESULTS``
     - When the command line flag is specified, test results are generated. Otherwise, results are compared against existing results (unless the environment variable is specified).
       The environment variable can be be set to ``"true"`` to generate test results or ``"false"`` to compare with existing results.
   * - ``--grackle-input-data-dir``
     - ``GRACKLE_INPUT_DATA_DIR``
     - points to the directory where ``Grackle`` input files are installed.
       If not specified, then all tests involving ``Grackle`` will be skipped.

Earlier versions of the tests also required the ``"USE_DOUBLE"`` environment variable to be set to ``"true"`` or ``"false"`` to indicate whether the code had been compiled in double or single precision.

.. code-block:: bash

   $ export TEST_RESULTS_DIR=~/enzoe_tests
   $ export CHARM_PATH=~/local/charm-v7.0.0/bin

Generating Test Answers
^^^^^^^^^^^^^^^^^^^^^^^

First, check out the highest numbered gold standard tag and compile ``enzo-e``.

.. code-block:: bash

   # in the future, you will need to subsitute 004 for a higher number
   $ git checkout gold-standard-004
   $ ...compile enzo-e

Then, run the test suite by calling ``pytest`` with the answer test directory (make sure to configure behavior correctly with command-line arguments or environment variables).
In the following snippet, we assume you are currently at the root of the Enzo-E repository and that you will replace ``<build-dir>`` with the directory where you build enzo-e (this is commonly ``./build``)

.. code-block:: bash

   $ pytest test/answer_tests --local-dir=~/enzoe_tests --build-dir=<build-dir> --answer-store
   ========================== test session starts ===========================
   platform linux -- Python 3.9.13, pytest-7.1.2, pluggy-1.0.0
   rootdir: /home/circleci/enzo-e
   collected 1 item

   test/answer_tests/test_vlct.py .                                   [100%]

   =========================== 1 passed in 13.26s ===========================

Assuming there are no errors, this will run the simulations associated with the
tests, perform the analysis required to produce the answers, save the answers to
files, and report that all tests have passed.

Comparing Test Answers
^^^^^^^^^^^^^^^^^^^^^^

Once test answers have been generated, the above steps need not be repeated until
the gold standard tag has been updated. Now, any later version of the code can be
run with the test suite to check for problems. To configure the test suite to compare with existing answers, omit the ``--answer-store`` flag and ensure that the ``GENERATE_TEST_RESULTS`` variable is either unset or set to ``"false"``.

.. code-block:: bash

   $ git checkout main
   $ ...compile enzo-e
   $ pytest test/answer_tests --local-dir=~/enzoe_tests --build-dir=<build-dir>

Helpful Tips
^^^^^^^^^^^^

By default, most output printed by ``enzo-e`` or the test scripts will be swallowed
by ``pytest``. When tests fail, the Python traceback may be shown, but not much
else. There are various flags to increase the verbosity of ``pytest``, but the
``-s`` flag will show all output, including from the simulation itself. The
``enzo-e`` answer test suite will also print out the values of all configuration
variables when this flag is given.

.. code-block:: bash

   $ pytest -s test/answer_tests # other args...

When debugging an issue it's sometimes helpful to force pytest to run a subset of tests.
This can be accomplished with the ``-k`` flag.
For example, to only run a subset of tests with ``"grackle"`` in the test name, one might execute

.. code-block:: bash

   $ pytest test/answer_tests -k "grackle" # other args...

When investigating a failing test or prototyping a brand-new test, it can sometimes be helpful to run the tests against multiple versions of enzo-e.
Rather than rebuilding Enzo-E each time you want to do that, you can instead build the different versions of Enzo-E in separate build-directories, and direct ``pytest`` to use the different builds with the ``--build-dir`` flag.

Creating New Answer Tests
-------------------------

This section follows the example of ``TestHLLCCloud`` in
`test/answer_tests/test_vlct.py`. Answer tests can be created by making a new Python
file in the `test/answer_tests` directory with a name starting with 'test\_' or by
adding to an existing file if the test falls within the theme given by its name. If
your test requires configuring a new simulation parameter file, see
:ref:`new-test-simulation` for information on setting that up.

The answer testing framework exists in `test/answer_tests/answer_testing.py`. New
test files created in the same directory can directly import from this file.

Creating a New Test Class
^^^^^^^^^^^^^^^^^^^^^^^^^

To make a new test, one must create a new Python class that subclasses the
``EnzoETest`` class. Three attributes must be defined within the class:

 * parameter_file: the relative path to the simulation parameter file from within
   the input directory.
 * max_runtime: the maximum runtime of the simulation in seconds. The simulation
   will be stopped and the test marked as failed if this is exceeded. Set this to
   something a bit longer than the typical runtime to detect when new changes have
   significantly altered the runtime. If not given, the max runtime is infinity.
 * ncpus: the number of processes with which to run the simulation.

.. code-block:: python

   from answer_testing import EnzoETest

   class TestHLLCCloud(EnzoETest):
       parameter_file = "vlct/dual_energy_cloud/hllc_cloud.in"
       max_runtime = 30
       ncpus = 1

Tests involving ``Grackle``
###########################

If the class is associated with a test simulation that invokes ``Grackle``, you need to annotate the class declaration with the ``uses_grackle`` decorator.

.. code-block:: python

   from answer_testing import EnzoETest, uses_grackle

   @uses_grackle
   class TestGrackleGeneral(EnzoETest):
       ...

For all classes annotated with this decorator:

 * the framework knows that it must make symbolic links to all files in the directory run by ``GRACKLE_INPUT_DATA_DIR`` before it runs the simulation associated with this class.
 * the testing framwork also knows to skip the associated test(s) if the ``GRACKLE_INPUT_DATA_DIR`` environment variable is unset.
 
If you forget to add this label, ``Enzo-E`` will not be able to locate the data file needed for Grackle (in a portable way).
Thus, the associated simulation and test will fail.

Creating the Test Function
^^^^^^^^^^^^^^^^^^^^^^^^^^

The code above configures the simulation associated with the test. The next step
is to write a function which will be run after the simulation completes
successfully. This is done by creating a class method within the test class. This
function should only take the argument ``self`` (because it's a class method) and
nothing else. The function will be run from within the directory where the
simulation was run, so it will be able to load any files that were output.

.. code-block:: python

   def test_hllc_cloud(self):
       fn = "hllc_cloud_0.0625/hllc_cloud_0.0625.block_list"
       assert os.path.exists(fn)

Tests are typically implemented with an ``assert`` or related statement. In the
above example, we check for the existence of a file that should have been created
by the simulation. This is not specifically an answer test as we are not comparing
with results from another version of the code. However, these sorts of assertion
checks can be included in your test function if they are useful for verifying
proper running of the code.

Creating an Answer Test Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To create an answer test that will automatically save data to files and compare
with other files, we make use of the ``ytdataset_test`` Python decorator, also
located in `test/answer_tests/answer_testing.py`.

.. code-block:: python

   from answer_testing import \
       EnzoETest, \
       ytdataset_test, \
       assert_array_rel_equal

We also import an assertion function that will check for relative closeness of
values in an array.

The ``ytdataset_test`` decorator can then be put immediately above the definition
of a test function. This wraps the test function in additional code that will save
test files and run comparisons. With the ``ytdataset_test``, one must also provide
a function that will perform the comparison of results.

.. code-block:: python

   @ytdataset_test(assert_array_rel_equal, decimals=8)
   def test_hllc_cloud(self):
       ds = yt.load("hllc_cloud_0.0625/hllc_cloud_0.0625.block_list")
       ad = ds.all_data()

       wfield = ("gas", "mass")
       data = {field[1]: ad.quantities.weighted_standard_deviation(field, wfield)
                for field in ds.field_list}

       return data

When using ``ytdataset_test`` decorator, **a test function must return a dictionary
of values.** The values in the dictionary can be anything, e.g., numbers, string,
arrays, etc. In the above example, we load a snapshot with ``yt`` and compute the
weighted average and standard deviation (the ``weighted_standard_deviation`` function
returns both) of all the fields on disk. We now only need to return that and the
``ytdataset_test`` wrapper will save a file named after the test function (in this
case, 'test_hllc_cloud.h5' and will use the ``assert_array_rel_equal`` function to
check that results agree to within 8 decimal places. Note, the NumPy
`testing <https://numpy.org/doc/stable/reference/routines.testing.html>`__ module
defines several other assertion functions which may be useful.

Including Additional Configuration Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to communicate additional configuration options is through
environment variables. Once an environment variable is set (i.e., with ``export``
in bash), it can be seen by your test using the ``os.environ`` dict. Below, we
use the USE_DOUBLE environment variable to determine whether ``enzo-e`` was compiled
in single or double precision, and adjust the tolerance on the tests accordingly.

.. code-block:: python

   import os

   use_double = os.environ.get("USE_DOUBLE", "false").lower() == "true"
   if use_double:
       decimals = 12
   else:
       decimals = 6

   # inside the TestHLLCCloud class
   @ytdataset_test(assert_array_rel_equal, decimals=decimals)
   def test_hllc_cloud(self):
       ...

.. note::

   The above code is primarily for the sake of example.
   In practice, we now automatically detect the code's precision from the enzo-e executable.

Alternatively, additional configuration options can be configured through new command-line flags, which are introduced and parsed by the `conftest.py` file in the `answer_test` directory.
This is generally more robust than adding environment variables (since the flags are more easily discovered and are more explicit).
But, in practice it's made slightly more complicated by the fact that flags are parsed with pytest hooks.
Flags added in this way work best with ``pytest`` fixtures, while our tests mostly leverage features from `Python's unittest module <https://docs.python.org/3/library/unittest.html>`_.


Caveats
^^^^^^^

Below are a few things to keep in mind when designing new tests.

Defining Multiple Test Functions within a Class
###############################################

Multiple test functions can be implemented within the same answer test class.
However, the test simulation will be run **for each test**. If you want to
perform multiple checks on a long running simulation, it is a better idea to
implement them all with separate asserts inside a single function.

Answer Test Functions Must Have Unique Names
############################################

Answer test functions that use the ``ytdataset_test`` wrapper must all have unique
names. This is because each results file will be named with the name of the function
itself.
