.. _pytest:

--------------------
The pytest Framework
--------------------

The ``pytest`` framework is primarily used for running answer tests, where a
simulation is run with two versions of ``enzo-e`` and their results are compared.
This is useful for testing problems with no analytical solution or generally
verifying that results from commonly run simulations don't drift.

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

Before the answer tests can be run, a few environment variables must be set to
configure behavior.

 * TEST_RESULTS_DIR: points to a directory in which answers will be stored
 * CHARM_PATH: points to the directory in which ``charmrun`` is located
 * GENERATE_TEST_RESULTS: "true" to generate test results, "false" to compare
   with existing results.

.. code-block:: bash

   $ export TEST_RESULTS_DIR=~/enzoe_tests
   $ export CHARM_PATH=~/local/charm-v7.0.0/bin

Generating Test Answers
^^^^^^^^^^^^^^^^^^^^^^^

First, check out the highest numbered gold standard tag and compile ``enzo-e``.

.. code-block:: bash

   $ git checkout gold-standard-1
   $ ...compile enzo-e

Then, configure the test suite to generate answers by setting 
GENERATE_TEST_RESULTS to true.

.. code-block:: bash

   $ export GENERATE_TEST_RESULTS=true

Finally, run the test suite by calling ``pytest`` with the answer test directory.

.. code-block:: bash

   $ pytest test/answer_tests
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
run with the test suite to check for problems. Set the GENERATE_TEST_RESULTS
environment variable to false to configure the test suite to compare with existing
answers.

.. code-block:: bash

   $ git checkout main
   $ ...compile enzo-e
   $ export GENERATE_TEST_RESULTS=false
   $ pytest test/answer_tests

Getting More Output from Pytest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, most output printed by ``enzo-e`` or the test scripts will be swallowed
by ``pytest``. When tests fail, the Python traceback may be shown, but not much
else. There are various flags to increase the verbosity of ``pytest``, but the
``-s`` flag will show all output, including from the simulation itself. The
``enzo-e`` answer test suite will also print out the values of all configuration
variables when this flag is given.

.. code-block:: bash

   $ pytest -s test/answer_tests

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
