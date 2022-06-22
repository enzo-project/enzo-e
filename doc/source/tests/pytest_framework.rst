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
