-------
Testing
-------

Enzo-E currently uses two different testing frameworks for its test suite,
:ref:`ctest <ctest>` and :ref:`pytest <pytest>`. Both frameworks are run on
all pull requests and each has its advantages. You are welcome to use either
when creating new tests. The ``ctest`` framework is used primarily for unit
tests (i.e., tests that run from self-contained binaries), but can also be
used for tests that run the ``enzo-e`` binary. The ``pytest`` framework is
mainly useful for answer tests, where results from two versions of an
``enzo-e`` simulation are compared.

.. toctree::
   ctest_framework
   ctest_new
   pytest_framework
   existing_tests
