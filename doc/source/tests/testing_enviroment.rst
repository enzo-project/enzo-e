----------------------
The Testing Enviroment
----------------------

The testing enviroment is does both unit and integration testing. In order to build a new test follow how to build a new test instructions and then call make test. In order to run the new test separately write ``charmrun +p4 bin/enzo-p input/NewTestDirectory/new_test.in`` in the main enzo directory. If charmrun command is not in your path charmrun's path must also be included:  ``~/Charm/bin/charmrun +p4 bin/enzo-p input/NewTestDirectory/new_test.in``.

The testing infrastructure does both unit and integration testing. Unit testing tests the components of enzo-e and integration testing tests the functionality of units working together.

The make test outputs first the unit test results in a count of passes and fails, and then a count of the passes and fails of the integration tests.
