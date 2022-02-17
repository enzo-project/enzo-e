------------------
Merge Stars Tests
------------------

Two tests of the enzo-e "merge stars" method (the serial test and the parallel test). All the files necessary for running the tests can be found in input/MergeStars,
and the tests can be run by executing input/MergeStars/run.sh from the top-level directory of the repository.

Each test involves setting up initial conditions with 1000 star particles randomly distributed in a sphere, with velocities directed towards the centre of the sphere. The sphere as a whole also has an
additional constant drift velocity. For the serial test, Enzo-E runs in serial mode with the parameter file merge_stars_test_serial.in. For the parallel test, Enzo-E is run using charm with
4 processes with the parameter file merge_stars_test_parallel.in. The methods used are "pm update" and "merge stars", so that the particles all move with constant velocity until they are merged
together. The python script input/MergeStars/mass_momentum_conservation.py then tests if mass and momentum are conserved across all the output snapshots, to within some tolerance,
which depends on the precision specified by the $CELLO_PREC environment variable (1.0e-4 for single precision, 1.0e-6 for double precision).
If all quantities are indeed conserved, the test passes, if not, the test fails. See input/MergeStars/README for further details on running the tests.
