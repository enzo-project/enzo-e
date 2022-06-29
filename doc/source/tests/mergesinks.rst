------------------
Merge Sinks Tests
------------------

Four tests of the enzo-e "merge sinks" method.

- the stationary serial test

- the stationary parallel test
    
- the drifting serial test

- the stationary parallel test
    
All the files necessary for running the tests can be found in input/merge_sinks.

Each test involves setting up initial conditions with 1000 sink particles randomly
distributed in a sphere, with velocities directed towards the centre of the sphere.
In the "drifting" tests, an additional uniform velocity is added to all the
particles, and in the "stationary" tests, there is no additional velocity.

For the serial tests, Enzo-E runs in serial mode.

For the parallel tests, Enzo-E is run in parallel on four PEs.

The methods used are "pm update" and "merge sinks", so that the particles all
move with constant velocity until they are merged together.

The python script input/merge_sinks/mass_momentum_conservation.py then tests if
mass and momentum are conserved across all the output snapshots, to within some
tolerance, which depends on the floating-point precision with which Enzo-E was compiled
(1.0e-4 for single precision, 1.0e-6 for double precision).
If all quantities are indeed conserved, the test passes, if not, the test fails.

These tests can be executed by running "ctest -R merge_sinks" in the build directory.
See input/merge_sinks/README for further information on these tests.
