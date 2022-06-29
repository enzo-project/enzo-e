------------------
Accretion Tests
------------------

Six tests of the enzo-e "accretion" method.

- the threshold accretion serial test

- the threshold accretion parallel test

- the Bondi-Hoyle accretion serial test

- the Bondi-Hoyle accretion parallel test

- the flux accretion serial test

- the flux accretion parallel test

  All the files necessary for running the tests can be found in input/accretion.

Each test involves setting up initial conditions with a sink particle with some initial mass,
position, and velocity, in a background medium of gas with a uniform density and pressure.
In the "threshold" and "Bondi-Hoyle" tests, the gas has an initial velocity of zero everywhere,
and inthe "flux" tests, the gas in every cell has a velocity with a uniform magnitude, directed
towards the initial position of the sink particle. In the "Bondi-Hoyle" tests, the parameters
are given values so that the initial Bondi-Hoyle radius is approximately equal to one cell width.

For the serial tests, Enzo-E runs in serial mode.

For the parallel tests, Enzo-E is run in parallel on four PEs.

The methods used are "mhd_vlct", "pm update", "merge sinks", and "accretion".

The python script input/accretion/mass_momentum_conservation.py tests if
mass and momentum are conserved across all the output snapshots, to within some
tolerance, which depends on the floating-point precision with which Enzo-E was compiled  (1.0e-6 for single precision, 1.0e-12 for double precision).
If all quantities are indeed conserved, the test passes, if not, the test fails.

These tests can be executed by running "ctest -R accretion" in the build directory.
See input/accretion/README for further information on these tests.
