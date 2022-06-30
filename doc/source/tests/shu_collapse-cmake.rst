-----------------
Shu Collapse Test
-----------------

Runs the Shu Collapse problem as described in Federrath et al 2010, ApJ, 713, 269, and tests
for mass conservation.

The initial conditions include static isothermal gas with an inverse square density profile within
a truncation radius, with the gas having a constant small "external density" outside of the
truncation radius. When the problem is left to run, a sink particle forms at the
center, and the gas collapses towards the center to be accreted by the sink particle.

This test can be executed by running "ctest -R shu_collapse" in the build directory.

All the files required for running the test can be found in input/shu_collapse. See
input/shu_collapse/README for further information.

