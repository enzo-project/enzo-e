-------
BB Test
-------

Runs the "BB Test" problem as described in Federrath et al 2010, ApJ, 713, 269, and tests
for mass conservation.

The initial conditions include isothermal gas, with the gas having a constant small
"external density" outside of the truncation radius. Within the truncation radius, the gas
density has the following form:

:math:`\rho(\phi) = \rho_0 (1 + A \cos(2 \phi)),`

where :math:`\rho` is the gas density, :math:`\phi` is the azimuthal angle in the spherical polar
coordinate system, :math:`\rho_0` is the mean density and :math:`A` is the (small) fluctuation amplitude. In
addition, the gas rotates around the z-axis as a solid-body.

When the problem is left to run, the gas collapses into a disk, and sink particles form near the
center, which accrete gas and merge together.

This test can be executed by running "ctest -R bb_test" in the build directory.

All the files required for running the test can be found in input/bb_test. See
input/bb_test/README for further information.

