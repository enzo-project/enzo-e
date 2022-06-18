----------------
Grackle Tests
----------------

Tests for the method that invokes Grackle. These tests set up a
cooling test without hydrodynamics to run many one-zone models in
Grackle, fully sampling the density, temperature, and metallicity
parameter space over which the chemistry and cooling/heating tables
are valid.

method_grackle_general
======================

This test compares the summary statistics computed for several grackle
fields after a certain period of time to previously archived values.

The simulation timesteps are much larger that the cooling/heating.
This makes it more likely that separate processing elements will
execute grackle routines at the same time (thus increasing the chances
of exposing hypothetical problems related to Grackle & SMP mode).

This test is somewhat fragile given that upgrading Grackle versions could
conceivably alter the field values. In the future it would be better to
replace this with a test that:

1. checks out and builds a previous commit of Enzo-E
2. runs the simulation and saves the exact field values after running the simulations
3. checks out and builds a newer commit of Enzo-E (while leaving the build of Grackle unchanged)
4. runs the simulation and confirms that the Grackle related field values are identical to the field values from the earlier simulation.

method_grackle_cooling_dt
=========================

This test runs Grackle for a fixed number of cycles, and compares the
final simulation time to a reference value. Each simulation timestep
is set fraction of the minimum magnitude of the cooling/heating
timestep.
