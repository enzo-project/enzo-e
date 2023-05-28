---------------
Cosmology Tests
---------------

Tests of the enzo-e cosmology methods using simple cosmology problems.

.. _cosmo-dd-multispecies-short:

cosmo-dd-multispecies-short
===========================
A simple test problem that runs the PPM, Cosmology, Gravity (with DD-solver), and Grackle (with multiple species) for a little more than 100 cycles (so that the domain has the opportunity for refinement).


method_cosmology-1
==================

Runs PPM and Cosmology methods with cosmology variables set in the Physics parameter.
This outputs image files for: density field, x velocity, y velocity, z velocity, x acceleration, y acceleration, z acceleration, total density and potential image files.
It also outputs an hdf5 data file.

This is a deprecated problem that will be removed in the near future.

method_cosmology-8
==================

This is the same as the preeceding problem, but uses 8 root blocks instead of 1.

This is a deprecated problem that will be removed in the near future.
