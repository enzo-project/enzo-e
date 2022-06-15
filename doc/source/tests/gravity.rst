-------------
Gravity Tests
-------------

Tests for gravity methods.

Completion-Time Tests
=====================

Runs gravity method using cg solver. Outputs mesh, phi, rho, x
acceleration, y acceleration image files. Tests succeed or fail based
on the final completion time of the simulation.

method_gravity_cg-1
~~~~~~~~~~~~~~~~~~~

Tests ``EnzoMethodGravityCg`` at P=1 in 2D.


method_gravity_cg-8
~~~~~~~~~~~~~~~~~~~

Tests ``EnzoMethodGravityCg`` at P=8 in 2D.

Analytic-Tests
==============

These tests are associated with a dedicated python-script that runs a
test problem with a known analytic solution.

stable_jeans_wave
~~~~~~~~~~~~~~~~~

The basic setup for this test is to initialize a 1D stable Jeans wave
in a 3D box that is inclined with respect all 3 Cartesian axes. Enzo-E
writes the initial snapshot to disk, evolves the wave over a complete
period, and then writes the final snapshot result to disk. Then
the python script computes the :math:`L_1` error norm computed by
comparing values in the initial and final snapshots (the analytic
solution expects them to be identical) and compares it with a
precomputed value.

For this test to pass, the gravity solver and hydro solver's gravity
source terms must each be implemented correctly. At this time, this
test has only been defined to use the VL+CT solver (without magnetic
fields) and the cg-solver. The python script invokes the test at 2
resolutions (16 cells per wavelength and 32 cells per wavelength).

To invoke this test, invoke the following command from the root
directory of the repository:

.. code-block:: bash

   python input/Gravity/run_stable_jeans_wave_test.py \
     --launch_cmd <path-to-enzoe-binary>

where ``<path-to-enzoe-binary>`` is replaced with the path to the Enzo-E
binary that you want to test. The test prints a summary of the results to
stdout and an exit code of of 0 indicates that all subtests passed.
