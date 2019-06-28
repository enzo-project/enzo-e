.. _using-methods:

**************
Enzo-E Methods
**************

*[ This page is under development ]*
  
This section decribes all Methods available in Enzo-E: what they do,
what method-specific parameters they have, what fields / particles
they access or update, what solvers if any used, and how different
Methods are intended to be used with each other.

Current Enzo-E methods available include those listed below.

``"ppm"``: hydrodynamics
========================

This implements the modified piecewise parabolic method (PPM) in Enzo.

parameters
----------

.. list-table:: Method ``ppm`` parameters
   :widths: 10 5 1 30
   :header-rows: 1
   
   * - Parameter
     - Type
     - Default
     - Description
   * - ``"density_floor"``
     - `float`
     - `1.0e-6`
     - `Lower limit on density`
   * - ``"diffusion"``
     - `logical`
     - `false`
     - `PPM diffusion parameter`
   * - ``"dual_energy"``
     - `logical`
     - `false`
     - `Whether to use dual-energy formalism`
   * - ``"dual_energy_eta_1"``
     - `float`
     - `0.001`
     - `Dual energy parameter eta 1`
   * - ``"dual_energy_eta_2"``
     - `float`
     - `0.1`
     - `Dual energy parameter eta 2`
   * - ``"flattening"``
     - `integer`
     - `3`
     - `PPM flattening parameter`
   * - ``"minimum_pressure_support_parameter"``
     - `integer`
     - `100`
     - `Enzo's MinimumPressureSupportParameter`
   * - ``"number_density_floor"``
     - `float`
     - `1.0e-6`
     - `Lower limit on number density`
   * - ``"pressure_floor"``
     - `float`
     - `1.0e-6`
     - `Lower limit on pressure`
   * - ``"pressure_free"``
     - `logical`
     - `false`
     - `Pressure-free flag`
   * - ``"steepening"``
     - `logical`
     - `false`
     - `PPM steepening parameter`
   * - ``"temperature_floor"``
     - `float`
     - `1.0e-6`
     - `Lower limit on temperature`
   * - ``"use_minimum_pressure_support"``
     - `logical`
     - `false`
     - `Minimum pressure support`
   * - ``"mol_weight"``
     - `float`
     - `0.6`
     - `Mean molecular mass`

fields
------

.. list-table:: Method ``ppm`` fields
   :widths: 5 5 1 30
   :header-rows: 1

   * - Field
     - Type
     - Read/Write
     - Description   
   * - ``"density"``
     - ``enzo_float``
     - [rw]
     -    
   * - ``"velocity_x"``
     - ``enzo_float``
     - [rw]
     -
   * - ``"velocity_y"``
     - ``enzo_float``
     - [rw]
     - if rank ≥ 2
   * - ``"velocity_z"``
     - ``enzo_float``
     -  [rw]
     - if rank ≥ 3
   * - ``"total_energy"``
     - ``enzo_float``
     - [rw]
     -
   * - ``"acceleration_x"``
     - ``enzo_float``
     - [r]
     - if gravity
   * - ``"acceleration_y"``
     - ``enzo_float``
     - [r]
     - if gravity and rank ≥ 2
   * - ``"acceleration_z"``
     - ``enzo_float``
     - [r]
     - if gravity and rank ≥ 3
   * - ``"pressure"``
     - ``enzo_float``
     - [w]
     - computed from ``total_energy``

``"ppml"``: MHD
===============

PPML ideal MHD solver

``"mhd_vlct"``: MHD
===================

VL + CT (van Leer + Constrained Transport) MHD solver based on the
unsplit Godunov method described by
`Stone & Gardiner (2009) <http://adsabs.harvard.edu/abs/2009NewA...14..139S>`_
.

The algorithm combines attributes similar to the MUSCL-Hancock with
the constrained transport method.  In short, the method includes two
main steps. First the values are integrated to the half
time-step. Then, the values at the half time-step are used to
caluclate the change in the quantities over the full time-step. The
primary representation of the magnetic fields are 3 face-centered
fields (one for each component), which lie on the cell face matching
the dimension of the component (e.g. the x component lies on the
x-face). The method also tracks cell-centered components of the bfield
which are just the averages of the face-centered components. The
implementation has been modified to perform reconstruction using
internal energy inplace of pressure (in typical Enzo fashion).

This method can effectively function as a hydrodynamics method if
the magnetic fields are all initialized to zero (if they start at
zero, they will never be modified to non-zero values).

Note that the courant factor should be less than 0.5.

This method only support 3 dimensions.


parameters
----------

.. list-table:: Method ``mhd_vlct`` parameters
   :widths: 10 5 1 30
   :header-rows: 1
   
   * - Parameter
     - Type
     - Default
     - Description
   * - ``"density_floor"``
     - `float`
     - `none`
     - `Lower limit on density (must exceed 0)`
   * - ``"pressure_floor"``
     - `logical`
     - `none`
     - `Lower limit on thermal pressure (must exceed 0)`
   * - ``"riemann_solver"``
     - `string`
     - `hlld`
     - `name of the Riemann Solver to use`
   * - ``"half_dt_reconstruct_method"``
     - `string`
     - `nn`
     - `Reconstruction method for half timestep`
   * - ``"full_dt_reconstruct_method"``
     - `string`
     - `plm`
     - `Reconstruction method for full timestep`


fields
------

.. list-table:: Method ``mhd_vlct`` fields
   :widths: 5 5 1 30
   :header-rows: 1

   * - Field
     - Type
     - Read/Write
     - Description   
   * - ``"density"``
     - ``enzo_float``
     - [rw]
     -    
   * - ``"velocity_x"``
     - ``enzo_float``
     - [rw]
     -
   * - ``"velocity_y"``
     - ``enzo_float``
     - [rw]
     -
   * - ``"velocity_z"``
     - ``enzo_float``
     - [rw]
     -
   * - ``"total_energy"``
     - ``enzo_float``
     - [rw]
     - specific total energy
   * - ``"bfield_x"``
     - ``enzo_float``
     - [rw]
     - Cell-centered bfield (average of the corresponding ``bfieldi_x``)
   * - ``"bfield_y"``
     - ``enzo_float``
     - [rw]
     - Cell-centered bfield (average of the corresponding ``bfieldi_y``)
   * - ``"bfield_z"``
     - ``enzo_float``
     - [rw]
     - Cell-centered bfield (average of the corresponding ``bfieldi_z``)
   * - ``"bfieldi_x"``
     - ``enzo_float``
     - [rw]
     - Primary representation of x-component of bfield (lies on x-faces).
   * - ``"bfieldi_y"``
     - ``enzo_float``
     - [rw]
     - Primary representation of y-component of bfield (lies on y-faces).
   * - ``"bfieldi_z"``
     - ``enzo_float``
     - [rw]
     - Primary representation of z-component of bfield (lies on z-faces).
   * - ``"pressure"``
     - ``enzo_float``
     - [w]
     - computed from ``total_energy``
       
At initialization the face-centered magnetic fields should be divergence
free. For non-trivial configurations, the values should be initialized
from the magnetic vector potential and then each component of the
cell-centered bfields should then be initialized by averaging the
face-centered values of the corresponding component. To help facillitate
the latter step we provide the ``vlct_bfield`` initializer.


.. _using-vlct-reconstruction:

reconstruction
--------------

This subsection details the available interpolation methods for
reconstructing the left and right states of the cell-centered
interfaces. Presently, all available methods perform reconstruction
on cell-centered primitive quantites,
:math:`{\bf w} = (\rho, {\bf v}, p, {\bf B})`

To simplify the determination of the necessary number of ghost
zones for a given combination of reconstruction algorithms on
a unigrid mesh, we define the concepts of "stale depth" and
"staling rate". We define a "stale" value as a value that needs
to be refreshed. "Stale depth" indicates the number of field
entries, starting from the outermost field values on a block,
that the region including "stale" values extends over. Every
time quantities are updated over a (partial/full) timestep,
the stale depth increases. We define the amount by which it
increases as the "staling rate" (which depends on the choice
of interpolation method).

For a unigrid simulation, the number of required ghost zones
is given by the sum of the staling rates for each selected
reconstruction method.

We provide the names used to specify each available method in
the input file, the associated staling depth, and a brief
description.

.. list-table:: Available ``mhd_vlct`` reconstructors
   :widths: 3 1 30
   :header-rows: 1
   
   * - Name
     - Staling Depth
     - Description
   * - ``"nn"``
     - `1`
     - `Nearest Neighbor - (1st order) reconstruction of primitives`
   * - ``"plm"``
     - `2`
     - `Piecwise Linear Method - (2nd order) reconstruction of primitives (usinga minmod limiter)`

We provide a few notes about the choice of interpolator for this algorithm:

   * The recommended choices of reconstruction algorithms are ``"nn"`` for the
     half-timestep and then ``"plm"`` for the full-timestep. Using ``"nn"``
     both times also works, however errors tests show that errors arise when
     ``"plm"`` is used both times.
   * It is supposed to be possible to reconstruct the characteristic quantities
     for this method or to use higher order reconstruction in place of ``"plm"``
   * Reconstruction is always performed on the cell-centered magnetic fields.
     After reconstructing values along a given axis, the values of the
     reconstructed magnetic field component for that axis are replaced by the
     face-centered magnetic field values.

.. _using-vlct-riemann-solver:

riemann solvers
---------------

This subsection details the available Riemann Solvers. Currently all
available Riemann Solvers are defined to use magnetic fields, however,
they all appropriately handle the cases where the magnetic fields are
unformly 0. We provide a list of the names used to specify each
Riemann Solver in the input file, and a brief description for each of
them:

  * ``"hlle"`` (or equivalently ``"enzo_hlle"``). The HLLE
    approximate Riemann solver with wavespeeds in the same way as
    Enzo. Currently this raises an error as it is not tested.
  * ``"athena_hlle"`` The HLLE approximate Riemann solver with
    wavespeeds computed using the procedure from 
    `Stone et al. (2008)
    <http://http://adsabs.harvard.edu/abs/2008ApJS..178..137S>`_ 
    (the minimum and maximum eigenvalues of Roe's matrix are allowed
    to be wavespeeds).
  * ``"hlld"`` The HLLD approximate Riemann solver.



``"pm_deposit"``: particle-mesh
===============================

Particle-mesh ("PM") method component to deposit of field and particle
mass into a "total density" field

parameters
----------

.. list-table:: Method ``ppm`` parameters
   :widths: 10 5 1 30
   :header-rows: 1
   
   * - Parameter
     - Type
     - Default
     - Description
   * - ``"alpha"``
     - `float`
     - `0.5`
     - `Deposit mass at time t + alpha * dt`

fields
------

particles
---------
   
``"pm_update"``: particle-mesh
==============================

Particle-mesh ("PM") method component to update particle positions
given acceleration fields
   
``"heat"``: heat equation
=========================

A sample Method for implementing forward-euler to solve the heat equation.   
   
``"grackle"``: chemistry/cooling
================================

Calls methods provided by the external Grackle 3.0 chemistry and
cooling library.
   
``"comoving_expansion"``: comoving expansion
============================================

Adds the comoving expansion terms to the physical variables.
   
``"turbulence"``: driving
=========================

Turbulence driving.

``"gravity"``: particle-mesh
============================

Particle-mesh ("PM") method component to compute gravitational
potential given a total density field, and calculate associated
acceleration fields.
   
``"trace"``: tracer particles
=============================

Moves tracer particles given the velocity field.    
   

