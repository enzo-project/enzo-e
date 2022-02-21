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

.. _vlct_overview:

``"mhd_vlct"``: hydrodynamics/MHD
=================================

This implements the VL + CT (van Leer + Constrained Transport) unsplit
Godunov method described by `Stone & Gardiner (2009)
<https://adsabs.harvard.edu/abs/2009NewA...14..139S>`_
. This solver operates in 2 modes (designated by the required
``Method:mhd_vlct:mhd_choice`` parameter):

  1. a pure hydrodynamical mode (it cannot handle magnetic fields)
  2. a MHD mode.

The algorithm is a predictor-corrector scheme, with attributes similar
to the MUSCL-Hancock method. For both modes, the method includes two
main steps. First the values are integrated to the half
time-step. Then, the values at the half time-step are used to
caluclate the change in the quantities over the full time-step.

In the MHD mode, the algorithm is combined with the constrained
transport method. The primary representation of the magnetic field,
:math:`\vec{B}`, components are stored in 3 face-centered Cello fields.
In more detail:

  - :math:`B_x` is stored on the x-faces
  - :math:`B_y` is stored on the y-faces
  - :math:`B_z` is stored on the z-faces

The method also tracks components of the magnetic fields at the cell-centers,
which just store the average of the values from the cell's faces.

Currently, this method only support 3-dimensional problems.  In the
future, alternative modes supporting MHD could easily be implemented
that use divergence-cleaning schemes instead of constrained transport.


parameters
----------

Note that the courant factor (specified by the ``"courant"``
parameter) should be less than 0.5.

.. list-table:: Method ``mhd_vlct`` parameters
   :widths: 10 5 1 30
   :header-rows: 1
   
   * - Parameter
     - Type
     - Default
     - Description
   * - ``"mhd_choice"``
     - `string`
     - `none`
     - `Specifies handling of magnetic fields (or lack thereof)`
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
   * - ``"theta_limiter"``
     - `float`
     - `1.5`
     - `controls dissipation of the "plm"/"plm_enzo" reconstruction
       method.`
   * - ``"dual_energy"``
     - `logical`
     - `false`
     - `Whether to use dual-energy formalism`
   * - ``"dual_energy_eta"``
     - `float`
     - `0.001`
     - `Dual energy parameter eta`


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
     - computed from ``total_energy`` (``internal_energy`` if dual-energy)
   * - ``internal_energy``
     - ``enzo_float``
     - [rw]
     - if dual-energy

In hydro-mode, none of the 6 fields used to store the magnetic field should
be defined.
       
At initialization the face-centered magnetic field should be
divergence free. Trivial configurations (e.g. a constant magnetic
field everywhere) can be provided with the ``"value"``
initializer. For non-trivial configurations, we have provide the
``"vlct_bfield"`` initializer which can initialize the magnetic fields
(face-centered and cell-centered) from expression(s) given in the
parameter file for component(s) of the magnetic vector potential.

.. _using-vlct-de:

dual-energy formalism
---------------------

The implementation of the dual-energy more closely resembles the
implementation employed in Enzo's Runge–Kutta integrator than the original
conception used by Enzo's ppm integrator, (for a description of that
implementation, see `Bryan et al (1995)
<https://ui.adsabs.harvard.edu/abs/1995CoPhC..89..149B>`_ ). There are 3
main differences from the original conception:

  1. internal energy is always used to compute pressure. In the original
     conception, pressure could be computed from ``total_energy`` or
     ``internal_energy`` (the decision was independent of synchronization).
  2. ``pressure`` and ``internal_energy`` are not separately reconstructed.
     Instead, just the pressure is reconstructed. The ``internal_energy``
     is computed at the left and right interfaces from the reconstructed
     quantities.
  3. Synchronization of the total and internal energies is a local
     operation that doesn't require knowledge of cell neighbors. In the
     original conception, knowledge of the immediate neighbors had been
     required (each synchronization incremented the stale depth - 3 extra
     ghost zones would have been required).

For clarity, the conditions for synchronization are provided below. The
specific ``internal_energy``, :math:`e`, is set to
:math:`e'= E - (v^2 + B^2/\rho)/2` (where :math:`E` is the specific
``total_energy``) when the following conditions are met:

  * :math:`c_s'^2 > \eta v^2`, where :math:`c_s'^2=\gamma(\gamma - 1) e'`.
  * :math:`c_s'^2 > \eta B^2/\rho` (this is always satisfied in hydro mode)
  * :math:`e' > e /2`

If the above condition is not met, then ``total_energy`` is set to
:math:`e + (v^2 + B^2/\rho)/2` in MHD mode (in hydro mode, it's set to
:math:`e + v^2/2`).
    
When ``"dual_energy_eta"``, is set to ``0``, :math:`e` is always set to
``e'``. This is done to provide support for Grackle (in the future)
without the dual-energy formalism.

*Note: in the future, the behavior described in difference 2, may change
to achieve better compatibility with Grackle.*

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
that the region encompassing "stale" values extends over. Every
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

.. list-table:: Available ``mhd_vlct`` reconstructors (and slope
		limiters)
   :widths: 3 1 30
   :header-rows: 1
   
   * - Name
     - Staling Depth
     - Description
   * - ``"nn"``
     - `1`
     - `Nearest Neighbor - (1st order) reconstruction of primitives`
   * - ``"plm"`` or ``"plm_enzo"``
     - `2`
     - `Piecwise Linear Method - (2nd order) reconstruction of
       primitives using the slope limiter from Enzo's Runge–Kutta
       integrator. This is tuned by the` ``"theta_limiter"``
       `parameter, which must satisfy` ``1 <= "theta_limiter" <=
       2``. `As in Enzo, the default value is 1.5. A value of 1 is the
       most dissipative and it is equivalent to the traditional minmod
       limiter. A value of 2 is the least dissipative and it
       corresponds to an MC limiter (monotized central-difference
       limiter).`
   * - ``"plm_athena"``
     - `2`
     - `Piecwise Linear Method - (2nd order) reconstruction of
       primitives using the slope limiter from Athena (& Athena++).
       For some primitive variable`, :math:`{\bf w}_{i}`, `the limited
       slope is defined in terms of the left- and right-differences:`
       :math:`\delta{\bf w}_{L,i}={\bf w}_{i}-{\bf w}_{i-1}` `and`
       :math:`\delta{\bf w}_{R,i}={\bf w}_{i-1}-{\bf w}_{i+1}`.  `If
       the signs of the differences don't match (or at least 1 is 0),
       then the limited slope is 0. Otherwise the limited slope is the
       harmonic mean of the differences.`

We provide a few notes about the choice of interpolator for this algorithm:

   * The recommended choices of reconstruction algorithms are ``"nn"`` for the
     half-timestep and then piecewise-linear reconstruction for the
     full-timestep (most test problems have been run using ``plm`` with
     ``theta_limiter=2``, matching the integrator description in
     `Stone & Gardiner 2009
     <https://adsabs.harvard.edu/abs/2009NewA...14..139S>`_ ). Using ``"nn"``
     both times also works, however tests show that errors arise when
     piecewise linear reconstruction is used both times.
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

  * ``"hll"`` The HLL approximate Riemann solver with wavespeeds
    bounds estimated as :math:`S_L = \min(u_L - a_L, u_R - a_R)` and
    :math:`S_R = \max(u_L + a_L, u_R + a_R)`. This is one of the
    proposed methods from Davis, 1988, SIAM J. Sci. and Stat. Comput.,
    9(3), 445–473. The same wavespeed estimator was used in MHD HLL
    solver implemented for Enzo's Runge Kutta solver. Currently, this
    has only been implemented for MHD mode and it will raise an error
    as it isn't tested.

  * ``"hlle"`` The HLLE approximate Riemann solver - the HLL solver
    with wavespeed bounds estimated according to
    Einfeldt, 1988, SJNA, 25(2), 294–318. This method allows the
    min/max eigenvalues of Roe's matrix to be wavespeed estimates. For a
    description of the procedure for MHD quantities, see
    `Stone et al. (2008)
    <https://adsabs.harvard.edu/abs/2008ApJS..178..137S>`_ .
    If using an HLL Riemann Solver, this is the recommended choice.
    Currently, this has only been implemented for MHD mode.

  * ``"hllc"`` The HLLC approximate Riemann solver.
    For an overview see Toro, 2009, *Riemann Solvers and Numerical
    Methods for Fluid Dynamics*, Springer-Verlag. This is a solver for
    hydrodynamical problems that models contact and shear waves. The
    wavespeed bounds are estimated according to the Einfeldt approach.
    This can only be used in hydro mode.
    
  * ``"hlld"`` The HLLD approximate Riemann solver described in
    Miyoshi & Kusano, 2005. JCP, 315, 344. The wavespeed bounds are
    estimated according to eqn 67 from the paper. This reduces to an
    HLLC Riemann Solver when magnetic fields are zero (the wavespeed
    bounds will differ from ``"hllc``). This can only be used in MHD
    mode.


.. note::

      When the dual-energy formalism is in use, all of the solvers treat
      the internal energy as a passively advected scalar.

      This is not completely self-consistent with the assumptions made by the
      HLLD solver. Unlike the other HLL-solvers which assume constant
      pressure in the intermediate regions of the Riemann Fan the HLLD solver
      assumes constant total pressure. It is unclear whether this causes any
      problems.

   


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


``"merge_stars"``: merge stars
==============================

Merges together star particles which are separated by less than a given
"merging radius". This is done by copying all star particles to / from
all neighbouring blocks. A Friend-of-Friends algorithm is used to
partition particles into groups, where all particles within a given group
are separated by less than a merging radius. If a group has more than one
particle, one of the particles has its properties changed: its position
becomes that of the centre-of-mass of the group, and it takes the total
mass, momentum and mass fraction of the whole group.
In addition, its 'lifetime' attribute is set to be the maximum of the group,
its 'creation_time' attribute is set to be the minimum of the group, and its
'id' attribute is set to the minimum of the group. Other particles in the
group are marked for deletion. The final step is for each block to delete
all the remaining star particles which are 'out-of-bounds' of the block.

This procedure cannot handle the case where particles originally
from non-neighbouring blocks are put into the same FoF group. If this is
found to occur, the program stops and prints an error message. This situation
is unlikely to happen, unless the merging radius is too large relative
to the block size.

Currently this will only run in unigrid mode. This is because when a block
calls this method, it assumes that all neighbouring blocks have the same size,
or equivalently, on the same refinement level.
For this reason, there is a check in the constructor of EnzoMethodMergeStars
for whether ``"Adapt: max_level"`` is equal to zero. In the future, we plan to
implement an accretion method, which will require a refinement condition that
any block containing an accreting particle, or neighbouring such a block, needs
to be on the highest level of refinement. In this case, the assumption that
blocks containing accreting particles (which we want to merge together), and
blocks neighbouring such blocks, are all on the same level of refinement
would be valid.

WARNING: there is currently a memory leak issue when running with this method
which can cause Enzo-E to crash in mysterious ways. If this problem is
encountered, it is advised to increase the batch size parameter
(``"Particle:batch_size"``) by a factor of a few
before attempting to run again. To be completely safe, the user can set a
batch size larger than the total number of star particles in the whole
simulation, which should be feasible for small test problems.

parameters
----------

.. list-table:: Method ``merge_stars`` parameters
   :widths: 10 5 1 30
   :header-rows: 1
   
   * - Parameter
     - Type
     - Default
     - Description
   * - ``"merging_radius_cells"``
     - `float`
     - `8.0`
     - `The merging radius relative to the cell-width`

   

