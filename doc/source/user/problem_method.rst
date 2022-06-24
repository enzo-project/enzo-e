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
   * - ``"diffusion"``
     - `logical`
     - `false`
     - `PPM diffusion parameter`
   * - ``"flattening"``
     - `integer`
     - `3`
     - `PPM flattening parameter`
   * - ``"minimum_pressure_support_parameter"``
     - `integer`
     - `100`
     - `Enzo's MinimumPressureSupportParameter`
   * - ``"pressure_free"``
     - `logical`
     - `false`
     - `Pressure-free flag`
   * - ``"steepening"``
     - `logical`
     - `false`
     - `PPM steepening parameter`
   * - ``"use_minimum_pressure_support"``
     - `logical`
     - `false`
     - `Minimum pressure support`

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

``fluid_props`` compatability
-----------------------------

This method is also compatible with the ``"bryan95"`` dual-energy formalism.
See :ref:`using-fluid_props-de` for additional details.

This method currently ignores all of the floor parameters that are set in the ``physics:fluid_props:floors`` section of the parameter file.


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

``fluid_props`` compatability
-----------------------------

This method makes use of the ``density`` and ``pressure`` floor parameters that are set in the ``physics:fluid_props:floors`` section of the parameter file.
See :ref:`using-fluid_props-floors` for more details about specifying these parameters.
This method requires that both parameters are specified and that they have positive values.

This method is also compatible with the ``"modern"`` dual-energy formalism.
See :ref:`using-fluid_props-de` for additional details.

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

For a given particle type to be deposited to the total density field,
it must be part of the ``"is_gravitating"`` group, and must have either
an attribute called ``"mass"``, or a constant called ``"mass"``, but
not both.

If ``"mass"`` is an attribute, we loop through the mass attribute array
to get the mass of each particle; and if ``"mass"`` is a constant with a
value specified in the input parameter file, the mass of each particle is
equal to this value. In either case, the value of the divided by the cell
volume to get a density quantity, which is deposited on to the grid via
a CIC interpolation scheme.
   
``"pm_update"``: particle-mesh
==============================

Particle-mesh ("PM") method component to update particle positions
given acceleration fields. Only particle types in the ``"is_gravitating"``
group are updated.
   
``"heat"``: heat equation
=========================

A sample Method for implementing forward-euler to solve the heat equation.   
   
``"grackle"``: chemistry/cooling
================================

Calls methods provided by the external Grackle 3.0 chemistry and
cooling library.

Compatability with hydro/mhd solvers
------------------------------------

The ``"grackle"`` method is compatible with both the ``"ppm"`` and the
``"mhd_vlct"`` methods. The convention is to list the hydro method
before ``"grackle"`` in the ``Field:list`` parameter.  This
configuration performs advection and radiative cooling in an
operator-split manner (*Note: there isn't currently support for
performing radiative cooling during the predictor step of the
VL+CT solver*).

Integration with hydro-solvers is self-consistent when
``Field:Grackle:primordial_chemistry`` has values of ``0`` or ``1``.
However, the integration is somewhat inconsistent when the parameter
exceeds ``1``. While users shouldn't be too concerned about this
latter scenario unless they are simulating conditions where
:math:`{\rm H}_2` makes up a significant fraction of the gas density,
we describe the inconsistencies in greater detail below.

When ``Field:Grackle:primordial_chemistry > 1``, the Grackle library
explicitly models chemistry involving :math:`{\rm H}_2` and how it
modifies the adiabtic index. Grackle's routines treat
:math:`\gamma_0`, the "nominal adiabatic index" specified by
``Field:gamma``, as the adiabatic index for all monatomic species
(this should be ``5.0/3.0``). To that end, Grackle supplies functions
that can effectively be represented as :math:`\gamma(e, n_{{\rm H}_2},
n_{\rm other})` and :math:`p(\rho, e, n_{{\rm H}_2}, n_{\rm
other})`. In these formulas:

- :math:`p`, :math:`\rho` and :math:`e` correspond to the quantities
  held by the ``pressure``, ``density`` and ``internal_energy``
  fields.  *(Note: the* :math:`\gamma` *function's dependence on*
  :math:`e` *accounts for the dependence of* :math:`\gamma_{{\rm
  H}_2}` *on temperature)*

- :math:`n_{{\rm H}_2}` specifies the number density of
  :math:`{\rm H}_2`. :math:`n_{\rm other}` specifies a selection of
  the other primordial species (that roughly approximate the total
  number density). In practice, these are computed from passively
  advected species fields.

There are a handful of locations within the ``"ppm"`` and
``"mhd_vlct"`` methods where this treatment is relevant:

1. **Computing the timestep:** each hydro/mhd
   method uses the :math:`p(\rho, e, n_{{\rm H}_2}, n_{\rm other})`
   function for the pressure values. However, they both use
   :math:`\gamma_0` in other places.

2. **Pre-reconstruction pressure calculation:** each hydro/mhd
   solver internally computes the pressure that is to be reconstructed
   with :math:`p=(\gamma_0 - 1)e\rho`.

3. **Riemann Solver:** in each hydro/mhd solver, the Riemann Solver
   completely ignore the grackle supplied functions.

4. **VL+CT Energy floor and DE synchronization:** the internal energy
   floor is computed from the pressure floor using: :math:`e_{\rm
   floor} = \frac{p_{\rm floor}}{(\gamma_0 - 1)\rho}` (thus,
   :math:`p_{\rm floor}` may exceed :math:`p(\rho, e_{\rm floor},
   \ldots)`). Additionally, synchronizing the internal energy with
   total energy relies on :math:`\gamma_0`.

5. **PPM reconstruction:** uses :math:`\gamma_0`.

   
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


``"merge_sinks"``: merge sinks
==============================

Merges together sink particles which are separated by less than a given
"merging radius". This is done by copying all sink particles to / from
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
all the remaining sink particles which are 'out-of-bounds' of the block.

Star particles must have an attribute called ``"mass"`` if this method
is used.

This procedure cannot handle the case where particles originally
from non-neighbouring blocks are put into the same FoF group. If this is
found to occur, the program stops and prints an error message. This situation
is unlikely to happen, unless the merging radius is too large relative
to the block size.

Currently this will only run in unigrid mode. This is because this method
will only work correctly if all blocks containing sink particles are of the
same size, or equivalently, on the same refinement level.
For this reason, there is a check in the constructor of EnzoMethodMergeSinks
for whether ``"Adapt: max_level"`` is equal to zero. In future, we plan to
implement an accretion method, which will require a refinement condition that
any block containing a sink particle needs to be on the highest level of
refinement. In this case, the assumption that
blocks containing sink particles are all on the same level of refinement
would be valid.

WARNING: there is currently a memory leak issue when running with this method
which can cause Enzo-E to crash in mysterious ways. If this problem is
encountered, it is advised to increase the batch size parameter
(``"Particle:batch_size"``) by a factor of a few
before attempting to run again. To be completely safe, the user can set a
batch size larger than the total number of sink particles in the whole
simulation, which should be feasible for small test problems.

parameters
----------

.. list-table:: Method ``merge_sinks`` parameters
   :widths: 10 5 1 30
   :header-rows: 1
   
   * - Parameter
     - Type
     - Default
     - Description
   * - ``"merging_radius_cells"``
     - `float`
     - `8.0`
     - `The merging radius in units of the minimum cell width (i.e.,
       the minimum across all 3 dimensions), at the highest refinement
       level.`

   

