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

.. _using-grackle-gamma-with-HD:

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
``Method:Grackle:primordial_chemistry`` has values of ``0`` or ``1``.
However, the integration is somewhat inconsistent when the parameter
exceeds ``1``. While users shouldn't be too concerned about this
latter scenario unless they are simulating conditions where
:math:`{\rm H}_2` makes up a significant fraction of the gas density,
we describe the inconsistencies in greater detail below.

When ``Method:Grackle:primordial_chemistry > 1``, the Grackle library
explicitly models chemistry involving :math:`{\rm H}_2` and how it
modifies the adiabtic index. Grackle's routines treat
:math:`\gamma_0`, the "nominal adiabatic index" specified by
``Physics:fluid_props:eos:gamma``, as the adiabatic index for all
monatomic species (this should be ``5.0/3.0``). To that end, Grackle
supplies functions that can effectively be represented as
:math:`\gamma(e, n_{{\rm H}_2}, n_{\rm other})` and :math:`p(\rho, e,
n_{{\rm H}_2}, n_{\rm other})`. In these formulas:

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
   :math:`\gamma_0` in other places (such as the occurence of
   adiabatic index in the sound speed formula).

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

This method requires sink particles to have the following attributes: ``"mass"``, ``"x"``,
``"y"``, ``"z"``, ``"vx"``, ``"vy"``, ``"vz"``, ``"is_copy"``, ``"id"``, ``"lifetime"``,
and ``"creation_time"``. All these attributes must be of type ``"default"``, except for
``"is_copy"`` and ``"id"`` which must be of type ``"int64"``. Furthermore, ``"is_copy"``
must be initialized to 0 for all particles.

This method also requires that the number of root blocks across all axes is greater than
2, i.e., that ``"Mesh:root_blocks" = [a,b,c]``, where ``a``, ``b``, and ``c`` are all
greater than 2.

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
implement a refinement condition that
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


``"accretion"``: accretion
==============================

For cells within a spherical accretion zone around a sink particle, mass is removed
(i.e., the values of the density field are reduced) and added to the sink particle.
The momentum change of gas is in the accretion zone due to the mass loss is accounted
for by changing the momentum of the sink particle, so that total momentum is
conserved. The amount of mass removed is determined by which "flavor" of accretion is
chosen (specified by the ``"accretion:flavor"`` parameter), as well as the values
of the "density threshold" (specified by ``"accretion:physical_density_threshold_cgs"``) and the
"maximum mass fraction" (specified by ``"accretion:max_mass_fraction"``).

In ``"threshold"`` flavor accretion, the change in density of each cell is zero if the current
density is below the density threshold. If the current density is above the density threshold,
the change in density is the current density minus the density threshold, or the maximum mass
fraction times the current density, whichever is smaller.

In ``"bondi_hoyle"`` flavor accretion, the density change in each cell is calculated according
to the method described in Mark R. Krumholz et al 2004, ApJ, 611, 399. Furthermore, the
density change is limited in the same way as in ``"threshold"`` accretion.

In ``"flux"`` flavor accretion, the density change in each cell is calculated according to the
method described in Andreas Bleuler & Romain Teyssier 2004, MNRAS, 445, 4015-4036.
Furthermore, the density change is limited in the same way as in ``"threshold"`` accretion.

In ``"dummy"`` flavor accretion, no accretion is done (essentially, the accretion rate is zero).
This can be useful for testing purposes.

This method can only be used if ``"merge_sinks"`` is also used, with ``"merge_sinks"`` preceding
``"accretion"``. In addition, this method requires the use of three spatial dimensions.

This method requires the following fields (in addition to the fields required by the hydro
method): ``"density_source"``, ``"density_source_accumulate"``, ``"mom_dens_x_source"``,
``"mom_dens_x_source_accumulate"``, ``"mom_dens_y_source"``, ``"mom_dens_y_source_accumulate"``,
``"mom_dens_z_source"``, and ``mom_dens_z_source_accumulate"``. In addition, if sink particles
have a ``"metal_fraction"`` attribute, there must be a ``"metal_density"`` field.

This method also requires sink particles to have the following attributes: ``"mass"``, ``"x"``,
``"y"``, ``"z"``, ``"vx"``, ``"vy"``, ``"vz"``, and ``"accretion_rate"``, which must all be
of type ``"default"``.

parameters
----------

.. list-table:: Method ``accretion`` parameters
   :widths: 10 5 1 30
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``"accretion_radius_cells"``
     - `float`
     - `4.0`
     - `The accretion radius (i.e., the radius of the spherical accretion zone)
       in units of the minimum cell width (i.e., if the cell width along all the x, y, and
       z-axes are hx, hy, and hz, then the minimum cell width is the minimum of hx, hy, and hz),
       at the highest refinement level. Its value must be less than one fewer than the minimum
       ghost depth  for "flux" accretion, and less than the minimum ghost depth
       for other flavors of accretion. The ghost depth is 4 (along all axes) by default.`
   * - ``"flavor"``
     - `string`
     - ``""``
     - `The flavor of accretion used, which can be either "threshold", "bondi_hoyle", or "flux".
       If this parameter is not set in the parameter file, or if some other string is
       provided, then the accretion method will be called but will do nothing.`
   * - ``"physical_density_threshold_cgs"``
     - `float`
     - `1.0e-24`
     - `The value of the physical density threshold in cgs units. The density in each cell in
       the accretion zone cannot go below this value during the accretion process. The value of
       the density threshold in code units must be greater than or equal to the value of the density
       floor imposed by the hydro method.`
   * - ``"max_mass_fraction"``
     - `float`
     - `0.25`
     - `This parameter specifies the maximum fraction of mass which can be accreted from a cell
       in one timestep. This value of this parameter must be between 0 and 1.`


``"sink_maker"``: sink maker
==============================

This method runs on blocks at the highest level of refinement, and forms sink particles in cells
which satisy certain criteria.

First, the gas density in the cell must be larger than the
density threshold, which is specified by the ``"sink_maker:physical_density_threshold_cgs"``
parameter. If so, the mass of the potential sink particle is
:math:`V_{cell} \times \max(\rho - \rho_{thresh}, f_{max} \, \rho)`, where :math:`V_{cell}` is the
cell volume, :math:`\rho` is the cell gas density, :math:`\rho_{thresh}` is the density
threshold, and :math:`f_{max}` is the maximum fraction of the cell mass which can be turned
into a sink particle in one timestep, which is specified by the
``"sink_maker:maximum_mass_fraction"`` parameter. This mass must be greater
than the minimum sink mass, which is specified by ``"sink_maker:min_sink_mass_solar"`` (in solar
mass units).

Next, the local Jeans length :math:`\lambda_J` is calculated, where
:math:`\lambda_J = \frac{\pi c_s^2}{G \rho}`, where :math:`c_s` is the sound speed of the gas
in the cell, and :math:`G` is the gravitational constant. It is then checked whether
:math:`\lambda_J < N_J \times h_{max}`, where :math:`N_J` is specified by
``"sink_maker:jeans_length_resolution_cells"``, and :math:`h_{max} = max(h_x, h_y, h_z)`, where
:math:`h_x`, :math:`h_y`, and :math:`h_z` are the cell widths along the x-, y- and z-axes,
respectively.

The next check is that the flow is converging. This is done by computing the strain tensor,
given by :math:`A_{ij} = \frac{1}{2} \, \left( \frac{dv_i}{dx_j} + \frac{dv_j}{dx_i} \right)`.
Since this tensor / matrix is real and symmetric, it has three real eigenvalues, and the check
is equivalent to checking that all three eigenvalues are negative.

The final check is optional, i.e., it is only done if ``"sink_maker:check_density_maximum"``
is "true", and a cell will pass this check if it is a local density maximum, that is, its
density is larger than the density in all 26 neighboring cells.

If a cell passes all the checks that are performed, a sink particle is created. Its position
is the coordinates of the center of the cell, plus a small random offset. The maximum size
of the random offset is controlled by ``"sink_maker:max_offset_cell_fraction"``.

This method requires sink particles to have the following attributes: ``"mass"``, ``"x"``,
``"y"``, ``"z"``, ``"vx"``, ``"vy"``, ``"vz"``, and ``"creation_time"``, which must all be
of type ``"default"``; and ``"id"`` and ``"is_copy"``, which must be of type ``"int64"``.
If sink particles have a ``"metal_fraction"`` attribute, there must be a
``"metal_density"`` field.


parameters
----------

.. list-table:: Method ``sink_maker`` parameters
   :widths: 10 5 1 30
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``"jeans_length_resolution_cells"``
     - `float`
     - `4.0`
     - `If the local Jeans length in a cell is less than this quantity multiplied by the maximum
       cell width, then the cell is a candidate for forming a sink. The maximum cell width is
       maximum value out of hx, hy, and hz, where hx, hy, and hz are the cell widths across the
       x-, y- and z-axes, respectively.`
   * - ``"physical_density_threshold_cgs"``
     - `float`
     - `1.0e-24`
     - `The value of the physical density threshold in cgs units. The density in a cell must be
       greater than the density threshold to be able to form a sink. The density in a cell after
       sink formation will be no less than the density threshold. The value of
       the density threshold in code units must be greater than or equal to the value of the
       density floor imposed by the hydro method.`
   * - ``"max_mass_fraction"``
     - `float`
     - `0.25`
     - `The mass of a newly-formed sink is bounded above by this parameter multiplied by the cell
       density multiplied by the cell volume. The value of this parameter must be between
       0 and 1.`
   * - ``"min_sink_mass_solar"``
     - `float`
     - `0.0`
     - `The minimum mass of a newly-formed sink particle, in solar mass units.`
   * - ``"check_density_maximum"``
     - `logical`
     - `true`
     - `Determines whether a cell is required to be a local density maximum in order to form a
       sink particle.`
   * - ``"max_offset_cell_fraction"``
     - `float`
     - `0.0`
     - `When a cell creates a sink particle, the x/y/z coordinate of its initial position will be
       the x/y/z coordinate of the center of the cell, plus a random value generated from a
       uniform distribution on the interval [-A,A], where A is equal to
       this parameter multiplied by the cell width along the x/y/z axis.`
   * - ``"offset_seed_shift"``
     - `integer`
     - `0`
     - `When computing the random offset for the initial position of a sink particle, we compute
       an unsigned 64 bit integer value from the cycle number, the block index, and the cell
       index, and then add on this value to give the seed for the random number generator.`
