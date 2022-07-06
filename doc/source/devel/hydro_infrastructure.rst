****************************
Hydro/MHD C++ Infrastructure
****************************

*[This page is under development]*

In this section we discuss some of the C++ infrastructure provided in
the Enzo Layer that can be optionally used to implement other
hydro/MHD methods. The infrastructure was used to implement other the
VL + CT MHD solver.

*Note: Currently, barotropic equations of state are not yet implemented
within the infrastucture. However they are mentioned throughout this
guide and slots have been explicitly left open for them to be
implemented within the framework.*

*Note: Currently brief summaries of the interfaces of each of the
objects are provided below. More detailed descriptions are provided
in the header files using doxygen documentation. If we end up
producing reference documentation from doxygen (via breathe), then
the interface summaries should be reduced or deleted.*

===============
Shorthand Terms
===============

We briefly define a few terms that are used throughout the
documentation and codebase

Integration/Primitive Quantities
-------------------------------------

Throughout this guide and the relevant sections of the codebase, we
categorize quantities as integration quantities and primitives.

The integration quantities are the cell-centered quantities that
Enzo-E, as a whole, uses to describe the state of the fluid. In other
words, these are the hydro/MHD quantities that Enzo-E expects to be
integrated from one time-step to the next. The integration quantities
include all passive scalars in their "conserved" form (i.e. as
densities).  All integration quantities can basically be
subcategorized as either "conserved" or "specific" (a specific
quantity like velocity that becomes conserved after multiplication by
density). In cases using the dual energy formalism, we treat the
specific internal energy as an integration quantity even though the
internal energy density is not technically conserved (it requires
source terms).

The primitive quantities follow the normal textbook definition. We use
the primitives internally within the hydro/MHD solver for
reconstruction. The primitives include all passive scalars in their
"specific" form (i.e. as mass fractions). Note that some quantities
(like density or velocity) are categorized as both an integration
quantity and a primitive.

To provide a more concrete example, we categorize the quantities related
to an ideal, adiabatic gas:

  * density - integration and primitive

  * velocity - integration and primitive

  * pressure - only primitive

  * (specific) total energy - only integration

  * magnetic field - integration and primitive

If using the dual energy formalism, we would categorize the (specific)
internal energy as an integration quantity. As of now, there wouldn't
be a primitive counterpart to the internal energy since it can be
computed from the reconstructed density and pressure.

stale depth
-----------

To help simplify the implementation of several operations, we
introduce the concept of "stale depth". At the start of a time-step,
there are never any "stale" values ("stale depth" is zero). However,
every time the flux divergence get's added to any quantities, the
outermost layer of up-to-date quantities (in the ghost zone) becomes
invalid; this happens because the fluxes on the exterior faces of
that layer are not accurately known. We refer to these invalid values
as "stale" and say that the stale depth increased by 1. The stale
depth can also be incremented by other operations (e.g. piecewise
linear reconstruction). At the end of a time-step, the stale depth
should be equal to or less than the ghost depth so that the all of the
"stale" values will be refreshed (resetting the "stale depth" to
zero).

We formally define "stale depth" as the number of the layers of the
outermost field entries that include "stale" values (for cell-centered
fields this is the number of cells from the edge that include stale
values). For a given stale depth:

  * A face-centered field that has values on the exterior of a mesh
    block will always have one more unstaled value along the axis of
    face-centering than a cell-centered field.

  * A face-centered field that doesn't have values on the exterior of
    a mesh block will always have one less unstaled value along the
    axis of face-centering than a cell-centered field.

The introduction of this formalism has 2 key benefits:

  1. Simplifies calculation of the required ghost depth.

  2. When used alongside ``CelloArray``, it drastically simplifies the
     determination of which indices to iterate over. The
     ``EnzoFieldArrayFactory`` can take the stale depth as an argument
     in its constructor and then all arrays that an instance builds
     will have the stale values clipped off.  This allows the bounds of 
     for-loops to be written as though the only reconstruction algorithm
     is nearest-neighbor interpolation and as though there never any
     preceeding partial timesteps.


Finally, for each reconstruction algorithm we define a staling
rate. This specifies the amount that the stale depth increases each
time values are updated over a (partial) timestep using the
algorithm. This is subdivided into

  1. "immediate staling rate," which species the amount by whch the
     stale-depth increases immediately after reconstruction (e.g. this is
     0 for nearest-neighbor interpolation and 1 for piecewise-linear
     interpolation).

  2. "delayed staling rate," which specifies the amount by which the
     stale depth increase after adding the flux divergence (e.g. 1 for
     both nearest-neighbor and piecewise linear interpolation).


.. _Centered-Field-Registry:

=======================
Centered Field Registry
=======================

The Hydro/MHD infrastructure helped motivate the creation of
``EnzoCenteredFieldRegistry``, to encapsulate a static (at runtime)
registry of all known fields used by the Enzo section of the
codebase and track some basic meta-data about the fields. Although it's
functionallity is presently limited, the ``EnzoCenteredFieldRegistry``
has the potential to be a general purpose tool that can be used for
other purposes in the Enzo layer of the codebase.

The idea is to maintain a list of all quantities, represented by
cell-centered fields, that are used by Enzo in the ``FIELD_TABLE``
macro. For each quantity, the table currently tracks its name, whether
its fundamentally a scalar or vector quantity, if the quantity can be
classified as "conserved", "specific", or "other", and whether or not
there is any circumstance in the codebase where the quantity is
"actively" advected. An entry about a scalar quantity registers a
field of the same name. An entry for a vector quantity registers 3
fields with the names ``"{qname}_x"``, ``"{qname}_y"``,
``"{qname}_z"``, where ``{qname}`` is the name of the quantity (e.g.
the row for the "velocity" quantity registers the ``"velocity_x"``,
``"velocity_y"``, and ``"velocity_z"`` fields).

At present, the registry currently provides operations:

  * to access quantity properties registered in ``FIELD_TABLE`` at
    runtime
  * to provide a list of known groups that can be used in the input file
    to identify fields as passively advected scalars (as of now, the
    only such group is ``"color"``).

==============
General Design
==============

    .. _GeneralDesignOverview-section:

Overview
--------

The hydrodynamic/MHD C++ toolkit can be summarized as a series of
classes that encapsulate various operations that performed by
hydrodynamic/MHD integrators. In most cases an abstract base class
exists to provide the interface for each operation. The main operation
classes include:

  * ``EnzoEquationOfState`` - encapsulates many of the operations
    related to the fluid's equation of state (e.g. computing pressure,
    converting the integration quantities or primitives)

  * ``EnzoReconstructor`` - encapsulates interpolation algorithms to
    reconstruct left/right interface states of cell-centered values

  * ``EnzoRiemann`` - encapsulates various Rimann Solver algorithms

  * ``EnzoIntegrationQuanUpdate`` - encapsulates the operation of
    updating integration quantities after a (partial) time-step.

  * ``EnzoBfieldMethod`` - encapsulates operations related to integrating
    magnetic fields that are not performed by the other operation classes.
    For example, a subclass exists for supporting Constrained Transport.

Each of these operation classes are fairly modular (to allow for
selective usage of the frame work components). However, all of the
classes require that an instance of ``EnzoEquationOfState`` get's
passed. The operation classes are also provided with ``PUP`` methods
to allow for easy serialization alongside the ``Method`` class that
makes use of them.

Each of the operation classes are designed to be configured upon
initialization. The instances can then be used multiple times per
time-step (along multiple dimensions if the operation is directional)
and in other time-steps. Lists (excluding passive scalars) of the
expected primitives and integration keys are respectively
*registered* during the construction of ``EnzoReconstructor`` and
``EnzoIntegrationQuanUpdate``. These keys must each share a name
with the registered quantities in ``FIELD_TABLE``. In contrast,
configuration of ``EnzoRiemann``, is less flexible and instances
actually specify the non-passive integration quantities and
non-passive primitives that they require. This difference exists
because the operations encapsulated by ``EnzoReconstructor`` and
``EnzoIntegrationQuanUpdate`` can be applied to individual quantities
in a far more independent manner.

Because all fields storing passively advected scalars are not
necessarily known when initializing a hydro/MHD integrator (i.e.
they could be initialized by a different Method or an initializer),
the passively advected scalars don't need to be registered when
constructing these classes. Instead, a ``std::vector<std::string>``
specifying the names of the passive scalars is often passed to the
method(s) of the class that perform(s) the encapsulated operation.

The implementation of these operation classes aims to avoid the
traditional approach in which field data is directly accessed from
a large array using macros or globally defined unscoped enums that
maps quantity component names to indices. This traditional approach
makes the introuduction of optional fields that are related to active
advection somewhat difficult (e.g. cosmic ray energy/fluxes, internal
energy for dual energy formalism, phi for dedner divergence cleaning).
Instead, our toolkit largely operates on maps/dictionaries containing
``EFlt3DArray`` instances (stored in ``EnzoEFltArrayMap``).

Use of ``EnzoEFltArrayMap``
---------------------------

The basic unit that get's operated on by these operation classes are
instances of the ``EnzoEFltArrayMap`` class. As the name may suggest,
these classes serve as a map/dictionary of instances of
``EFlt3DArray`` (or equivalently, instances of
``CelloArray<enzo_float,3>``). For more details about how to use
``EnzoEFltArrayMap``, see :ref:`EnzoEFltArrayMap-Description`

In the context of this toolkit, the keys of an ``EnzoEFltArrayMap``
are usually the names of a scalar quantity (like ``"density"``) or
component of a vector quantity (like ``"velocity_x"``). Each key is
paired with an instance of ``EFlt3DArray`` that stores associated
data. To simplify logic, arrays are not aliased between separate maps.
Below, we provide a description of the main uses of
``EnzoEFltArrayMap`` by the provided operation classes:

  1. Map of cell-centered integration quantities.

     * This has keys named for all integration scalar quantities and
       components of integration vector quantities. The associated
       arrays hold the values of the cell-centered quantities at a
       given time.

     * This also contains key-value pairs for passively advected
       scalars. In this context, the passive scalars are stored in
       "conserved" form.

     * In a predictor-corrector scheme (like VL+CT), we might have
       multiple maps used to store values at different partial
       timesteps.

  2. Map of cell-centered primitive quantities.

     * This map is used to temporarily store the cell-centered
       primitive quantities for use in reconstruction.

     * This also contains key-value pairs for passively advected
       scalars. In this context, the passive scalars are stored in
       "specific" form.

     * Quantities in both the primitive map and integration map should
       NOT be aliases of each other. They should be deepcopies instead.

  3. Map of temporary cell-centered values for tracking the total
     change in a quantity over a timestep.

     * This map holds key-array pairs named for all integration
       quantities. For each (partial) timestep, these arrays are used
       to accumulate the total change in the conserved form of each
       quantity. This includes the flux divergence and the
       contributions from source terms. At the end of the (partial)
       timestep, these are used to actually update the values of the
       integration quantities

  4. Map of reconstructed left/right primitive quantites

     * 2 instances of ``EnzoEFltArrayMap`` are used to respectively
       hold the reconstructed left and right interface primitive
       quantities. This should share have the same keys that are
       described for the second category of maps.
     * These maps are frequently passed to instances of
       ``EnzoReconstructor`` to store the reconstructed passively
       advected scalars and primitive quantities. Then, these are
       frequently passed to ``EnzoRiemann`` to compute fluxes for
       the integration quantities and passively advected scalars.

  5. Maps of Riemann Flux fields

     * An instance of this kind of map is required for each
       dimension and is used to hold the face-centered fluxes along
       that dimension. The contained arrays should all be defined with the
       appropriate shape for holding data stored on the mesh face along the
       dimension corresponding to the flux. In other words, if a block
       normally holds ``n`` elements (including ghost zones) along axis
       ``i``, then an array used to store fluxes along axis ``i`` should
       hold ``n-1`` elements along axis ``i``.
     * This should have all of the same keys that are in the the first
       category of maps.
     * This kind of map should contain keys named for all passively advected
       scalars and registered integration quantities. The set of keys in these
       maps should be identical to the set of keys in the first category of
       maps, regardless of whether a quantity is "specific" or "conserved"
       (e.g. the map will hold a "velocity_x" key even though the associated
       array stores the x-component of the momentum density flux).

In general, the use of ``EnzoEFltArrayMap`` objects with common sets
of keys helps simplify the implementation of various methods (e.g. the
cell-centered array associated with "density" is used to reconstruct
values that are stored in the fields of the "density"
array in the primitive map).


=================
Equation Of State
=================

All of the operations related to the equation of state are handled by
subclasses of the abstract base class, ``EnzoEquationOfState``. The
class has a number of responsibilities. Currently the only concrete
subclass of ``EnzoEquationOfState`` is the ``EnzoEOSIdeal`` class
which encapsulates the properties of an ideal, adiabatic gas. This
class can optionally support use of the dual-energy formalism (For
details about the currently expected implementation of the
dual-energy formalism see the ``"modern"`` section of
:ref:`using-fluid_props-de` ).

The ``EnzoEquationOfState`` has the following interface:

.. code-block:: c++

   bool is_barotropic();

Returns whether the equation of state is barotropic or not.

*Currently, no barotropic equations of state have been implemented and
none of the wavespeed calculations for the Riemann solvers currently
support barotropic equations of state.*

.. code-block:: c++

   bool uses_dual_energy_formalism();

Returns whether the dual energy formalism is in use.

.. code-block:: c++

   enzo_float get_gamma();

Returns the ratio of the specific heats. This is only required to
yield a reasonable value if the gas is not barotropic.

.. code-block:: c++

   enzo_float get_isothermal_sound_speed();

Returns the isothermal sound speed. This is only required to yield a
reasonable value for barotropic equations of state.

.. code-block:: c++

   enzo_float get_density_floor();

Returns the density floor.

.. code-block:: c++

   enzo_float get_pressure_floor();

Returns the thermal pressure floor.

.. code-block:: c++

   apply_floor_to_energy_and_sync(EnzoEFltArrayMap &integration_map,
                                  int stale_depth);

This method applies the applies the pressure floor to the total_energy
array specified in ``integration_map``. If using the dual-energy formalism
the floor is also applied to the internal energy (also specified in 
``integration_map``) and synchronizes the internal energy with the total
energy. If the equation of state is barotropic, this should do nothing.

.. code-block:: c++

   void pressure_from_integration(const EnzoEFltArrayMap &integration_map,
                                  const CelloArray<enzo_float, 3> &pressure,
                                  int stale_depth,
                                  bool ignore_grackle = false);

This method computes the pressure from the integration quantities
(stored in ``integration_map``) and stores the result in ``pressure``.
This wraps the ``EnzoComputePressure`` object whose default behavior
is to use the Grackle-supplied routine for computing pressure when the
simulation is configured to use ``EnzoMethodGrackle``. The
``ignore_grackle`` parameter can be used to avoid using that routine (the
parameter is meaningless if the Grackle routine would not otherwise
get used). This parameter's primary purpose is to provide the option
to suppress the effects of molecular hydrogen on the adiabatic index
(when Grackle is configured with ``primordial_chemistry > 1``).

.. code-block:: c++

   void primitive_from_integration
     (const EnzoEFltArrayMap &integration_map, EnzoEFltArrayMap &primitive_map,
      int stale_depth, const std::vector<std::string> &passive_list,
      bool ignore_grackle = false);

This method is responsible for computing the primitive quantities (to
be held in ``primitive_map``) from the integration quantities (stored
in ``integration_map``).  Non-passive scalar quantities appearing in
both ``integration_map`` and ``primitive_map`` are simply deepcopied
and passive scalar quantities are converted from conserved-form to
specific form. For a non-barotropic EOS, this also computes pressure
(by calling ``EnzoEquationOfState::pressure_from_integration``).

*In the future, it might be worth considering making this into a subclass
of Cello's ``Physics`` class. If that is done, it may be advisable to
allow for switching between different dual-energy formalism
implementations.*


How to extend
-------------

New equations of state can be added by subclassing and providing the
subclass with implementations for the pure virtual functions
``EnzoEquationOfState``. *Once a second concrete subclass of*
``EnzoEquationOfState`` *is provided, it may be worthwhile to introduce
a factory method.*

=============
Reconstructor
=============

The reconstruction algorithms have been factored out to their own
classes. All implementation of reconstruction algorithms are derived
from the ``EnzoReconstructor`` abstract base class.

To get a pointer to an instance of a concrete implementation of
``EnzoReconstructor``, use the
``EnzoReconstructor::construct_reconstructor`` static factory method:

.. code-block:: c++

   EnzoReconstructor* construct_reconstructor
    (const std::vector<std::string> active_primitive_keys,
     std::string name, enzo_float theta_limiter);

The factory method requires that we register the keys of the
non-passive scalar primitive quantities that are are to be
reconstructed via ``active_primitive_keys``. We specify
the name of the reconstruction algorithm, ``name``. Note that the
primitive keys should correspond to quantities specified in
``FIELD_TABLE`` ; for more details about ``FIELD_TABLE``, see
:ref:`Centered-Field-Registry`

Public Interface
----------------
The main interface function provided by this class is:

.. code-block:: c++

    void reconstruct_interface
      (const EnzoEFltArrayMap &prim_map, EnzoEFltArrayMap &priml_map,
       EnzoEFltArrayMap &primr_map, int dim, EnzoEquationOfState *eos,
       int stale_depth, const std::vector<std::string>& passive_list);

This function takes the cell-centered primtive quantities (specified
by the contents of ``prim_map``) and computes the left and right
reconstructed states (the results are stored in ``priml_map`` and
``primr_map``) along the dimension specifed by ``dim``. If dim has a
value of ``0``/ ``1``/ ``2`` then the values are reconstructed along
the x-/y-/z-axis. ``stale_depth`` indicates the current stale_depth
for the supplied cell-centered quantities (prior to
reconstruction). ``priml_map`` and ``primr_map`` should have the same
shapes as ``prim_map``, except along the reconstruction axis; along that
axis ``prim_map`` should be able to hold 1 more value.
``passive_list`` is used to specify the
names (keys) of the passively advected quantities that are to be
reconstructed.

The ``int EnzoReconstructor::immediate_staling_rate()`` method is
provided to determine the amount by which the stale depth increases
immediately after reconstruction, for a given algorithm. The
``int EnzoReconstructor::delayed_staling_rate()`` method returns how much
the stale depth increases after adding flux divergence, computed from
the reconstructed values, to the integration quantities  (this is
normally 1). Finally ``int EnzoReconstructor::total_staling_rate()``
gives the sum of the results yielded by the prior 2 methods.

How to extend
-------------

To add a new reconstructor, subclass ``EnzoReconstructor`` and provide
definitions for the virtual methods.  The implementations of the
``immediate_staling_rate()`` and ``total_staling_rate()`` virtual
methods must also be provided. Additionally, the factory method
``EnzoReconstructor::construct_reconstructor`` must also be modified
to return pointers to instances of the new class when the appropriate
name is passed as an argument, and the name of the new reconstructor
should be added to :ref:`using-vlct-reconstruction`

Currently, to add new slope limiters for existing reconstruction
algorithms new classes are effectively defined. The piecewise linear
reconstruction algorithm is implemented as a class template
``EnzoReconstructorPLM<Limiter>`` where ``Limiter`` is a functor that
implements a specific slope limiter. ``Limiter`` must be default
constructible and provide a function call operation, `operator()`. The
function call operation must have a signature matching:

.. code-block:: c++

   enzo_float Functor::operator()(enzo_float vm1, enzo_float v, enzo_float vp1,
                                  enzo_float theta_limiter);

Give three contiguous primitive values along the axis of
interpolation, (``vm1``, ``v``, and ``vp1``) the method should compute the
limited slope. The ``theta_limiter`` parameter that can be optionally
used to tune the limiter (or ignored).

When a new a ``Limiter`` functor is defined to be used to specialize
``EnzoReconstructorPLM``, the new specialization must be added to
enzo.CI. The other steps mentioned at the start of this subsection for
implementing new reconstruction algorithms must also be followed.

*The use an enum with a switch statement was considered for switching
between different slope limiters. However we determined that the compiler
would not pull the switch statement outside of the loop.
Therefore templates are used to avoid executing the switch statement on
every single iteration.*

*Having multiple slope limiters available at runtime may be
unnecessary (or not worth the larger binary size). It might be worth
considering using preprocessor macros to allow for specification of
the slope limiter at compile time.*

==============
Riemann Solver
==============

All implementations of (approximate) Riemann solver algorithms are
derived from the ``EnzoRiemann`` abstract base class.

    .. _Riemann-Usage-section:

Usage Notes
-----------

To get a pointer to a concrete implemenation of ``EnzoRiemann``, call the
static factory method:

.. code-block:: c++

   EnzoRiemann* EnzoRiemann::construct_riemann
   (const EnzoRiemann::FactoryArgs& factory_args) noexcept;

The above signature looks a little intimidating. In reality, you could write
something more like the following snippet

.. code-block:: c++

   EnzoRiemann* ptr = EnzoRiemann::construct_riemann({solver, mhd,
                                                      internal_energy});

in which
  * ``solver`` is a ``std::string`` specifying the name of the solver
  * ``mhd`` is a boolean specifying whether magnetic fields are present
  * ``internal_energy`` is a boolean specifying if the internal energy flux
    must be computed.

This weird indirection only currently exists to accomodate ``charm++``\'s
``pup`` functionality. Once Enzo-E transitions to its custom restart
functionallity, we'll remove this indirection.

An instance of ``EnzoRiemann`` specifies the expected non-passive keys
(and key-order) that the ``flux_map`` argument should have when passed to its
``solve`` method (these keys correspond to integration quantities).

.. code-block:: c++

   const std::vector<std::string> integration_quantity_keys() const;

The following method specifies the expected non-passive keys (and key-order)
that the ``priml_map`` and ``primr_map`` arguments should have when passed
an ``EnzoRiemann``\'s ``solve`` method (these keys correspond to primitive
quantities).

.. code-block:: c++

   const std::vector<std::string> primitive_quantity_keys() const;


The main interface function of ``EnzoRiemann`` is:

.. code-block:: c++

   void solve(const EnzoEFltArrayMap &prim_map_l,
              const EnzoEFltArrayMap &prim_map_r,
              EnzoEFltArrayMap &flux_map, int dim, EnzoEquationOfState *eos,
              int stale_depth, const str_vec_t &passive_list,
              const CelloArray<enzo_float,3> * const interface_velocity) const;

In this function, the ``prim_map_l`` and ``prim_map_r`` arguments are
references to the ``EnzoEFltArrayMap`` objects holding the arrays of
reconstructed left/right primitive quantities. The ``flux_map``
argument holds the face-centered arrays where the computed fluxes for
each integration quantity are written. ``dim`` indicates the dimension
along which the flux should be computed (0,1,2 corresponds to x,y,z).
``interface_velocity`` is an optional argument used to specify a
pointer to an array that can be used to store interface velocity
values computed by the Riemann Solver (this is primarily used for
computing internal energy source terms when the dual energy formalism
is in use).

Some additional notes:

  *  The first ``EnzoRiemann::primitive_quantity_keys().size()`` keys of
     ``prim_map_l`` and ``prim_map_r`` should match the values and order of
     ``EnzoRiemann::primitive_quantity_keys()``.

  * Likewise, the first ``EnzoRiemann::integration_quantity_keys().size()``
    keys of ``flux_map`` should match the values and order of
    ``EnzoRiemann::integration_quantity_keys()``.

  * ``prim_map_l``, ``prim_map_r``, and ``flux_map`` should also each contain
    keys for each of the passive scalars in ``passive_list`` (the order of
    these is not currently enforced).

  * All of the arrays in ``prim_map_l``, ``prim_map_r``, and ``flux_map``
    should have the same shape. If ``interface_velocity`` is specified, it
    should also have that shape.

  * Calling the ``contiguous_arrays()`` instance method for ``prim_map_l``,
    ``prim_map_r``, and ``flux_map`` must return ``true`` in each case (in
    other words, each array in the array maps should be stored in
    nearly-contiguous 4D arrays).

Implementation Notes: Kernels
-----------------------------

At the time of writing, the calculations specific to each Riemann
Solver are implemented as *kernels* (inspired by *Kokkos*). More
precise requiements are detailed in :ref:`KernelReq-section`, but in
short each kernel:

  * defines a ``operator()(int iz, int iy, int ix) const`` method for
    computing the flux at the specified cell-interface

  * specifies a specialization of the template class
    ``EnzoRiemannLUT<InputLUT>``. This both indicates the expected
    actively advected integration (and primitive) quantities **AND**
    acts as a compile-time lookup table (that maps quantity names to
    unique array indices). See :ref:`EnzoRiemannLUT-section` for
    more details.

  * Specifies the expected equation of state, by specifying the expected
    type of EOS Struct objects. See :ref:`EOSStructObject-section` for
    more details about EOS Struct objects.

  * are configured by an instance of ``KernelConfig<EOSStructT>`` (see
    :ref:`KernelConfig-section` for more details).

The ``EnzoRiemannImpl<Kernel>`` class template is used to wrap each
kernel. This class template subclasses the ``EnzoRiemann`` abstract
base class and is used to implement the interface described above,
in :ref:`Riemann-Usage-section`, for each kernel (i.e. in other words,
we use it to implement "type erasure"). In more detail,
``EnzoRiemannImpl<Kernel>::integration_quantity_keys()`` and
``EnzoRiemannImpl<Kernel>::primitive_quantity_keys()`` methods specify
the fields (and field ordering) required by a given kernel.

The steps of ``EnzoRiemannImpl<Kernel>::solve`` are fairly
straight-forward:

  1. An instance of ``KernelConfig<Kernel::EOSStructT>`` is created.

  2. The ``Kernel`` is constructed and executed at each cell-interface

  3. The fluxes for passively advected scalars are computed
     (this step is completely independent of the choice of ``Kernel``).


    .. _EnzoRiemannLUT-section:

EnzoRiemannLUT
~~~~~~~~~~~~~~

As described above in the :ref:`GeneralDesignOverview-section`\, we
sought to avoid the common approach of
hydro codes that map actively advected quantities indices with macros
or globally defined unscoped enums. The ``EnzoRiemannLUT<InputLUT>``
template class basically serves as a compromise between this traditional
approach approach and using a hash table (which introduce unacceptable
overhead) for organizing quantities in the main loop of
``EnzoRiemannImpl<Kernel>``. Alternatively it can be thought of as a
scoped version of the traditional approach.

This is a template class that provides the following features at compile
time:

    * a lookup table (LUT) that maps the names of components of a subset
      of the actively advected integration quantities defined in
      ``FIELD_TABLE`` to unique, contiguous indices.

    * the number of integration quantity components included in the table

    * a way to iterate over just the conserved or specific integration
      quantities values that are stored in an array using these mapping

    * a way to query which of the actively advected integration quantities
      in FIELD_TABLE are not included in the LUT

These feature are provided via the definition of publicly accessible
integer constants in every specialization of the template class. All
specializations have:

    * a constant called ``num_entries`` equal to the number of integration
      quantity components included in the lookup table

    * a constant called ``specific_start`` equal to the number of components
      of conserved integration quantities included in the lookup table

    * ``qkey`` constants, which include constants named for the components
      of ALL actively advected integration quantities in FIELD_TABLE. A
      constant associated with a SCALAR quantity, ``{qname}``, is simply
      called ``{qname}`` while constants associated with a vector quantity
      ``{qname}`` are called ``{qname}_i``, ``{qname}_j``, and ``{qname}_k``.

The ``qkey`` constants serve as both the keys of the lookup table and a
way to check whether a component of an actively advected quantity is
included in the table. Their values are satisfy the following conditions:

    * All constants named for values corresponding to quantities NOT
      included in the lookup table have values of ``-1``

    * All constants named for conserved integration quantities have unique
      integer values in the internal ``[0,specific_start)``

    * All constants named for specific integration quantities have unique
      integer values in the interval ``[specific_start, num_entries)``

    * Constants must be defined for all 3 components (or none of the
      components) of a vector quantity (e.g. velocity or magnetic
      fields).  Additionally, the ``k``\th component of a vector
      quantity is expected to have a value that is 1 larger than that
      of the ``j``\th component and 2 larger than the ``i``\th
      component.

The lookup table is always expected to include density and the 3 velocity
components.

This template class also provides a handful of helpful static methods to
programmatically probe the table's contents at runtime and validate that
the above requirements are specified.

For the sake of providing some concrete examples about how the code works,
let's assume that we have a class ``MyIntegLUT`` that is defined as:

.. code-block:: c++

   struct MyIntegLUT {
     enum vals { density=0, velocity_i, velocity_j, velocity_k,
                 total_energy, num_entries, specific_start = 1};
   };

The template specialization ``EnzoRiemannLUT<MyIntegLUT>`` assumes that
all undefined ``qkey`` constants omitted from ``MyIntegLUT`` are not included
in the lookup table and will define them within the template specialization
to have values of ``-1``.

To access the index associated with density or the jth component of
velocity, one would evaluate:

.. code-block:: c++

   int density_index = EnzoRiemannLUT<MyIntegLUT>::density; //=0
   int vj_index = EnzoRiemannLUT<MyIntegLUT>::velocity_j;   //=2

Additionally, the value of ``EnzoRiemannLUT<MyIntegLUT>::bfield_k`` would be
``-1``.

It makes more sense to talk about the use of this template class when we
have a companion array. For convenience, the alias template
``lutarray<LUT>`` type is defined. The type,
``lutarray<EnzoRiemannLUT<InputLUT>>`` is an alias of the type
``std::array<enzo_float, EnzoRiemannLUT<InputLUT>::num_entries>;``.

As an example, imagine that the total kinetic energy density needs to be
computed at a single location from an values stored in an array, ``integ``,
of type ``lutarray<EnzoRiemannLUT<MyIntegLUT>>``:

.. code-block:: c++

   using LUT = EnzoRiemannLUT<MyIntegLUT>;
   enzo_float v2 = (integ[LUT::velocity_i] * integ[LUT::velocity_i] +
                    integ[LUT::velocity_j] * integ[LUT::velocity_j] +
                    integ[LUT::velocity_k] * integ[LUT::velocity_k]);
   enzo_float kinetic = 0.5 * integ[LUT::density] * v2;


``EnzoRiemannLUT<InputLUT>``, makes it very easy to
write generic code that can be reused for multiple different lookup table
by using by passing its concrete specializations as a template argument
to other template functions/classes. Consider the case where a single
template function is desired to compute the total non-thermal energy
density at a single location for an arbitrary lookup table:

.. code-block:: c++

   template <class LUT>
   enzo_float calc_nonthermal_edens(lutarray<LUT> prim)
   {
     enzo_float v2 = (prim[LUT::velocity_i] * prim[LUT::velocity_i] +
     prim[LUT::velocity_j] * prim[LUT::velocity_j] +
     prim[LUT::velocity_k] * prim[LUT::velocity_k]);

     enzo_float bi = (LUT::bfield_i >= 0) ? prim[LUT::bfield_i] : 0;
     enzo_float bj = (LUT::bfield_j >= 0) ? prim[LUT::bfield_j] : 0;
     enzo_float bk = (LUT::bfield_k >= 0) ? prim[LUT::bfield_k] : 0;
     enzo_float b2 = bi*bi + bj*bj + bk*bk;

     return 0.5(v2*prim[LUT::density] + b2);
   }


.. _EOSStructObject-section:

EOSStruct Objects
~~~~~~~~~~~~~~~~~~

These objects are supposed to be lightweight struct/classes that
encapsulate an equation of state.  It's also important that these
objects are cheap to copy. It's our intention to keep these
independent of the other Riemann Solver tooling (particularly
``EnzoRiemannLUT``).

Our only concrete example at this time is ``EOSStructIdeal``. But
this very much expresses the basic idea. The struct provides methods
that just require density and pressure values to compute:

  * specific internal energy
  * internal energy density
  * sound speed
  * fast magnetosonic speed (this requires magnetic field values to
    also be specified)

In the future, we may need to slightly revisit the expected signature
if/when we add support for an isothermal fluid.

.. note::

   At this time, ``EOSStructIdeal`` is completely independent of the
   ``EnzoEOSIdeal`` subclass of ``EnzoEquationOfState``.

   In the immediate future, there are plans to unify these
   implementations (the ``EnzoEquationOfState`` machinery will be
   refactored to make use of these ``EOSStruct`` objects). When that
   comes to pass, we will flesh out this section in greater detail.


    .. _KernelConfig-section:

``KernelConfig<EOSStructT>``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``KernelConfig`` template class simply groups the configuration
parameters for kernel (the alternative would be to require that the
members of ``KernelConfig`` are listed as members for every class
implementing a kernel).

The relationship between this class and a kernel is modelled after the
relationship between captured values and a lambda function. We made
this class as lightweight as possible to encourage the compiler to
make similar optimizations. The declaration of ``KernelConfig`` is
reproduced below:

.. literalinclude:: ../../../src/Enzo/EnzoRiemann/EnzoRiemannImpl.hpp
   :language: c++
   :start-after: SPHINX-SNIPPET-KERNELCONFIG-START-INCLUDE
   :end-before: SPHINX-SNIPPET-KERNELCONFIG-END-INCLUDE

.. note::
   In the near-future, the above snippet will be replaced with nicer looking
   documentation that will be auto-generated using doxygen and breathe
   
More details will be given below about internal energy calculations.

    .. _KernelReq-section:

Kernel Requirements
~~~~~~~~~~~~~~~~~~~

The basic skeleton for defining a kernel is defined below. For concreteness, we've assumed that this kernel uses the ``MHDLUT`` and an ideal EOS (but those could easily be changed).

.. code-block:: c++

   struct MyKernel
   {

   public: // typedefs
     using LUT = EnzoRiemannLUT<MHDLUT>;
     using EOSStructT = EOSStructIdeal;

   public: // fields
     const KernelConfig<EOSStructT> config;

   public: // methods
     FORCE_INLINE void operator()(const int iz,
                                  const int iy,
                                  const int ix) const noexcept
     {
       // compute the fluxes at (iz, iy, ix) & store result in config.flux_arr
       //
       // For concreteness:
       //   - the left and right reconstructed density values (for the same
       //     cell-interface) can be retrieved using:
       //       config.prim_arr_l(LUT::density, iz, iy, ix)
       //       config.prim_arr_r(LUT::density, iz, iy, ix)
       //   - we want to write computed density flux to:
       //       config.flux_arr(LUT::density, iz, iy, ix)
       //
       // It may be helpful to point out that:
       //   when (config.dim == 0) -> this should computes on x-faces
       //   when (config.dim == 1) -> this should computes on y-faces
       //   when (config.dim == 2) -> this should computes on z-faces
       //
       // For more information, including how exactly indices map to spatial
       // locations, check the docstrings for KernelConfig. In practice, it's
       // sufficient to understand that (iz,iy,ix) refers to the same spatial
       // location in ALL of the arrays.
       //
       // assuming it makes sense for the EOS (it's non-barotropic), this should
       // also compute values needed for the dual energy formalism and store
       // results in config.internal_energy_flux_arr(iz,iy,ix) &
       // config.velocity_i_bar_arr(iz,iy,ix)

     }
   };

There are a couple of things to note:

- As explained in :ref:`KernelConfig-section`, when you access vector
  quantities from ``config.prim_arr_l``, ``config.prim_arr_r``, or
  ``config.flux_arr`` using ``LUT`` (e.g. ``LUT::velocity_i``,
  ``LUT::velocity_j``, ``LUT::bfield_k``), the ``i``, ``j``, & ``k``
  components **always** map to the ``x``, ``y``, & ``z`` components,
  respectively.

    - To be concrete: ``config.flux_arr(LUT::velocity_j, iz, iy, ix)``
      always refers to the flux of the ``y`` velocity component.

    - To try preserve symmetry, we often use ``config.dim`` to remap
      the ``i``, ``j``, and ``k`` components (like what is done in
      Athena++). The current convention is to include the following
      code block at the start of ``operator()``:

      .. code-block:: c++

        const int external_velocity_i = config.dim + LUT::velocity_i;
        const int external_velocity_j = ((config.dim+1)%3) + LUT::velocity_i;
        const int external_velocity_k = ((config.dim+2)%3) + LUT::velocity_i;
        const int external_bfield_i = config.dim + LUT::bfield_i;
        const int external_bfield_j = ((config.dim+1)%3) + LUT::bfield_i;
        const int external_bfield_k = ((config.dim+2)%3) + LUT::bfield_i;

      and to use ``external_velocity_i``, ``external_velocity_j``,
      ``external_velocity_k`` instead of ``LUT::velocity_i``,
      ``LUT::velocity_j``, ``LUT::velocity_k`` when accessing values
      from ``prim_arr_l``, ``prim_arr_r``, and ``flux_arr``.

   - We have given some thought to trying to abstract this vector
     permutation, but it unfortunately remains unclear how to do
     this in a way that preserves performance across different
     compilers, without requiring template specializations of
     ``operator()`` for computing fluxes along different
     dimensions.

- It's also possible to make ``config`` a private field or to
  introduce additional fields. However, if you do either of those
  things, you will need to define a constructor that accepts
  ``config`` as a single argument.

- ``FORCE_INLINE`` is a macro that uses compiler-specific extensions to force
  inlining of functions. Please use this sparingly outside, only in cases where
  you see a noticable speedup (inlining everything can actually slow code down
  from inflating the binary's size).


Dual Energy Formalism Treatment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dual-energy formalism requires that a Riemann Solver computes fluxes for the internal energy, and a velocity component at the cell-interfaces.

Computing these quantities doesn't change the Riemann Solver
calculation (or the required reconstructed primitives), it just
involves an extra calculation. To avoid unnecessary template code generation
or introducing branching:

- ``EnzoRiemannLUT`` generally does not include ``internal_energy`` as an
  entry (there are other reasons why this is convenient)

- Kernels **always** compute the necessary quantities for the dual-energy
  formalism (assuming the EOS is compatible with the dual-energy formalism)

Instead, ``EnzoRiemannImpl<Kernel>`` manages the dual-energy configuration:

- ``EnzoRiemannImpl<Kernel>::solve`` will **ALWAYS** make sure that
  ``config.internal_energy_flux_arr`` and
  ``config.velocity_i_bar_arr`` are valid arrays. If arrays aren't
  provided (i.e. the dual-energy formalism is not in use), scratch
  data will be used for this explicit purpose.

- ``EnzoRiemannImpl<Kernel>`` is responsible for including
  ``"internal_energy"`` in the output of
  ``integration_quantity_keys()``, based on how it's configured.


Adding new quantites
--------------------

To add support for new actively advected integration cell-centered
quantities (e.g. cosmic ray energy/flux), the table of cell-centered
quantities (``FIELD_TABLE``) must be updated. See
:ref:`Centered-Field-Registry` for more details.  To add support for
computing fluxes for such quantities, modifications must be made to
either ``EnzoRiemannImpl`` or the ``Kernel`` of an existing
solver. Alternatively, for certain quantities, a brand new solver
may need to be introduced.

When adding a new integration vector quantity, you also need to add a
few lines to the main for-loop of ``EnzoRiemannImpl`` for copying
values to ``wl``/``wr`` and from ``fluxes`` (The existing code doing
this for the velocity and magnetic fields should be used as a guide).

Adding new solvers
------------------

New Riemann Solvers can currently be added to the infrastructure by
either subclassing ``EnzoRiemann`` or defining a new specialization
of ``EnzoRiemannImpl<Kernel>``. In either case, the
``EnzoRiemann::construct_riemann`` factory method must be modified to
return the new solver and :ref:`using-vlct-riemann-solver`
should be updated.

The additional steps for implementing a new Riemann solver by speciallizing
``EnzoRiemannImpl<Kernel>`` are as follows:

  1. Define a new ``Kernel`` class (e.g. ``HLLDKernel``)

  2. *(optional)* define an alias name for the specialization of
     ``EnzoRiemannImpl`` that uses the new ``Kernel`` class
     (e.g. ``using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDKernel>;``).


.. note::

   Due to the template-heavy nature of our implementation, the
   Riemann Solvers been separated from the rest of the Enzo-layer into
   their own sub-library.

   This choice was made because the convention in the Enzo-layer (at
   least before we made this separation) is to effectively include
   every header file in every source file. Since template code usually
   needs to be written in header files, changes to the Riemann Solver
   algorithms used to frequently trigger rebuilds of the entire
   Enzo-layer. This separation has significantly sped up incremental
   rebuilds during the development workflow for the Riemann Solvers.

===============================
Updating integration quantities
===============================

The ``EnzoIntegrationQuanUpdate`` class has been provided to encapsulate
the operation of updating integration quantities after a (partial)
time-step. The operation was factored out of the ``EnzoMethodMHDVlct``
class since it appear in all Godunov solvers.

The constructor for ``EnzoIntegrationQuanUpdate`` has the following
signature:

.. code-block:: c++

   EnzoIntegrationQuanUpdate(std::vector<std::string> integration_quantity_keys,
                             bool skip_B_update)

The function requires that we:

  * register the keys of the integration quantities (with
    ``integration_quantity_keys``)
  * indicate whether the update to the magnetic field should
    be skipped.

The integration quantity keys should match the names specified
in ``FIELD_TABLE``; see :ref:`Centered-Field-Registry` for more
details. The update to the magnetic field should be skipped when
Constrained Transport is in use (since the magnetic field update is
handled separately). If the magnetic field is not specified as an
integration quantity, then the value specified for ``skip_B_update`` is
unimportant

The following method is used to compute the change in (the conserved
form of) the integration and passively advected quantites due to the
flux divergence along dimension ``dim`` over the (partial) timestep
``dt``. The arrays in ``dUcons_map`` are used to accumulate the total
changes in these quantities. ``passive_list`` lists the names (keys)
of the passively advected scalars.

.. code-block:: c++

   void accumulate_flux_component
     (int dim, double dt, enzo_float cell_width,
     const EnzoEFltArrayMap &flux_map,
      EnzoEFltArrayMap &dUcons_map, int stale_depth,
      const std::vector<std::string> &passive_list) const;

The method used to clear the values of the arrays used for accumulation is
provided below. This sanitization should be performed before starting
to accumulate flux divergence or source terms. The ``passive_list``
argument is used in the same way as the previous function.

.. code-block:: c++

    void clear_dUcons_group(EnzoEFltArrayMap &dUcons_map, enzo_float value,
                            const std::vector<std::string> &passive_list) const;

The method used to actually add the accumulated change in the integration
(specified in ``dUcons_map``) to the values of the
integration quantities from the start of the timestep (specificed by
``initial_integration_map``) has the following signature:

.. code-block:: c++

   void update_quantities
     (EnzoEFltArrayMap &initial_integration_map,
      const EnzoEFltArrayMap &dUcons_map,
      EnzoEFltArrayMap &out_integration_map,
      EnzoEFltArrayMap &out_conserved_passive_scalar,
      EnzoEquationOfState *eos, int stale_depth,
      const std::vector<std::string> &passive_list) const;

The fields included in ``dUcons_map`` should include contributions
from both the flux divergence AND source terms. The results for the
actively advected quanties are stored in ``out_integration_map`` and
the results for the passively advected scalars are stored in conserved
form in the arrays held by ``out_conserved_passive_scalar`` (note that
the initial values of the passive scalars specified in
``initial_integration_map`` are in specific form).

==========================
Magnetic Field Integration
==========================

Subclasses of the abstract base class, ``EnzoBfieldMethod`` are used
to implement magnetic field integration-related operations. While
operations like reconstruction and flux calculations of relevant
quantities are expected to be carried out with ``EnzoReconstructor``
and ``EnzoRiemann``, all other magnetic field integration-related
operations should be encapsulated by ``EnzoBfieldMethod``.

Currently, the only subclass is ``EnzoBfieldMethodCT``, which
implements operations related to Constrained Transport. Other
subclasses could be implemented in the future that encapsulate other
integration methods (e.g. divergence cleaning).

From the perspective of an integrator that employs
``EnzoBfieldMethod``, the primary result of each operation is to
modify the values cell-centered/reconstructed quantities, since that's
all the integrator directly needs to know about. In reality, side
effects performed by these operations can be equally as important. For
example, ``EnzoBfieldMethodCT`` implicitly needs to update
face-centered magnetic field values (given that the face-centered
values serve as the primary representation, and the cell-centered
values are derived directly from them).


To accomplish these goals, ``EnzoBfieldMethod``, basically implements
a state machine. It basically provides 3 classes of methods: (i) state
machine-methods, (ii) physics methods, and (iii) descriptor methods.

State Machine Methods
---------------------

When the ``EnzoBfieldMethod`` is first constructed, it has an
uninitialized state. During construction the number of partial
timesteps (``num_partial_timesteps``) involved per cycle must be
specified.

At the beginning of an integration cycle (when an
``EnzoBfieldMethod`` object is unitialized), the cello ``Block`` that
is going to be integrated block that must be specified using the
following method.

.. code-block:: c++

   void register_target_block(Block *block) noexcept;

This method will correctly set the internal state and will invoke the
virtual ``register_target_block_`` method, which is used by subclasses
to preload relevant data from ``block`` and for the delayed
initialization of scratch arrays (since the shapes may not be known at
construction).

Once a target block has been registerred, the ``EnzoBfieldMethod`` object
is now ready to perform integration-related operations for the first partial
timestep (the physics methods can now be called). The following method is
used to increment the partial timestep:

.. code-block:: c++

   void increment_partial_timestep() noexcept;

The target block is unregistered once this method to has been called
``num_partial_timesteps`` times. Any calls to ``register_target_block``
while a block is still registered will currently cause an error.


There are couple of things to keep in mind:

   * Any calls to physics methods or other state machine or other when
     no target block is registered are not allowed.
   * It's EXTREMELY important that ``increment_partial_timestep`` is
     always invoked ``num_partial_timesteps`` after a target block is
     registered and before there is chance for blocks to migrate
     between nodes. In other words, the a target block should always
     be registerred and unregistered during a single call to the cello
     ``Method`` object that represents the integrator.


Physics Methods
---------------

These methods are actually used to perform the relevant magnetic field
integration operations. Each method is a pure virtual method that must
be implemented by a subclass (even if the method immediately
returns). These methods were all written and named based on the
operations of Constrained Transport (CT). In the future, additional
methods may need to be introduced to facillitate the implementation of
other magnetic field integration schemes.

These methods are listed below with brief description. For more details,
please see the docstring. The methods are expected to generally
be called in the general order that they are listed. While this isn't
currently enforced, incorrect results may arise if they aren't called
in the proper order.

In the context of CT, the following method is used to overwrite the
reconstructed value magnetic field component that corresponds to the
axis of reconstruction with the (internally tracked) face-centered
value.

.. code-block:: c++

   void correct_reconstructed_bfield(EnzoEFltArrayMap &l_map,
                                     EnzoEFltArrayMap &r_map, int dim,
                                     int stale_depth) noexcept;

The following method is used by ``EnzoBfieldMethodCT`` to take note of
the upwind direction after computing the Riemann Fluxes along a
dimension ``dim``.

.. code-block:: c++

   void identify_upwind(const EnzoEFltArrayMap &flux_map, int dim,
                        int stale_depth) noexcept;

Finally, the following method is used to actually update the cell-centered
magnetic field values.

.. code-block:: c++

   void update_all_bfield_components
     (EnzoEFltArrayMap &cur_prim_map, const EnzoEFltArrayMap &xflux_map,
      const EnzoEFltArrayMap &yflux_map, const EnzoEFltArrayMap &zflux_map,
      EnzoEFltArrayMap &out_centered_bfield_map, enzo_float dt,
      int stale_depth) noexcept;

In ``EnzoBfieldMethodCT`` this will also update the face-centered
magnetic field values (it assumes that ``identify_upwind`` was called
once for each dimension and uses the stored data). When using this
alongside ``EnzoIntegrationQuanUpdate``, care needs to be taken about the
order in which this method is called relative to
``EnzoIntegrationQuanUpdate::update_quantities`` that accounts for the time
when floors are applied to the total energy.

Descriptor Methods
------------------

These are virtual methods that can be invoked at any time after the
``EnzoBfieldMethod`` object has been constructed. These are used to
describe requirements of the given magnetic field integration method.

Currently, only one such method exists:

.. code-block:: c++

   void check_required_fields() const noexcept;

These may change in the future.

How to extend
-------------

Implementing a new method for magnetic field integration is fairly
straight-forward. Basically all you have to do is implement a subclass
of ``EnzoBfieldMethod``. In addition to providing implementations for
each each physics and descriptor method, the subclass also needs to
implement:

.. code-block:: c++

   void register_target_block_(Block *target_block,
                               bool first_initialization) noexcept;

As mentioned earlier, this method is called by
``register_target_block`` while registering a new target block. In
this call the subclass should preload any data it will need from the
``target_block``. The ``first_initialization`` argument indicate
whether this is the first time a ``target_block`` is being registered
after the instance has been constructed (this includes the first time
following deserialization after a restart). It can be used to help with
lazy intialization of scratch space.

*Once a second concrete subclass of* ``EnzoBfieldMethod`` *is
provided, it may be worthwhile to introduce a factory method.*
