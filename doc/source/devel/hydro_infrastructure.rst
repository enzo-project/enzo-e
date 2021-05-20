****************************
Hydro/MHD C++ Infrastructure
****************************

*[This page is under development]*

In this section we discuss some of the C++ infrastructure provided in
the Enzo Layer that can be optionally used to implement other
hydro/MHD methods. The infrastructure was used to implement other the
VL + CT MHD solver.

*Note: Currently, barotropic equations of state and compatibility with
Grackle (for* ``primordial_chemistry > 0`` *) are not yet implemented
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

Integrable/Reconstructable Quantities
-------------------------------------

Throughout this guide and the relevant sections of the codebase, we
refer to quantities as reconstructable (i.e. they are used for
reconstructing left/right interface states) and integrable (i.e. the
primary quantities that are evolved by the integrator). There is a
high degree of overlap between these categories and the precise
categorization depends on the equation of state. As of now all
integrable quantities are "conserved" or "specific" (a quantity like
velocity that when multiplied by density becomes conserved).

As an example, the categorization of the quantities for an ideal,
adiabatic gas are:

  * density - reconstructable and integrable

  * velocity - reconstructable and integrable

  * pressure - only reconstructable

  * (specific) total energy - only integrable

  * magnetic field - reconstructable and integrable

*Note: Both the reconstructable and integrable quantities are
frequently referred to as primitives throughout the codebase as
primitives. This is mostly historical and mainly refers to the fact
that the collection of quantities are not all of the quantities are
"conserved" (at least some subset of them are "specific").*

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

  * for building ``Grouping`` objects that contain registered quantity
    fields.
  * to access quantity properties registerred in ``FIELD_TABLE`` at
    runtime
  * provide a list of known groups that can be used in the input file
    to identify fields as passively advected scalars (as of now, the
    only such group is ``"colour"``).

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
    converting reconstructable quantities to integrable quantities and
    vice-versa)

  * ``EnzoReconstructor`` - encapsulates interpolation algorithms to
    reconstruct left/right interface states of from cell-centered
    values

  * ``EnzoRiemann`` - encapsulates various Rimann Solver algorithms

  * ``EnzoIntegrableUpdate`` - encapsulates the operation of updating
    integrable quantities after a (partial) time-step.

  * ``EnzoBfieldMethod`` - encapsulates operations related to integrating
    magnetic fields that are not performed by the other operation classes.
    For example, a subclass exists for supporting Constrained Transport.

Each of these operation classes are fairly modular (to allow for
selective usage of the frame work components). However, all of the
classes require that an instance of ``EnzoEquationOfState`` get's
passed.

Each of the operation classes are designed to be configured upon
initialization. The instances can then be used multiple times per
time-step (along multiple dimensions if the operation is directional)
and in other time-steps. The operation classes are also provided with
``PUP`` methods to allow for easy migration alongside the ``Method``
class that makes use them.

For each operation class (other than ``EnzoEquationOfState``), the
expected integrable or reconstructable quantities (other than passively
advected scalars) are *registered* at construction.

  * The names of all reconstructable quantites that get registered
    in the construction of ``EnzoRiemann`` must share a name
    with the registered quantities in ``FIELD_TABLE``.

  * All registered integrable quantity names in the construction of
    ``EnzoRiemann`` or ``EnzoIntegrableUpdate`` must be specified in
    ``FIELD_TABLE`` as quantities that are actively advected in some
    contexts.

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

Overview
~~~~~~~~
The basic unit that get's operated on by these operation classes
are instances of the ``EnzoEFltArrayMap`` class. As the name may
suggest, these classes serve as a map/dictionary of instances of
``EFlt3DArray``.


Specific Usage
~~~~~~~~~~~~~~


In the context of this toolkit, the keys of an ``EnzoEFltArrayMap``
are usually the names of a scalar quantity (like ``"density"``)
or component of a vector quantity (like ``"velocity_x"``). Each key
is paired with an instance of ``EFlt3DArray`` that stores data
associated data. Below, we provide a description of the main
uses of ``EnzoEFltArrayMap`` by the provided operation classes:

  1. Map of cell-centered quantities.

     * This has quantity/quantity-component keys named for all
       integrable and reconstructable quantities used by the
       integrator. The associated arrays hold the values of the
       cell-centered quantities at a given time. We currently store
       integrable and reconstructed quantities together due to the
       high degree of overlap between each category.

     * This also contains key-value pairs for passively advected
       scalars. Note that in this context, the passive scalars are
       usually represented in "specific" form. For reference, the
       general convention throughout Enzo-E is to represent primarily
       passively advected scalars in "conserved" form (as mass
       densities) outside of hydrodynamic integrator methods and to
       convert them to "specific" form (mass fractions) within
       hydrodynamic integrator methods.

  2. Map of temporary cell-centered for tracking the total change
     in a quantity over a timestep.

     * This map holds key-array pairs named for all integrable
       quantities and groups of passively advected scalars. For each
       (partial) timestep, these arrays are used to accumulate the
       total change in the conserved form of each quantity. This
       includes the flux divergence and the contributions from source
       terms. At the end of the (partial) timestep, these are used to
       actually update the values of the integrable quantities

  3. Map of reconstructed left/right quantites

     * 2 instances of ``EnzoEFltArrayMap`` are used to respectively hold
       the reconstructed left and right interface quantities. This should
       share have the same keys that are described for the first category
       of maps.
     * These maps are frequently passed to instances of ``EnzoReconstructor``
       to store the reconstructed passively advected scalars and
       reconstructable quantities. They are then usually passed to
       ``EnzoEquationOfState`` to compute and store the reconstructed
       integrable quantities and reconstructed pressure. Then, these are
       frequently passed to ``EnzoRiemann`` to compute fluxes from the
       integrable quantities and the passively advected scalars.
     * Although this inherently represents data centered on the faces of
       the mesh, the contained arrays should formally have the shape required
       to hold cell-centered data. This is done to facillitate the reuse of
       these maps to hold reconstructed fields along each dimension. This
       means that there is always some unused allocated memory at the end of
       one of the dimensions.

  4. Maps of Riemann Flux fields

     * An instance of this kind of map is required for each
       dimension and is used to hold the face-centered fluxes along
       that dimension. The contained arrays should all be defined with the
       appropriate shape for holding data stored on the mesh face along the
       dimension corresponding to the flux. In other words, if a block
       normally holds ``n`` elements (including ghost zones) along axis
       ``i``, then an array used to store fluxes along axis ``i`` should
       hold ``n-1`` elements along axis ``i``.
     * This kind of map should contain keys named for all passively advected
       scalars and registered integrable quantities. The set of keys in these
       maps should be identical to the set of keys in the first category of
       maps, regardless of whether a quantity is "specific" or "conserved"
       (e.g. the map will hold a "velocity_x" key even though the associated
       array stores the x-component of the momentum density flux).

Note that the ``EnzoEquationOfState`` and ``EnzoIntegrableUpdate``
classes additionally require a ``EnzoEFltArrayMap`` object that hold the
passively advected scalars in conserved form.

In general, the use of ``EnzoEFltArrayMap`` objects with common sets of keys
helps simplify the implementation of various methods (e.g. the
cell-centered array associated with "density" is used to compute the
reconstruct values that are stored in the fields of the "density"
group in the reconstructed grouping).


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
dual-energy formalism see :ref:`using-vlct-de` ).

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

*In the future, the interface will need to be revisited once Grackle
is fully supported and it will be possible for gamma to vary
spatially.*

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

   apply_floor_to_energy_and_sync(EnzoEFltArrayMap &integrable_map,
                                  int stale_depth);

This method applies the applies the pressure floor to the total_energy
array specified in ``integrable_map``. If using the dual-energy formalism
the floor is also applied to the internal energy (also specified in 
``integrable_map``) and synchronizes the internal energy with the total
energy. If the equation of state is barotropic, this should do nothing.

.. code-block:: c++

   void pressure_from_integrable(EnzoEFltArrayMap &integrable_map,
                                 const EFlt3DArray &pressure,
                                 EnzoEFltArrayMap &conserved_passive_map,
                                 int stale_depth);

This method computes the pressure from the integrable quantities
(stored in ``integrable_map``) and stores the result in ``pressure``.
``conserved_passive_map`` should include the passive scalars in
conserved form.  The last argument currently doesn't do anything
and will only be important if Grackle is in use (with
``primordial_chemistry>1``). 

*In principle this should wrap* ``EnzoComputePressure``, *but
currently that is not the case. Some minor refactoring is needed to
allow EnzoComputePressure to compute Pressure based on arrays
specified in a* ``EnzoEFltArrayMap`` *object and we are holding off on
this until we implement full support for Grackle. Currently, when the
dual-energy_formalism is in use, pressure is simply computed from
internal energy.*

.. code-block:: c++

   void pressure_from_reconstructable(EnzoEFltArrayMap &reconstructable,
                                      EFlt3DArray &pressure,
                                      int stale_depth);

This method computes the pressure from the reconstructable quantities
(stored in ``reconstructable``) and stores the result in ``pressure``.

Note: for a non-barotropic equation of state, pressure is considered a
reconstructable quantity. In that case, if the pressure array in
``reconstructable`` is an alias of ``pressure``, nothing happens.
However, if the arrays aren't aliases, then the values are simply
copied between arrays.

.. code-block:: c++

   void reconstructable_from_integrable
     (EnzoEFltArrayMap &integrable, EnzoEFltArrayMap &reconstructable,
      EnzoEFltArrayMap &conserved_passive_map, int stale_depth,
      const std::vector<std::string> &passive_list);

This method is responsible for computing the reconstructable
quantities (to be held in ``reconstructable``) from the
integrable quantities (stored in ``integrable``). Note that
because of the high degree of overlap between the quantities in each
category, the overlapping quantities are assumed to be represented by
aliased arrays. Entries in ``passive_list`` specifies the passively
advected scalars that should be present in both maps.
The conserved form of the passively advected scalars
must be provided (stored in ``conserved_passive_map``) in case the
equation of state is barotropic and Grackle is in use. 

For a barotropic equation of state, this nominally does nothing while
for a non-barotropic equation of state, this just computes pressure by
calling ``EnzoEquationOfState::pressure_from_integrable``.

.. code-block:: c++

   void integrable_from_reconstructable
     (EnzoEFltArrayMap &reconstructable, EnzoEFltArrayMap &integrable,
      int stale_depth, const std::vector<std::string> &passive_list);

This method computes the integrable quantities (to be held in
``integrable``) from the reconstructable quantities (stored in
``reconstructable``). Again, because of the high degree of
overlap between the quantities in each category, the overlapping
quantities are assumed to be represented by aliased quantities.
Entries in ``passive_list`` specifies the passively
advected scalars that should be present in both maps.

For a barotropic equation of state, this nominally does nothing, while
for a non-barotropic equation of state, this nominally just computes
specific total energy. If the dual-energy formalism is in use this
also computes the internal energy.

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
    (const std::vector<std::string> active_reconstructed_quantities,
     std::string name, enzo_float theta_limiter);

The factory method requires that we register the names of the
reconstructable quantities via ``active_reconstructable_quantities``
and specify the name of the reconstruction algorithm, ``name``. Note
that the names of the reconstructable quantites should match
quantities specified in ``FIELD_TABLE`` ; for more details about
``FIELD_TABLE``, see :ref:`Centered-Field-Registry`

Public Interface
----------------
The main interface function provided by this class is:

.. code-block:: c++

    void reconstruct_interface
      (EnzoEFltArrayMap &prim_map, EnzoEFltArrayMap &priml_map,
       EnzoEFltArrayMap &primr_map, int dim, EnzoEquationOfState *eos,
       int stale_depth, const std::vector<std::string>& passive_list);

This function takes the cell-centered reconstructable primtive
quantities (specified by the contents of ``prim_map``) and computes
the left and right reconstructed states (the results are stored in
``priml_map`` and ``primr_map``) along the dimension specifed by
``dim``. If dim has a value of ``0``/ ``1``/ ``2`` then the values are
reconstructed along the x-/y-/z-axis. ``stale_depth`` indicates the
current stale_depth for the supplied cell-centered quantities (prior
to reconstruction). Note that the arrays in ``priml_map`` and
``primr_map`` should have arrays that are large enough to store
cell-centered quantitites so that they can be reused to hold the
face-centered fields along each dimension. ``passive_list`` is used to
specify the names (keys) of the passively advected quantities that are
to be reconstructed.

The ``int EnzoReconstructor::immediate_staling_rate()`` method is
provided to determine the amount by which the stale depth increases
immediately after reconstruction, for a given algorithm. The
``int EnzoReconstructor::delayed_staling_rate()`` method returns how much
the stale depth increases after adding flux divergence, computed from
the reconstructed values, to the integrable quantities  (this is
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

The Riemann Solvers have been factored out to their own classes. All
implementation of (approximate) Riemann solver algorithms are derived
from the ``EnzoRiemann`` abstract base class.


Usage Notes
-----------

To get a pointer to a concrete implemenation of ``EnzoRiemann``, call the
static factory method:

.. code-block:: c++

   EnzoRiemann* EnzoRiemann::construct_riemann
     (std::vector<std::string> integrable_quantities,
      std::string solver);

The factory method requires that we both register the names of the
integrable quantities (excluding passively advected scalars), with
``integrable_quantities``, and specify the name of the solver
``solver``. Note that the names of the integrable quantites should
match quantities specified in ``FIELD_TABLE`` that are identified as
being actively advected. For more details about ``FIELD_TABLE``, see
:ref:`Centered-Field-Registry`

The main interface function of ``EnzoRiemann`` is:

.. code-block:: c++

   void solve
     (EnzoEFltArrayMap &prim_map_l, EnzoEFltArrayMap &prim_map_r,
      const EFlt3DArray &pressure_array_l, const EFlt3DArray &pressure_array_r,
      EnzoEFltArrayMap &flux_map, int dim, EnzoEquationOfState *eos,
      int stale_depth, const std::vector<std::string> &passive_lists,
      EFlt3DArray *interface_velocity)

In this function, the ``prim_map_l`` and ``prim_map_r`` arguments are
references to the ``EnzoEFltArrayMap`` objects holding the arrays of
reconstructed left/right integrable quantities and passively advected
scalars. The ``pressure_array_l``/ ``pressure_array_r`` arguments
specify arrays holding the left/right reconstructed pressure. The
``flux_map`` argument holds the face-centered arrays where the
computed fluxes for each integrable quantity and passively advected
scalar will be stored. ``dim`` indicates the dimension along which the
flux should be computed (0,1,2 corresponds to x,y,z).
``interface_velocity`` is an optional argument used to specify a
pointer to an array that can be used to store interface velocity
values computed by the Riemann Solver (this is primarily used for
computing internal energy source terms when the dual energy formalism
is in use).


Implementation Notes: ``EnzoRiemannImpl``
-----------------------------------------

Historically, in many hydro codes (including Enzo) there is a lot of code
duplication between implementations of different types of Riemann Solvers
(e.g. converting left/right primitives to left/right conserved quantities
and computing left/right fluxes). To try to reduce some of this
duplication without sacrificing speed, we have defined the
``EnzoRiemannImpl<ImplFunctor>`` class template (which is a subclass of
``EnzoRiemann``).

The class template factors out common code shared by many approximate
Riemann Solvers (e.g. HLLE, HLLC, HLLD and possibly LLF & Roe solvers).
The template argument, ``ImplFunctor``, is a functor that implements
solver-specific calculations and is called at every cell-interface.
Additionally, the functor also specifies a specialization of the
template class ``EnzoRiemannLUT<InputLUT>`` that primarily

  * Specifies the exact set of actively advected integrable quantities
    that a given solver expects
  * Serves as a compile-time lookup table. It statically maps the names
    of the all of the components of the relevant actively advected
    quantities to unique array indices.

See :ref:`EnzoRiemannLUT-section`
for a more detailed description of ``EnzoRiemannLUT<InputLUT>`` and
examples of how it is used.

*Note: a more traditional inheritance-based approach that uses a
virtual method to implement solver-specific code. Calling a virtual
method in the main loop introduces overhead and prevents inlining.*

``EnzoRiemannImpl`` Control flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A brief overview of the ``EnzoRiemannImpl<ImplFunctor>::solve``
control flow is provided below. Basically the function loops over all
cell interfaces, along a given dimension, where the flux should be
computed. At each location, the following sequence of operations are
performed:

  1. Retrieve the left and right primitives at the given location from
     the input arrays and stores them in stack-allocated arrays of
     ``enzo_float`` elements called ``wl`` and ``wr``. As mentioned
     above, the values are organized according to the specialization
     of ``EnzoRiemannLUT<InputLUT>`` provided by the ``ImplFunctor``
     (hereafter, ``ImplFunctor::LUT``)
  2. The left and right pressure values are retrieved from the
     temporary fields holding the values that were precomputed from
     the reconstructed quantities (presumably using a concrete
     subclass of ``EnzoEquationOfState``). The values are stored in
     local variables ``pressure_l`` and ``pressure_r``.
  3. The conserved forms of the left and right reconstructed
     primitives and stored in the arrays called ``Ul`` and
     ``Ur``. Primitives that are always in conserved form (e.g.
     density or magnetic field). The elements of ``Ul`` / ``Ur``
     are also ordered by ``ImplFunctor::LUT`` (e.g. the index for a
     given component of the velocity in ``wl`` / ``wr`` matches the
     index for the same component of the momentum in ``Ul`` / ``Ur``).
  4. The standard left and right hydro/MHD fluxes are computed using
     the above quantities and stored in ``Fl`` and ``Fr`` (organized by
     ``ImplFunctor::LUT``)
  5. These quantities are all passed to the static public
     ``operator()`` method provided by ``ImplFunctor`` that returns the
     array of interface fluxes in the array, ``fluxes``. (It also
     computes the interface velocity)
  6. The interface fluxes and interface velocity are then copied into the
     output fields.

A separate method is provided to compute the fluxes for the passively
advected quantities. This method will also be compute the fluxes of any
specified quantities that are nominally actively advected, but can fall
back to using passive advection when the solver doesn't explictly support
it (the main example is ``"internal_energy"``)
     
*Note: Currently EnzoRiemannImpl has only been tested and known to
work for 3D problems. Additionally, no solvers (or more specifically,
wavespeed calculations) are currently implemented that explicitly
support barotropic equations of state (however, all of the machinery
is in place to support them).*

*Note: It might make sense to move calculation of conserved quantities
and fluxes into* ``ImplFunctor`` *. For some solvers, it may not be
necessary to compute all of this information. The template functions
that perform these operations have already been factored out into the*
``enzo_riemann_impl`` *namespace - so the transition would be easy to
accomplish.*

ImplFunctor template argument
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This subsection provides a brief description of the ``ImplFunctor``
template argument used to specialize ``EnzoRiemannImpl<ImplFunctor>``.
The class is expected to:

    * be default constructible

    * publically define the ``LUT`` type, which should be a specialization
      of the ``EnzoRiemannLUT<InputLUT>`` template class.
      ``ImplFunctor::LUT`` should indicate which actively advected
      quantities are expected by ``ImplFunctor`` and how they organized.
      For more details about how how ``EnzoRiemannLUT<InputLUT>`` is used,
      see :ref:`EnzoRiemannLUT-section`
           
    * provide the const-qualified function call method, ``operator()``.

The expected function signature of the ``operator()`` method is as follows:

.. code-block:: c++

   lutarray<ImplFunctor::LUT> ImplFunctor::operator()
     (const lutarray<ImplFunctor::LUT> flux_l,
      const lutarray<ImplFunctor::LUT> flux_r,
      const lutarray<ImplFunctor::LUT> prim_l,
      const lutarray<ImplFunctor::LUT> prim_r,
      const lutarray<ImplFunctor::LUT> cons_l,
      const lutarray<ImplFunctor::LUT> cons_r,
      enzo_float pressure_l, enzo_float pressure_r, bool barotropic_eos,
      enzo_float gamma, enzo_float isothermal_cs, enzo_float &vi_bar) const;

This function is called at every cell-interface and returns an array
holding the Riemann Flux at a given cell-interface. Note that
``lutarray<ImplFunctor::LUT>`` is actually an alias for
``std::array<enzo_float, ImplFunctor::LUT::NEQ>``. Each of these
arrays hold values associated with the components of each relevant
actively advected quantity and are organized according to
``ImplFunctor::LUT`` (again, see :ref:`EnzoRiemannLUT-section` for
more details about the ``LUT`` type).

``flux_l``/ ``flux_r``, ``prim_l``/ ``prim_r``, and ``cons_l``/
``cons_r`` store the left/right interface fluxes values, primitive
quantities, and conserved quantities (they are passed ``Fl``/ ``Fr``,
``wl``/ ``wr``, and ``Ul``/ ``Ur``, respectively).

The left and right reconstructed pressure values are passed as
``pressure_l`` and ``pressure_r``. ``barotropic_eos`` indicates
whether the fluid's equation of state is barotropic. If ``true``,
then ``isothermal_cs`` is expected to be non-zero and if ``false``,
then ``gamma`` is expected to be positive.

*Note: in the future, it would be worth experimenting with annotating the *
``operator()`` *method of ``ImplFunctor`` classes with the compiler
directive * ``__attribute__((always_inline))`` * to force inlining (this
works on g++, icc and clang).*

    .. _EnzoRiemannLUT-section:

EnzoAdvectionFieldLUT
~~~~~~~~~~~~~~~~~~~~~

As described above in the :ref:`GeneralDesignOverview-section` of the
General Design section, we sought to avoid the common approach of
hydro codes that map actively advected quantities indices with macros
or globally defined unscoped enums. The ``EnzoRiemannLUT<InputLUT>``
template class basically serves as a compromise between this traditional
approach approach and using a hash table (which introduce unacceptable
overhead) for organizing quantities in the main loop of
``EnzoRiemannImpl<ImplFunctor>``. Alternatively it can be thought of as a
scoped version of the traditional approach.

This is a template class that provides the following features at compile
time:

    * a lookup table (LUT) that maps the names of components of a subset
      of the actively advected quantities defined in ``FIELD_TABLE`` to
      unique, contiguous indices.

    * the number of quantity components included in the table

    * a way to iterate over just the conserved quantities or specific
      quantities values that are stored in an array using these mapping

    * a way to query which of the actively advected quantities in
      FIELD_TABLE are not included in the LUT

These feature are provided via the definition of publicly accessible
integer constants in every specialization of the template class. All
specializations have:

    * a constant called ``NEQ`` equal to the number of quantity components
      included in the lookup table

    * a constant called ``specific_start`` equal to the number of components
      of conserved quantities included in the lookup table

    * ``qkey`` constants, which include constants named for the components
      of ALL actively advected quantities in FIELD_TABLE. A constant
      associated with a SCALAR quantity, ``{qname}``, is simply called
      ``{qname}`` while constants associated with a vector quantity
      ``{qname}`` are called ``{qname}_i``, ``{qname}_j``, and ``{qname}_k``.

The `qkey` constants serve as both the keys of the lookup table and a
way to check whether a component of an actively advected quantity is
included in the table. Their values are satisfy the following conditions:

    * All constants named for values corresponding to quantities included
      in the table have values of ``-1``

    * All constants named for conserved quantities have unique integer
      values in the internal ``[0,specific_start)``

    * All constants named for specific quantities have unique integer
      values in the interval ``[specific_start, NEQ)``

The lookup table is always expected to include density and the 3 velocity
components. Although it may not be strictly enforced (yet), the lookup
table is also expected to include either all 3 components of a vector
quantity or None of them.

This template class also provides a handful of helpful static methods to
programmatically probe the table's contents at runtime and validate that
the above requirements are specified.

For the sake of providing some concrete examples about how the code works,
let's assume that we have a class ``MyInputLUT`` that is defined as:

.. code-block:: c++

   struct MyIntLUT {
     enum vals { density=0, velocity_i, velocity_j, velocity_k,
                 total_energy, NEQ, specific_start = 1};
   };

The template specialization ``EnzoRiemannLUT<MyIntLUT>`` assumes that
all undefined `qkey` constants omitted from ``MyIntLUT`` are not included
in the lookup table and will define them within the template specialization
to have values of ``-1``.

To access the index associated with density or the jth component of
velocity, one would evaluate:

.. code-block:: c++

   int density_index = EnzoRiemannLUT<MyInLUT>::density; //=0
   int vj_index = EnzoRiemannLUT<MyInLUT>::velocity_j;   //=2


It makes more sense to talk about the use of this template class when we
have a companion array. For convenience, the alias template
``lutarray<LUT>`` type is defined. The type,
``lutarray<EnzoRiemannLUT<InputLUT>>`` is an alias of the type
``std::array<enzo_float, EnzoRiemannLUT<InputLUT>::NEQ>;``.

As an example, imagine that the total kinetic energy density needs to be
computed at a single location from an values stored in an array, ``prim``,
of type ``lutarray<EnzoRiemannLUT<MyInLUT>>``:

.. code-block:: c++

   using LUT = EnzoRiemannLUT<MyInLUT>;
   enzo_float v2 = (prim[LUT::velocity_i] * prim[LUT::velocity_i] +
                    prim[LUT::velocity_j] * prim[LUT::velocity_j] +
   prim[LUT::velocity_k] * prim[LUT::velocity_k]);
   enzo_float kinetic = 0.5 * prim[LUT::density] * v2;


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


Adding new quantites
--------------------

To add support for new actively advected integrable cell-centered
quantities (e.g. cosmic ray energy/flux), the table of cell-centered
quantities (``FIELD_TABLE``) must be updated. See
:ref:`Centered-Field-Registry`
for more details.

To add support for computing fluxes for such quantities, modifications
must be made to ``EnzoRiemannImpl``. Currently, an abstract base class
called for ``EnzoFluxFunctor`` is provided for this purpose. The idea
is define a subclass to be defined for each additional set of flux
calculations and then in then have the factory method,
``EnzoRiemann::construct_riemann``, pass an array of the relevant
functors to ``EnzoRiemannImpl``.

*However, because the functors are called as pointers will probably
incur overhead. In reality, the better solution might be to hardcode
in the additonal flux calculation functions in some kind of helper
method of* ``EnzoRiemannImpl``.

Adding new solvers
------------------

New Riemann Solvers can currently be added to the infrastructure by
either subclasseding ``EnzoRiemann`` or defining a new specialization
of ``EnzoRiemannImpl<ImplFunctor>``. In either case, the
``EnzoRiemann::construct_riemann`` factory method must be modified to
return the new solver and :ref:`using-vlct-riemann-solver`
should be updated.

The additional steps for implementing a new Riemann solver by speciallizing
``EnzoRiemannImpl<ImplFunctor>`` are as follows:

  1. Define a new ``ImplFunctor`` class (e.g. ``HLLDImpl``)

  2. Add the new particlular specialization of ``EnzoRiemannImpl`` to
     enzo.CI (e.g. add the line:
     ``PUPable EnzoRiemannImpl<HLLDImpl>;``)

  3. *(optional)* define an alias name for the specialization of
     ``EnzoRiemannImpl`` that uses the new ``ImplFunctor`` class
     (e.g. ``using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDImpl>;``).

==============================
Updating integrable quantities
==============================

The ``EnzoIntegrableUpdate`` class has been provided to encapsulate
the operation of updating integrable quantities after a (partial)
time-step. The operation was factored out of the ``EnzoMethodMHDVlct``
class since it appear in all Godunov solvers.

The constructor for ``EnzoIntegrableUpdate`` has the following
signature:

.. code-block:: c++

   EnzoIntegrableUpdate(std::vector<std::string> integrable_groups,
		        bool skip_B_update)

The function requires that we:

  * register the names of the integrable quantities (with
    ``integrable_groups``)
  * indicate whether the update to the magnetic field should
    be skipped.

The names of the integrable quantites should match the names specified
in ``FIELD_TABLE``; see :ref:`Centered-Field-Registry` for more
details. The update to the magnetic field should be skipped when
Constrained Transport is in use (since the magnetic field update is
handled separately). If the magnetic field is not specified as an
integrable quantity, then the value specified for ``skip_B_update`` is
unimportant

The following method is used to compute the change in (the conserved
form of) the integrable and passively advected quantites due to the
flux divergence along dimension ``dim`` over the (partial) imestep
``dt``. The arrays in ``dUcons_map`` are used to accumulate the total
changes in these quantities. ``passive_list`` lists the names (keys)
of the passively advected scalars.

.. code-block:: c++

   void accumulate_flux_component
     (int dim, double dt, enzo_float cell_width, EnzoEFltArrayMap &flux_map,
      EnzoEFltArrayMap &dUcons_map, int stale_depth,
      const std::vector<std::string> &passive_list) const;

The method used to clear the values of the arrays used for accumulation is
provided below. This sanitization should be performed before starting
to accumulate flux divergence or source terms. The ``passive_list``
argument is used in the same way as the previous function.

.. code-block:: c++

    void clear_dUcons_group(EnzoEFltArrayMap &dUcons_map, enzo_float value,
                            const std::vector<std::string> &passive_list) const;

The method used to actually add the accumulated change in the integrable
(specified in ``dUcons_map``) to the values of the
integrable quantities from the start of the timestep (specificed by
``initial_integrable_map``) has the following signature:

.. code-block:: c++

   void update_quantities
     (EnzoEFltArrayMap &initial_integrable_map, EnzoEFltArrayMap &dUcons_map,
      EnzoEFltArrayMap &out_integrable_map,
      EnzoEFltArrayMap &out_conserved_passive_scalar,
      EnzoEquationOfState *eos, int stale_depth,
      const std::vector<std::string> &passive_list) const;

The fields included in ``dUcons_map`` should include contributions
from both the flux divergence AND source terms. The results for the
actively advected quanties are stored in ``out_integrable_map`` and
the results for the passively advected scalars are stored in conserved
form in the arrays held by ``out_conserved_passive_scalar`` (note that
the initial values of the passive scalars specified in
``initial_integrable_map`` are in specific form).

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
     (EnzoEFltArrayMap &cur_prim_map, EnzoEFltArrayMap &xflux_map,
      EnzoEFltArrayMap &yflux_map, EnzoEFltArrayMap &zflux_map,
      EnzoEFltArrayMap &out_centered_bfield_map, enzo_float dt,
      int stale_depth) noexcept;

In ``EnzoBfieldMethodCT`` this will also update the face-centered
magnetic field values (it assumes that ``identify_upwind`` was called
once for each dimension and uses the stored data). When using this
alongside ``EnzoIntegrableUpdate``, care needs to be taken about the
order in which this method is called relative to
``EnzoIntegrableUpdate::update_quantities`` that accounts for the time
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
