****************************
Hydro/MHD C++ Infrastructure
****************************

*[This page is under development]*

In this section we discuss some of the C++ infrastructure provided in
the Enzo Layer that can be optionally used to implement other
hydro/MHD methods. The infrastructure was used to implement other the
VL + CT MHD solver.

*Note: Currently, barotropic equations of state and compatibility with
Grackle are not yet implemented within the infrastucture. However they
are mentioned throughout this guide and slots have been explicitly left
open for them to be implemented within the framework. Additionally note
that while passively advected scalars are mostly supported, there is
not yet support for renomralizing the specific values of multiple
passively advected scalars to have a sum of 1.*

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
algorithm. This is subdivided into the "immediate staling rate" and
"delayed staling rate." The former specifies the amount by whch the
stale-depth increases after immediately after reconstruction
(e.g. this is 0 for nearest-neighbor interpolation and 1 for
piecewise-linear interpolation). The latter indicates that amount by
which the stale depth increase after adding the flux divergence
(e.g. 1 for both nearest-neighbor and piecewise linear interpolation).


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

The hydrodynamic/MHD C++ framework can be summarized as a series of
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
expected integrable or reconstructable quantities are *registered*
at construction.

  * The names of all reconstructable quantites that get registered
    in the construction of ``EnzoRiemann`` must share a name
    with the registered quantities in ``FIELD_TABLE``.

  * All registered integrable quantity names in the construction of
    ``EnzoRiemann`` or ``EnzoIntegrableUpdate`` must be specified in
    ``FIELD_TABLE`` as quantities that are actively advected in some
    contexts.

The expected groups of passively advected scalars are also specified at
construction. Known groups names of passive scalars can be retrieved
from the ``EnzoCenteredFieldRegistry`` class. For more information
about ``EnzoCenteredFieldRegistry`` and ``FIELD_TABLE`` see
:ref:`Centered-Field-Registry`

The implementation of these operation classes largely aims to avoid
employing following the traditional approach in which most field
data is directly accessed from a large array using macros or
globally defined unscoped enums that maps quantity component names
to indices. This traditional approach makes the introuduction of
optional fields that are related to active advection somewhat
difficult (e.g. cosmic ray energy/fluxes, internal energy for
dual energy formalism, phi for dedner divergence cleaning). In
place of this traditional approach leans heavily on Cello's
provided infrastructure for field data and make heavy usage of
of the ``Grouping`` class.

Use of ``Grouping``
-------------------

Overview
~~~~~~~~
The basic unit that get's operated on by these operation classes
are instances of Cello's ``Grouping`` class. We essentially use
them as containers of quantities (they hold the names of fields
related to those quantities).

The ``Grouping`` class was originally defined to organize field names
or particles types into named categories (or groups). A given field
name can be placed into more than one group.  The API primarily
supports adding fields (& particle types) to new or existing groups,
querying whether a field name belongs to a group, determining the
number of field names within a group and iterating over the fields of
a group. Note that the API does not currently provide a way to get the
names of all registered groups.

For the purposes of the hydrodynamic/MHD framework, the ``Grouping``
objects are used in more selective ways. Frequently, the groups are
used to serve as aliases for the fields that represent quantities
registered in operation classes. Aliases for scalar quantites (like
"density") are expected to hold a single related field. Aliases for
vector quantity (like velocity) are expected to hold a field for
each spatial component of the quantity. Instances of ``Grouping``
can also include groups that contain field names representing
passively advected scalars (e.g. you might have a collection of
fields for passively advected scalars in a group called "colour").

This described usage of ``Grouping`` objects amounts to crude
associative arrays (aka maps or dictionaries) that effectively maps
grop names to field data. Although they technically map the group
names to field names, frequently the field names are immediately
used to load the field data.

*Note: It would probably be beneficial to replace this usage of*
``Grouping`` *with an actual associative array* (*e.g.* ``std::map`` )
*that directly maps names to data. Doing so would reduce the
complexity of the code (and the amount of required documentation),
would reduce coupling of the hydro machinery to the cello block and
field machinery (making its usage moreflexible), and may even carry
some performance benefits.*

Specific Usage
~~~~~~~~~~~~~~

The names of groups expected in an instance of ``Grouping``
are the names of the quantities (and groups of passive scalars)
registered during the creation of the operation classes. Specific
instances of ``Grouping`` always contain fields that serve some
related kind of quantity. Below, we provide a description of the main
types of ``Groupings`` required for the provided operation classes:

  1. Primary grouping of cell-centered quantities.

     * This has groups named for all integrable and reconstructable
       quantities used by the integrator. We store them together due
       to the high degree of overlap between each categories. All
       groups named after integrable quantities should hold permanent
       fields that hold the values at the start of the time-step and
       get updated at the end of the time-step.
     * This also contains groups of passively advected scalars. Note
       that the fields contained within this group should all be
       temporary and they should all represent the passive scalars in
       "specific" form at the start of the timestep. The general
       convention (not just within this infrastructure) is for
       passively advected scalars to be primarily represented in
       "conserved" form (mass densities) outside of hydrodynamic
       integrator methods and to be converted to "specific" form (mass
       fractions) within the integrator methods
 
  2. Grouping of temporary cell-centered quantities

     * This grouping is identical to the above grouping (it must have
       all of the same groups of fields), except that the contained
       fields are used to temporarily hold quantities after partial
       time-steps. Based on the number of partial timesteps used by a
       particular method there might be 0 or multiple of these
       groupings.
     * Note that this type of grouping is used instead of the field
       history feature to avoid conflicts related to various
       ``Method`` objects (whether or not they directly implement
       hydro/MHD solvers) having different assumptions about the
       stored field history and to reduce the memory footprint.

  3. Grouping of temporary cell-centered for tracking the total change
     in a quantity over a timestep.

     * This grouping holds groups named for all integrable quantities and
       groups of passively advected scalars. For each (partial) timestep,
       the fields in the grouping are used to accumulate the total change
       in the conserved form of each quantity. This includes the flux
       divergence and the contributions from source terms. At the end of
       the (partial) timestep, these are used to actually update the
       values of the integrable quantities

  4. Groupings of reconstructed left/right quantites

     * 2 instances of groupings of this kind are used to respectively hold
       the reconstructed left and right interface quantities. This should
       contain all of the group names posessed in the above 2 groupings.
     * These groupings of fields are frequently passed to instances of
       ``EnzoReconstructor`` store the reconstructed passively
       advected scalars and reconstrutable quantities. They are then
       usually passed to ``EnzoEquationOfState`` to compute and store the
       reconstructed integrable quantities and reconstructed pressure.
       Then, these are frequently passed to ``EnzoRiemann`` to compute
       fluxes from the integrable quantities and the passively advected
       scalars.
     * Although this inherently represents face-centered data, the
       contained fields should be formally defined as
       cell-centered. This is done to allow for reuse of these fields
       to hold reconstructed fields along each dimension. This means
       that there is always some unused allocated memory at the end of
       the array allocated for each contained field.  The
       ``EnzoFieldArrayFactory::reconstructed_field`` method is
       provided to load the fields held by these groupings as
       ``CelloArray`` instances with the appropriate face-centered
       dimensions.

  5. Grouping of Riemann Flux fields

     * An instance of this kind of grouping is required for each
       dimension and is used to hold the face-centered fluxes along
       that dimension. The contained fields are all nominally
       temporary and should all be defined as face-centered along that
       dimension and they should not have space for values on the
       exterior of mesh blocks.
     * This kind of grouping should contain the names of all
       registered integrable fields and the registerred names of
       passively advected scalar groupings. The same integrable
       quantity names should be used here that are also used in the
       primary group, regardless of whether a quantity is "specific"
       or "conserved" (e.g. this kind of grouping always has
       "velocity" even though the contained fluxes are technically
       momentum density fluxes).

Note that the ``EnzoEquationOfState`` and ``EnzoIntegrableUpdate``
classes additionally require a ``Grouping`` object that hold the
passively advected scalars in conserved form.

In general, the use of ``Grouping`` objects with common sets of names
helps simplify the implementation of various methods (e.g. the
cell-centered field associated with "density" is used to compute the
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

   apply_floor_to_energy_and_sync(Block *block,
                                  Grouping &integrable_group,
                                  int stale_depth);

This method applies the applies the pressure floor to the total_energy
field specified by ``integrable_group``. If using the dual-energy formalism
the floor is also applied to the internal energy (also specified by the
``integrable_group``) and synchronizes the internal energy with the total
energy. If the equation of state is barotropic, this should do nothing.

.. code-block:: c++

   void pressure_from_integrable(Block *block,
                                 Grouping &integrable_group,
                                 std::string pressure_name,
                                 Grouping &conserved_passive_group,
                                 int stale_depth);

This method computes the pressure from the integrable quantities
(stored in ``integrable_group``) and stores the result in the field
specified by ``pressure_name``. The ``conserved_passive_group`` should
include the fields of passive scalars in conserved form.  This
currently doesn't do anything and will only be important if
Grackle is in use. 

*In principle this should wrap* ``EnzoComputePressure``, *but
currently that is not the case. Some minor refactoring is needed to
allow EnzoComputePressure to compute Pressure based on fields
specified in a* ``Grouping`` *object and we are holding off on this
until we implement full support for Grackle. Currently, when the
dual-energy_formalism is in use, pressure is simply computed from
internal energy.*

.. code-block:: c++

   void pressure_from_reconstructable(Block *block,
                                      Grouping &reconstructable_group,
                                      std::string pressure_name,
                                      int stale_depth,
                                      int reconstructed_axis);

This method computes the pressure from the reconstructable quantities
(stored in ``reconstructable_group``) and stores the result in the
field held by ``pressure_name``. ``reconstructed_axis`` is used to
specify if the fields are reconstructed. A value of -1 means that the
fields are cell-centered. A value of 0, 1, or 2 means that the fields
are reconstructed and they only contain valid values on the interior
x, y, or z faces of the mesh block.

Note: for a non-barotropic equation of state, pressure is considered a
reconstructable quantity. In that case, if the pressure field in
``reconstructable_group`` matches ``pressure_name``, nothing
happens. However, if the field names do not match, then the values are
simply copied between fields.

.. code-block:: c++

   void reconstructable_from_integrable (Block *block,
                                         Grouping &integrable_group,
                                         Grouping &reconstructable_group,
                                         Grouping &conserved_passive_group,
                                         int stale_depth);

This method is responsible for computing the reconstructable
quantities (to be held in ``reconstructable_group``) from the
integrable quantities (stored in ``integrable_group``). Note that
because of the high degree of overlap between the quantities in each
category, the overlapping quantities are assumed to be represented by
the same fields (this should get explicitly checked). The conserved
form of the passively advected scalars must be provided (stored in
``conserved_passive_group``) in case the equation of state is
barotropic and Grackle is in use.

For a barotropic equation of state, this nominally does nothing while
for a non-barotropic equation of state, this just computes pressure by
calling ``EnzoEquationOfState::pressure_from_integrable``.

.. code-block:: c++

   void integrable_from_reconstructable(Block *block,
                                        Grouping &reconstructable_group,
                                        Grouping &integrable_group,
                                        int stale_depth,
                                        int reconstructed_axis);

This method computes the integrable quantities (to be held in
``integrable_group``) from the reconstructable quantities (stored in
``reconstructable_group``). Again, because of the high degree of
overlap between the quantities in each category, the overlapping
quantities are assumed to be represented by the same fields (this
should again get explicitly checked). ``reconstructed_axis`` is used to
specify if the fields are reconstructed. A value of -1 means that the
fields are cell-centered. A value of 0, 1, or 2 means that the fields
are reconstructed and they only contain valid values on the interior
x, y, or z faces of the mesh block.

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
    (std::vector<std::string> reconstructable_groups,
     std::vector<std::string> passive_groups, std::string name);

The factory method requires that we register the names of the
reconstructable quantities (with ``reconstructable_groups``), register
the names of the groups containing the passively advected quantities
(with ``passive_groups``) and specify the name of the reconstruction
algorithm, ``name``. Note that the names of the reconstructable
quantites should match quantities specified in ``FIELD_TABLE`` ; For
more details about ``FIELD_TABLE``, see :ref:`Centered-Field-Registry`

Public Interface
----------------
The main interface function provided by this class is:

.. code-block:: c++

    void reconstruct_interface (Block *block, Grouping &prim_group,
                                Grouping &priml_group, Grouping &primr_group,
				int dim, EnzoEquationOfState *eos,
				int stale_depth)

This function takes the cell-centered reconstructable primtive
quantities (specified by the fields in ``prim_group``) and computes
the left and right reconstructed states (the results are stored in
``priml_group`` and ``primr_group``) along the dimension specifed by
``dim``. If dim has a value of ``0``/ ``1``/ ``2`` then the values are
reconstructed along the x-/y-/z-axis. ``stale_depth`` indicates the
current stale_depth for the supplied cell-centered quantities (prior
to reconstruction). Note that the fields in ``priml_group`` and
``primr_group`` should be formally defined as cell-centered so that
they can be reused to hold the face-centered fields along each
dimension.

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
implementation of reconstruction algorithms are derived from the
``EnzoRiemann`` abstract base class.


Usage Notes
-----------

To get a pointer to a concrete implemenation of ``EnzoRiemann``, call the
static factory method:

.. code-block:: c++

   EnzoRiemann* construct_riemann(std::vector<std::string> integrable_groups,
                                  std::vector<std::string> passive_groups,
                                  std::string solver);

The factory method requires that we register the names of the integrable
quantities (with ``integrable_groups``), register the names of the groups
containing the passively advected quantities (with ``passive_groups``)
and requires that we specify the name of the solver ``solver``. Note that
the names of the integrable quantites should match quantities specified in
``FIELD_TABLE`` that are identified as being actively advected. For more
details about ``FIELD_TABLE``, see
:ref:`Centered-Field-Registry`

The main interface function of ``EnzoRiemann`` is:

.. code-block:: c++

   void solve (Block *block,
               Grouping &priml_group, Grouping &primr_group, 
	       std::string pressure_name_l,
               std::string pressure_name_r,
               Grouping &flux_group, int dim,
               EnzoEquationOfState *eos,
               int stale_depth, std::string interface_velocity_name = "");

In this function, the ``priml_group`` and ``primr_group`` arguments
are references to the ``Grouping`` objects holding the reconstructed
left/right integrable quantity fields and groups of passively advected
scalar fields. The ``pressure_name_l``/ ``pressure_name_r`` arguments
hold the names of the left/right reconstructed pressure. The
``flux_group`` argument holds the face-centered fields where the
computed fluxes for each integrable quantity and passively advected
scalar will be stored. ``dim`` indicates the dimension along which the
flux should be computed (0,1,2 corresponds to x,y,z).
``interface_velocity_name`` is an optional argument used to specify
the name of a field that can be used to store interface velocity
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
    of the all of the components of the relevant actively advected fields
    to unique array indices.

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
     the fields and stores them in stack-allocated arrays of
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
		        bool skip_B_update,
		        std::vector<std::string> passive_groups)

The function requires that we:

  * register the names of the integrable quantities (with
    ``integrable_groups``)
  * indicate whether the update to the magnetic field should
    be skipped.
  * register the names of the groups containing the passively
    advected quantities (with ``passive_groups``).

The names of the integrable quantites should match the names specified
in ``FIELD_TABLE``; see :ref:`Centered-Field-Registry`
for more details. The update to the magnetic field should
be skipped when Constrained Transport is in use (since the magnetic
field update is handled separately). If the magnetic field is not
specified as an integrable quantity, then the value specified for
``skip_B_update`` is unimportant

The following method is used to compute the change in (the conserved
form of) the integrable and passively advected quantites due to the
flux divergence along dimension ``dim`` over the (partial) imestep
``dt``. These values are added to the the fields used to accumulate
the total changes in these quantities (specified by ``dUCons_group``).

.. code-block:: c++

   void accumulate_flux_component(Block *block, int dim, double dt,
                                  Grouping &flux_group, Grouping &dUcons_group,
                                  int stale_depth) const;

The method used to clear the values of the fields for accumulation is
provided below. This sanitization should be performed before starting
to accumulate flux divergence or source terms.

.. code-block:: c++

    void clear_dUcons_group(Block *block, Grouping &dUcons_group,
                            enzo_float value) const;

The method used to actually add the accumulated change in the integrable
(specified in ``dUcons_group``) to the values of the
integrable quantities from the start of the timestep (specificed by
``initial_integrable_group``) has the following signature:

.. code-block:: c++

   void update_quantities(Block *block, Grouping &initial_integrable_group,
                          Grouping &dUcons_group,
                          Grouping &out_integrable_group,
                          Grouping &out_conserved_passive_scalar,
                          EnzoEquationOfState *eos, int stale_depth) const;

The fields included in ``dUcons_group`` should include contributions from
both the flux divergence AND source terms. The results for the actively
advected quanties are stored in ``out_integrable_group`` and the results for
the passively advected scalars are stored in conserved form in the fields
held by ``out_conserved_passive_scalar`` (note that the initial values of
the passive scalars specified in ``initial_integrable_group`` are in
specific form).
