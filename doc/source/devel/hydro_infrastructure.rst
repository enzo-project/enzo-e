****************************
Hydro/MHD C++ Infrastructure
****************************

*[This page is under development]*

In this section we discuss some of the C++ infrastructure provided in the Enzo
Layer that was used to implement the VL + CT MHD solver that can be optionally
re-used to implement other hydro/MHD methods

*Note: Currently barotropic equations of state, compatibility with Grackle, and the dual energy formalism are not yet implemented within the infrastucture, however mentions are made about them in this guide and spots have been explicitly left for them to be implemented within the framework. Additionally note that while passively advected scalars are mostly supported, there is not yet support for renomralizing the specific values of multiple passively advected scalars to have a sum of 1.*

*Note: Currently brief summaries of the interfaces of each of the objects are provided below. More detailed descriptions are provided in the header files using doxygen documentation. If we can figure out how to generate reference documenation from doxygen (using breathe), then the summaries below should be deleted.*

===============
Shorthand Terms
===============

Here we briefly define a few terms that we defined and used throughout the
documentation and codebase

Integrable/Reconstructable Quantities
-------------------------------------

Throughout this guide and the relevant sections of the codebase, we
refer to quantities as reconstructable (i.e. they are used for
reconstructing left/right interface states) and integrable
quantities (i.e. the primary quantities that are evolved by the
integrator). There is frequently a high degree of overlap between these
categories and the precise categorization depends on the equation of
state. As of now all integrable quantities are "conserved" or
"specific" (a quantity like velocity that when multiplied by density
becomes conserved).

As an example, the categorization of the quantities for an ideal,
adiabatic gas are:

  * density - reconstructable and integrable

  * velocity - reconstructable and integrable

  * pressure - only reconstructable

  * (specific) total energy - only integrable

  * magnetic field - reconstructable and integrable

*Note: For historical reasons, both the reconstructable and integrable quantities are frequently referred to as primitives throughout the codebase as primitives. There mainly refers to the fact that the collection of quantities are not all of the quantities are conserved (at least some subset of them are "specific")*

stale depth
-----------

To help simplify the implementation of several operations, we
introduce the concept of "stale depth". At the start of a time-step,
there are never any "stale" quantities (stale depeth is
zero). However, every time flux divergence get's added to any
quantities, the outermost layer of up-to-date quantities (in the ghost
zone) becomes invalid; this happens because on the exterior faces of
the layer are not accurately known. We refer to these invalid values
as "stale" and say that the stale depth increased by 1. The stale
depth can also be incremented by other operations (e.g. piecewise
linear reconstruction). At the end of a time-step, the stale depth
should be equal to or less than the ghost depth so that the all of the
"stale" values can be refreshed (resetting the "stale depth" to zero).

We formally define "stale depth" as the number of the layers of the
outermost field entries that include "stale" values (for cell-centered
fields this is the number of cells from the edge that include stale
values). This means that for a given stale depth:

  * A face-centered field that has values on the exterior of a mesh
    block will always have one more unstaled value along the axis of
    face-centering than a cell-centered field.

  * A face-centered field that doesn't have values on the exterior of
    a mesh block will always have one less unstaled value along the
    axis of face-centering than a cell-centered field.

The introduction of this formalism has 2 key benefits:

  1. Simplifies calculation of correct number of ghost zones for unigrid
     simulations

  2. When used alongside ``CelloArray``, it drastically simplifies the
     calculation of which indices should be iterated over. The
     ``EnzoFieldArrayFactory`` can take the stale depth as an argument
     in its constructor and then any array it builds has the stale
     values clipped off.  This allows for-loop bounds to be written as
     though the only reconstruction algorithm is nearest-neighbor
     interpolation and as though there never any preceeding partial
     timesteps.


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
codebase and track some basic meta-data about the fields. The registry
was primarily created to help manage instances of the
``EnzoAdvectionFieldLUT`` struct (see below). Although it's
functionallity is presently limited, the ``EnzoCenteredFieldRegistry``
has the potential to be a general purpose tool that can be used for
other purposes in the Enzo layer of the codebase.

The idea is to maintain a list of all quantities, represented by
cell-centered fields, that are used by Enzo in the ``FIELD_TABLE``
macro. For each quantity, the table currently tracks its name,
whether its fundamentally a scalar or vector quantity, if the
quantity can be classified as "conserved", "specific", or "other",
and whether or not there is any circumstance where it is an
"actively" advected quantity. An entry about a scalar quantity
registers a field of the same name while an entry about a vector
quantity registers 3 fields with the names ``"{qname}_x"``, ``"{qname}_y"``,
``"{qname}_z"`` where ``{qname}`` is the name of the quantity (e.g.
the row for the "velocity" quantity registers the ``"velocity_x"``,
``"velocity_y"``, and ``"velocity_z"`` fields).

The registry was primarily motivated to perform operations related
to the ``EnzoAdvectionFieldLUT`` struct (see below). The registry
also tracks a list of known groups of passively scalars that can be
specified in the input file to indicate that the particular field
is passively advected (as of now, the only such group is
``"colour"``). Although it's functionallity is currently limited,
the ``EnzoCenteredFieldRegistry`` has the capacity to be a general
purpose tool useful other parts of the Enzo layer of the codebase.

.. _EnzoAdvectionFieldLUT-section:

EnzoAdvectionFieldLUT
---------------------

Motivation
~~~~~~~~~~

To implement a Riemann solver, most of the advected quantities need
to be known all at a given cell-interface to help compute
the wave speeds.

Unlike Enzo (and many other codes) we wanted to
avoid statically declaring which indices of an array correspond to
specific fields. Doing so complicates situations where different
combinations of fields can be optionally used (e.g. like optional
magnetic fields, dropping total energy for barotropics equations of
state, optionally advecting internal energy for dual energy
formalism or tracking cosmic ray energy density and cosmic ray
energy fluxes). 

It's also useful to be able to iterate over all of values at a given
interface after computing the wavespeeds to reduce the code necessary
to compute the Riemann Fluxes (and make it easier to add new fields).
It turns out that standard hash tables are too slow for these
operations.

Description
~~~~~~~~~~~

The above requirements motivate the creation of the
``EnzoAdvectionFieldLUT`` class is a barebones class designed to
be as close as possible to a C-style struct with members named for
each quantity registered in ``FIELD_TABLE`` that is actively advected
in Enzo in some context. Scalar quantities have members named
directly after them and vector quantities have 3 members named after
them: ``{qname}_i``, ``{qname}_j``, ``{qname}_k``, where ``qname``
is the name of the quantity (e.g. defining the velocity quantity in
``FIELD_TABLE`` causes the creation of the ``velocity_i``,
``velocity_j``, and ``velocity_k`` members). 

``EnzoCenteredFieldRegistry::prepare_advection_lut`` is used to
prepare an instance of ``EnzoAdvectionFieldLUT`` for a set of
integrable quantities. The function determines the length of an
array large enough to hold all of the relevant fields and
initializes each relevant struct member to have an integer value
corresponding to a unique index in the aforementioned array.
Any memebers of the class that don't correspond to a specified
integrable quantity has its value set to -1.

After being setup, instance of ``EnzoAdvectionFieldLUT`` can be
used as a lookup table. For example, imagine ``wl`` is an array of
reconstructed integrable primitives at a given cell interface
where the values are ordered according to a properly initialized
instance of ``EnzoAdvectionFieldLUT`` called ``lut``. In this case,
``wl[prim_lut.density]`` and ``wl[prim_lut.total_energy]`` give the
values of the density and (specific) total energy. 

As an added bonus, ``EnzoCenteredFieldRegistry::prepare_advection_lut``
orders the quantities based on whether they are "conserved",
"specific", or "other" (note: there shouldn't be any actively advected
quantities classified as "other"), and returns iteration limits of the
array that includes each class of value. This simplifies the conversion
of all integrable quantities to "conserved" form in a Riemann Solver
and when the integrable quantities are updated.

Additionally, ``EnzoCenteredFieldRegistry::prepare_advection_lut``,
allows for the specification of "flagged" quantities. When a quantity
is "flagged" the ordering of the indices assigned to the the members
of ``EnzoAdvectionFieldLUT`` and the iteration limits returned by the
function are modified such that the iteration limits don't include
the "flagged" quantities. This is useful when updating integrable
quantities in certain cases (e.g. if the dual energy formalism is in
use or if magnetic fields get updated by constrained transport).

``EnzoCenteredFieldRegistry::load_array_of_fields`` is
provided to assist with the loading of instances arrays of
``CelloArray`` encapsulating field data that are ordered according to
the ordering specified in an instance of ``EnzoAdvectionFieldLUT``.

*Note: Currently, the usage of this class is confined to the implementation of the* ``EnzoRiemannImpl`` *and* ``EnzoIntegrableUpdate`` *classes.*

==============
General Design
==============

Overview
--------

The hydrodynamic/MHD C++ framework can be summarized as a series of
classes that encapsulate various operations that are performed in the
hydrodynamic/MHD integrators. In most cases an abstract base class
exists to provide the interface for each operation. The main operation
classes include:

  * ``EnzoEquationOfState`` - encapsulate many of the operations
    related to the fluid's equation of state (e.g. computing pressure,
    converting reconstructable quantities to integrable quantities and
    vice-versa)

  * ``EnzoReconstructor`` - encapsulates interpolation algorithms to
    reconstruct left/right interface states of from cell-centered
    values

  * ``EnzoRiemann`` - encapsulates various Rimann Solver algorithms

  * ``EnzoIntegrableUpdate`` - encapsulates the operation of updating
    integrable quantities after a (partial) time-step.

Each of these operation classes are fairly modular (which means that
parts of the infrastructure can be selectively used); the only caveat
is that most classes use an instance of ``EnzoEquationOfState``.

Each of the operation classes are designed to be configured upon
initialization. The instances can then be used multiple times per
time-step (along multiple dimensions if the operation is directional)
and in other time-steps. The operation classes are also provided with
``PUP`` methods to allow for easy migration alongside the ``Method``
class that makes use them.

For each operation class (other than ``EnzoEquationOfState``), the
expected integrable or reconstructable quantities are specified at
construction.  The names of all reconstructable quantites that get
registered in the construction of ``EnzoRiemann`` must share a name
with the registered quantities in ``FIELD_TABLE``.  All registered
integrable quantity names in the construction of ``EnzoRiemann`` or
``EnzoIntegrableUpdate`` must be specified in ``FIELD_TABLE`` as a
quantity that is actively advected in some context. The expected
groups of passively advected scalars are also specified at
construction. Known groups names of passive scalars can be retrieved
from the ``EnzoCenteredFieldRegistry`` class. For more information
about ``EnzoCenteredFieldRegistry`` and ``FIELD_TABLE`` see
:ref:`Centered-Field-Registry`


Groupings
---------

The basic unit that get's operated on by these operation classes
are instances of Cello's ``Grouping`` class. We essentially use
them as containers of quantities (they hold the names of fields
related to those quantities).

The ``Grouping`` class was originally defined
to organize field names or particles types into named categories
(or groups). A given field name can be placed into more than one group.
The API primarily supports adding fields (& particle types) to existing
or new groups, querying whether the field name belongs to a group,
determining the number of field names within a group and iterating over
the fields in the group. Note that API does not currently provide a way
to get the names of all registered groups.

For the purposes of the hydrodynamic/MHD framework, the ``Grouping``
objects are used in more selective ways. Frequently, the groups are
used to serve as aliases for quantites used in integration. Aliases
for scalar quantites (like "density") are expected to hold a single
field related to the quantity while a vector quantity (like velocity)
is expected to hold a field for each component of the
quantity. Instances of ``Grouping`` also include groups that contain
field names representing passively advected scalars (e.g. you might
have a collection of fields in a group called "colour").

The required names of groups within an instance of ``Grouping`` are
given by the quantity names and names of groups of passively advected
scalars that are passed to the various operation classes at
construction. Specific instances of ``Grouping`` always fields that
with some related purpose. Below, we provide a description of the main
types of ``Groupings`` required for the provided operation classes:

  1. Primary group of cell-centered quantities.

     * This has groups named for all integrable and reconstructable quantities
       used by the integrator. We store them together due to the high degree of
       overlap between each categories. All groups named after integrable
       quantities should hold permanent fields that hold the values at the
       start of the time-step and get updated at the end of the time-step.

     * This also contains groups of passively advected scalars. Note that the
       fields contained within this group should all be temporary and they
       should all represent the passive scalars in "specific" form at the
       start of the timestep. The convention is for passively advected scalars
       to be primarily represented in "conserved" form (mass densities)
       outside of hydrodynamic integrator methods and to be converted to
       "specific" form (mass fractions) within integrator methods
 
  2. Group of temporary cell-centered quantities.

     * This grouping is identical to the above grouping (it must have
       all of the same groups of fields), except that the contained
       fields are used to temporarily hold quantities after partial
       time-steps. Based on the number of partial timesteps used by a
       method there might be 0 or multiple of these groupings.
     * Note that this type of grouping is used instead of the field history
       feature to avoid conflicts related to various ``Method`` objects
       (whether or not they directly implement hydro/MHD solvers) having
       different assumptions about the stored field history.

  3. Reconstructed left/right fields

     * 2 instances of groupings of this kind are used to respectively hold
       the reconstructed left and right interface quantities. This should
       contain all of the group names posessed in the above 2 groupings.
     * These sets of fields are frequently passed to ``EnzoReconstructor``
       to hold store the reconstructed passively advected scalars and
       reconstrutable quantities. They are then usually passed to
       ``EnzoEquationOfState`` to compute and store the reconstructed
       integrable quantities and reconstructed pressure. Then, these are
       frequently passed to ``EnzoRiemann`` to compute fluxes for each
       quantitiy.
     * Although this inherently represents face-centered data, the fields
       contained should be formally defined as cell-centered. This is done
       to allow for reuse of these fields to hold reconstructed fields along
       each dimension. This means that there is always some unused allocated
       memory at the end of the array allocated for each contained field.
       The ``EnzoFieldArrayFactory::reconstructed_field`` method is provided
       to load the fields held by these groupings as ``CelloArray`` instances
       with the appropriate face-centered dimensions.

  4. Riemann Flux fields

     * An instance of this kind of grouping is required for each
       dimension and is used to hold the face-centered fluxes along
       that dimension (the contained fields are all nominally
       temporary and should all be defined as face-centered along that
       dimension and they should not have space for values on the
       exterior of mesh blocks).
     * This kind of grouping should contain the names of all
       registered integrable fields and the registerred names of
       passively advected scalar groupings. The same integrable
       quantity names should be used here that are also used in the
       primary group, regardless of whether a quantity is specific
       or conserved (e.g. this kind of grouping always has "velocity"
       even though the contained fluxes are technically momentum
       density fluxes).

Note that the ``EnzoEquationOfState`` and ``EnzoIntegrableUpdate``
classes additionally require a ``Grouping`` object that hold the
passively advected scalars in conserved form.

In general, the use of ``Grouping`` objects with common sets of names
helps simplify the implementation of various methods (e.g. the
cell-centered field associated with "density" is used to compute the
reconstruct values that are stored in the fields of the "density"
group in the reconstructed grouping).

*Note: The use of the* ``Grouping``
*objects over a more traditional approach where an array of pointers (to the field data), where macros/enums are statically defined that map array indices to quantity names was originally motivated by the greater flexibility of adding new,optional quantities to the integrator (e.g. cosmic ray energy/fluxes).*



=================
Equation Of State
=================

All of the operations related to the equation of state are handled by
subclasses of the abstract base class, ``EnzoEquationOfState``. The
class has a number of responsibilities. Currently the only concrete
subclass of ``EnzoEquationOfState`` is the ``EnzoEOSIdeal`` class
which encapsulates the properties of an ideal, adiabatic gas.

The ``EnzoEquationOfState`` has the following interface:

.. code-block:: c++

   bool is_barotropic();

Returns whether the equation of state is barotropic or not.
*Currently, no barotropic equations of state have been implemented and
none of the wavespeed calculations for the Riemann solvers currently
support barotropic equations of state.*

.. code-block:: c++

   bool uses_dual_energy_formalism();

Returns whether the dual energy formalism is in use. *Currently, the dual energy formalism is not supported.*

.. code-block:: c++

   enzo_float get_gamma();

Returns the ratio of the specific heats. This is only required to yield a
reasonable value if the gas is not barotropic. *In the future, the interface will need to be revisited once Grackle is fully supported and it will be possible for gamma to vary spatially.*

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

   void apply_floor_to_total_energy(Block *block,
                                    Grouping &integrable_group,
				    int stale_depth);

This method makes sure that the field in the total_energy grouping of
the integrable group satisfies the pressure floor. This should do
nothing for a barotropic equation of state.

.. code-block:: c++

   pressure_from_integrable(Block *block,
                            Grouping &integrable_group,
                            std::string pressure_name,
                            Grouping &conserved_passive_group,
                            int stale_depth)

This method computes the pressure from the integrable quantities
(stored in ``integrable_group``) and stores the result in the field
specified by ``pressure_name``. The ``conserved_passive_group`` should
include the fields of passive scalars in conserved form.  This
currently doesn't do anything and will only be important if
Grackle is in use. *In principle this should wrap*
``EnzoComputePressure``,
*but currently that is not the case. (some small refactoring needs to be performed to allow EnzoComputePressure to compute Pressure based on fields specified in a* ``Grouping`` *object. We are holding off on this until we implement support for Grackle.* 

.. code-block:: c++

   void pressure_from_reconstructable(Block *block,
                                      Grouping &reconstructable_group,
                                      std::string pressure_name,
                                      int stale_depth,
                                      int reconstructed_axis)

This method computes the pressure from the reconstructable quantities
(stored in ``reconstructable_group`` and stores the result in the
field held by ``pressure_name``. ``reconstructed_axis`` is used to
specify if the fields are reconstructed. A value of -1 means that the
fields are cell-centered. A value of 0, 1, or 2 means that the fields
are reconstructed and they only contain valid values at x, y, or z
faces.

For a non-barotropic equation of state, pressure is considered a
reconstructable quantity. In that case, if the pressure field in
``reconstructable_group`` matches ``pressure_name``, nothing
happens. But if the field names do not match, then the values are
simply copied.

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
for a non-barotropic equation of state, this nominally just computes
the pressure.

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
are reconstructed and they only contain valid values at x, y, or z
faces.

For a barotropic equation of state, this nominally does nothing, while
for a non-barotropic equation of state, this nominally just computes
specific total energy.


How to extend
-------------
New equations of state can be added by subclassing and providing the subclass
with implementations for the pure virtual functions
``EnzoEquationOfState``

=============
Reconstructor
=============

The reconstruction algorithms have been factored out to their own classes. All
implementation of reconstruction algorithms are derived from the
``EnzoReconstructor`` abstract base class.

Pointers to concrete derived objects that encapsulate a particular object can
be retrieved using the ``EnzoReconstructor::construct_reconstructor`` static
factory method.

.. code-block:: c++

   EnzoReconstructor* construct_reconstructor
    (std::vector<std::string> reconstructable_groups,
     std::vector<std::string> passive_groups, std::string name);

The factory method requires that we register the names of the
reconstructable quantities (with ``reconstructable_groups``), register
the names of the groups containing the passively advected quantities
(with ``passive_groups``) and requires that we specify the name of the
reconstruction algorithm, ``name``. Note that the names of the
reconstructable quantites should match quantities specified in
``FIELD_TABLE`` ; For more details about ``FIELD_TABLE``, see
:ref:`Centered-Field-Registry`

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
provided to determine how much the stale depth increases immediately
after using a particular reconstruction algorithm. The
``int EnzoReconstructor::delayed_staling_rate()`` method returns how much
the stale depth increases after adding flux divergnece to the
integrable quantities computed using the reconstructed values (this is
normally 1). Finally ``int EnzoReconstructor::total_staling_rate()``
gives the sum of the results yielded by the prior 2 methods.

How to extend
-------------

To add a new reconstructor, subclass ``EnzoReconstructor`` and provide
definitions for the virtual methods. The factory method
``EnzoReconstructor::construct_reconstructor`` must also be modified
to return pointers to instances of the new class when the appropriate
name is passed as an argument. Additionally, implementations of the
``immediate_staling_rate()`` and ``total_staling_rate()`` virtual
methods must be provided

To take an existing reconstructor and make a new slope limiter available, a
different class should probably be declared. But, a system reminiscent of the
approximate RiemannSolvers could potentially be adopted to reduce redundant
code.

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

   EnzoRiemann* EnzoRiemann::construct_riemann
   (std::vector<std::string> integrable_groups,
    std::vector<std::string> passive_groups,
    std::string solver);

The factory method requires that we register the names of the integrable
quantities (with ``integrable_groups``), register the names of the groups
containing the passively advected quantities (with ``passive_groups``)
and requires that we specify the name of the solver ``solver``. Note that
the names of the integrable quantites should match quantities specified in
``FIELD_TABLE`` that are identified as being actively advected in some
context. For more details about ``FIELD_TABLE``, see
:ref:`Centered-Field-Registry`

To actually use the Riemann solver, the virtual solve method should be called.
The signature for this method is

.. code-block:: c++

   void EnzoRiemann::solve (Block *block,
                            Grouping &priml_group, Grouping &primr_group, 
			    std::string pressure_name_l,
			    std::string pressure_name_r,
			    Grouping &flux_group, int dim,
			    EnzoEquationOfState *eos,
			    int stale_depth);

In this function, the ``priml_group`` and ``primr_group`` arguments
are references to ``Grouping`` objects that have groups named for each
item in the ``integrable_group`` and ``passive_groups`` arguments
originally passed to the factory method. Within each group, the
``Grouping`` objects should contain the field names holding the left
and right reconstructed values that represent a quantity. For
``SCALAR`` quantites held in integrable_group, a group should only
have 1 field while for ``VECTOR`` quantites, a group should have 3
fields (labelled for the x, y, and z components).


Implementation Notes
--------------------

Historically, when Enzo (and many other codes) have implemented multiple
Riemann Solvers, there has been a large amount of code duplication
(e.g. converting left/right primitives to left/right conserved quantities
and computing left/right fluxes). To try to reduce some of the code
duplication without sacrificing speed, we have defined the
``EnzoRiemannImpl<ImplStruct>`` class template (which is a subclass of
``EnzoRiemann``).

Basically, the idea is that ``EnzoRiemannImpl<ImplStruct>`` class
template factors out duplicate code shared by many approximate Riemann
Solvers (e.g. HLLE, HLLC, HLLD and possibly LLF & Roe solvers). The
``ImplStruct`` is a simple struct/class that actually implements a
method that is responsible for the different code in each type of
solver and that gets called to compute the flux at every cell
interface. The more traditional object-oriented approach would have
been to make ``EnzoRiemannImpl`` an abstract class with a virtual
method reponsible for the solver-specific code. However, the act of
looking up the virtual method causes a performance hit and prevents
the code from being inlined within the main loop.

EnzoRiemannImpl Control flow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We note that at each location, arrays are constructed to hold
different integrable quantites and an instance of the
``EnzoAdvectionFieldLUT`` struct is used as a lookup table.
Basically, the struct has members named after all potential
advectable, integrable quantities and the members corresponding to the
registered integrable quantities assigned values that correspond to
indices in an array. See :ref:`EnzoAdvectionFieldLUT-section`
for a
more detailed description of the struct. This struct is used inplace
of the more tradional approach of defining global enums or macros to
map quantity names to array indices. Note that we when we tried to use
a built-in hash tables to perform this task, had a significant
performance cost.

Below, a brief overview of the ``EnzoRiemannImpl::solve`` control flow
is provided. Basically the function loops over all cell interfaces,
along a given dimension, where the flux should be computed. At each
location, the following sequence of operations are performed:

  1. Retrieve the left and right primitives at the given location from
     the fields and stores them in arrays of ``enzo_float`` elements
     called ``wl`` and ``wr``. The elements are ordered based on a
     preconfigured instance of ``EnzoAdvectionFieldLUT`` called
     ``lut_``.
  2. The left and right pressure values are retrieved from the
     temporary fields holding the values that were precomputed from
     the reconstructed quantities (presumably using a concrete
     subclass of ``EnzoEquationOfState``). The values are stored in
     ``pressure_l`` and ``pressure_r``.
  3. The conserved forms of the left and right reconstructed
     primitives and stored in the arrays called ``Ul`` and
     ``Ur``. Note that the primitives that are always in conserved
     form (e.g. density or magnetic field) are simply copied over. The
     elements of ``Ul`` / ``Ur`` maintain the same ordering as those
     of ``wl`` / ``wr`` (e.g. the index for a given component of the
     velocity in ``wl`` / ``wr`` is the index for the same component
     of the momentum in ``Ul`` / ``Ur``).
  4. The standard left and right hydro/MHD fluxes are computed using
     the above quantities and stored in ``Fl`` and ``Fr``; the
     elements are again ordered by ``lut_``.
  5. In principle, non-standard fluxes are then computed and stored in
     ``Fl`` and ``Fr`` (this might include quantities like cosmic ray
     energy density and flux density OR internal energy for the dual
     energy formalism)
  6. These quantities are all passed to the static public
     ``calc_riemann_fluxes`` method provided by ``ImplStruct``. This
     method then directly updates the fields provided to hold each
     Riemann Flux.

A separate method is provided to compute the fluxes for the passively
advected quantities.
     
*Currently EnzoRiemannImpl has only been tested and known to work for 3D problems. Additionally, no solvers are currently implemented that explicitly support barotropic equations of state, but all of the machinery is in place to support them.*

ImplStruct Class
~~~~~~~~~~~~~~~~

This subsection provides a brief description of the ``ImplStruct`` classes
used to specialize ``EnzoRiemannImpl<ImplStruct>`` to implement specific
Riemann solvers. Basically an ``ImplStruct`` must provide two static public
methods ``calc_riemann_fluxes`` and ``scratch_space_length``.

The ``calc_riemann_fluxes`` static method computes the Riemann Flux at a given
cell interface. The expected function signature should looks like:

.. code-block:: c++

   void ImplStruct::calc_riemann_fluxes
     (const enzo_float flux_l[], const enzo_float flux_r[],
      const enzo_float prim_l[], const enzo_float prim_r[],
      const enzo_float cons_l[], const enzo_float cons_r[],
      const enzo_float pressure_l, const enzo_float pressure_r,
      const EnzoAdvectionFieldLUT lut, const int n_keys,
      const bool barotropic_eos, const enzo_float gamma,
      const enzo_float isothermal_cs,
      const int iz, const int iy, const int ix,
      EFlt3DArray flux_arrays[], enzo_float scratch_space[]);

``flux_l``/ ``flux_r``, ``prim_l``/ ``prim_r``, and ``cons_l``/
``cons_r`` store the left/right interface fluxes values, primitive
quantities, and conserved quantities (they are passed ``Fl``/ ``Fr``,
``wl``/ ``wr``, and ``Ul``/ ``Ur``, respectively). The left and right
reconstructed pressure values are passed as ``pressure_l`` and
``pressure_r``. The ``lut`` maps the names of different quantities to
indices for each array and ``n_keys`` specifies the number of elements
in each array.

``barotropic_eos`` indicates whether the fluid equation of state is
barotropic. If ``true``, then ``isothermal_cs`` is expected to be non-zero and
if ``false``, then ``gamma`` is expected to be positive.

We note that the calculated Riemann Flux for a quantity stored at index ``i``
of the above arrays should be stored at ``flux_arrays[j](iz,iy,ix)``. Finally,
``scratch_space`` serves as a place to temporarily save quantites during the
calculation.

The length of ``scratch_space`` array expected for a given
``ImplStruct`` is calculated by its other required static method
``scratch_space_length``. The function signature for this method is:

.. code-block:: c++

   int ImplStruct::scratch_space_length(const int n_keys);

Here, ``n_keys`` is the number of elements that arrays like ``prim_l``
and ``prim_r`` hold.


Adding new quantites
--------------------

To add support for new actively advected integrable cell-centered
quantities (e.g. cosmic ray energy/flux), the table of cell-centered
quantities (``FIELD_TABLE``) must be updated.

To add support for computing fluxes for such quantities, modifications must be
made to ``EnzoRiemannImpl``. Currently, an abstract base class called for
``EnzoFluxFunctor`` is provided for this purpose. The idea is define a
subclass to be defined for each additional set of flux calculations and then
in then have the factory method, ``EnzoRiemann::construct_riemann``, pass an
array of the relevant functors to ``EnzoRiemannImpl``.

However, the fact that the functors will be pointers will probably
incur overhead. In reality, the better solution might be to hardcode
in the additonal flux calculation functions in some kind of helper method of
``EnzoRiemannImpl``.

Adding new solvers
------------------

New Riemann Solvers can currently be added to this infrastructure in 2 ways.
They can be subclassed from ``EnzoRiemann`` or ``EnzoRiemannImpl<ImplStruct>``
can be specialized. In either case, the ``EnzoRiemann::construct_riemann``
factory method must be modified to return the new solver and
:ref:`using-vlct-riemann-solver` should be updated.

The additional steps for implementing a new Riemann solver by speciallizing
``EnzoRiemannImpl<ImplStruct>`` are as follows:

  1. Define a new ``ImplStruct`` class (e.g. ``HLLDImpl``)

  2. Add the new particlular specialization of ``EnzoRiemannImpl`` to enzo.CI
     (e.g. add the line: ``PUPable EnzoRiemannImpl<HLLDImpl>;``)

  3. *(optional)* define an alias name for the specialization of
     ``EnzoRiemannImpl`` that uses the new ``ImplStruct`` class
     (e.g. ``using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDImpl>;``).


==============================
Updating integrable quantities
==============================

The ``EnzoIntegrableUpdate`` class has been provided to encapsulate
the operation of updating integrable quantities after a (partial)
time-step. It has been factored out of the ``EnzoMethodMHDVlct``
class because this operation will appear in all Godunov solvers.

The constructor for ``EnzoIntegrableUpdate`` has the following
signature:

.. code-block:: c++

   EnzoIntegrableUpdate(std::vector<std::string> integrable_groups,
		        bool skip_B_update,
		        bool dual_energy_formalism,
		        std::vector<std::string> passive_groups)

The requires that we register the names of the integrable
quantities (with ``integrable_groups``), indicate whether the update
to the magnetic field should be skipped, if the dual energy formalism
is in use or if we register the names of the groups
containing the passively advected quantities (with ``passive_groups``)
and requires that we specify the name of the solver ``solver``.

The update to the magnetic field should be skipped when Constrained
Transport is in use (which handles the magnetic field update
separately). The value specified for ``skip_B_update`` is unimportant
if the magnetic field is not specified as an integrable group. Note that
the ``dual_energy_formalism`` argument is purely for demonstration
purposes. It is not yet implemented and if a true value is specified,
then an error will be raised.

Note that the names of the integrable quantites should match the names
specified in ``FIELD_TABLE``; for more details see
:ref:`Centered-Field-Registry`

The main interface function has the signature 

.. code-block:: c++

   void update_quantities(Block *block, Grouping &initial_integrable_group,
			  Grouping &xflux_group, Grouping &yflux_group,
			  Grouping &zflux_group,
			  Grouping &out_integrable_group,
			  Grouping &out_conserved_passive_scalar,
			  EnzoEquationOfState *eos,
			  double dt, int stale_depth);

This function adds the flux divergence (computed from ``xflux_group``,
``yflux_group``, ``zflux_group``) to the values of the both the
actively and passively advected integrable quantities (from
``initial_integrable_group``). The results for the actively advected
quanties are stored in ``out_integrable_group`` and the results for
the passively advected scalars are stored in conserved form in the
fields held by ``out_conserved_passive_scalar`` (note that the initial
values of the passive scalars specified in ``initial_integrable_group``
are in specific form).

*Once source terms need to be added it may make sense to make the
consolidation of the fluxes and source terms into a separate step.*
