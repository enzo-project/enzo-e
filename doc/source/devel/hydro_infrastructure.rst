****************************
Hydro/MHD C++ Infrastructure
****************************

In this section we discuss some of the C++ infrastructure provided in the Enzo
Layer that was used to implement the VL + CT MHD solver that can be optionally
re-used to implement other hydro/MHD methods


===============
Shorthand Terms
===============

Here we briefly define a few terms that we defined and used throughout the
documentation and codebase

stale depth
-----------

Benefits:

  * Simplifies calculation of correct number of ghost zones for unigrid
    simulations

  * When used alongside ``CelloArray``, it drastically simplifies the
    calculation of which indices should be iterated over. The
    ``EnzoFieldArrayFactory`` can take the stale depth as an argument in its
    constructor. Then any array it builds has the stale depth clipped off.
    This allows you to write all for loop logic as though you are using
    nearest-neighbor interpolation and ignore any preceeding partial
    timesteps.


Integrable/Reconstructable Quantities
-------------------------------------

Throughout this guide and the relevant code we refer to quantities by the
terms "Integrable" and "Reconstructable".


==============
General Design
==============

Explain and motivate the use of Groupings everywhere (

Explain the idea of Registration of quantities within each "fixture"

Dimension invariance


=======================
Centered Field Registry
=======================

Talk about Centered Field Registry and the ``FIELD_TABLE``

Talk about ``EnzoAdvectionFieldLUT`` - compare it to alternatives

EnzoAdvectionFieldLUT
---------------------

Motivation
~~~~~~~~~~

*Cleanup*

  * The calculation of fluxes and wave speeds requires random access of
    specific fields. Having the ability to iterate over the entries
    of multiple arrays simultaneously is also very convenient.

      * Unlike Enzo (and many other codes) we wanted to avoid statically
	declaring which indices correspond to which fields. Doing so
	complicates situations when fields can be optionally added
	(e.g. like optional magnetic fields, dropping total energy for
	barotropics equations of state, optionally advecting internal
	energy for dual energy formalism or tracking cosmic ray energy density
	and cosmic ray energy fluxes)

      * Instead the ordering of fields are determined at registration. We cam
	sort the values by whether they are conserved/specific (and iterate
	over all conserved or all specific primitive quantities), which allows
	simplifies the challenges of having optional fields.

  * We settled on using the field_lut struct as a lookup table. The
    struct has members named for every quantity listed in FIELD_TABLE

      * For a ``SCALAR``, the member name directly matches the name in
        column 1

      * For a VECTOR, there are 3 members: {name}_i, {name}_j, {name}_k
        ({name} corresponds to the name appearing in column 1)

    Each struct contains members named advection-related (related to
    riemann fluxes) quantities in the table.

  * ``EnzoCenteredFieldRegistry::prepare_advection_lut`` is used to
    prepare the lookup table for pre-specified quantities and determines
    the length of an array large enough to hold fields representing all
    of the quantities. The function then initializes all members of the
    struct to be equal to quantities from 0 through 1 less than the
    length of an arrays and organizes them based on whether the quanitiy
    is conserved, specific, or other. The function also yields the
    indices required to just iterate over each category of field. All
    members of the struct that don't correspond to a specified quantity
    are set to -1

  * ``EnzoCenteredFieldRegistry::load_array_of_fields`` constructs an array of
    instances of ``EFlt3DArray`` where the ordering of the fields is
    determined by the supplied instance of ``EnzoAdvectionFieldLUT``


Example Usage
~~~~~~~~~~~~~

  * Example: If we have an array of reconstructed integrable primitives,
    ``wl``, and an initialized instance of ``EnzoAdvectionFieldLUT``, ``lut``,
    then ``wl[prim_lut.density]`` and ``wl[prim_lut.total_energy]`` indicates
    the entries reserved for the density and (specific) total energy

=================
Equation Of State
=================

talk about ``EnzoEquationOfState`` and ``EnzoEOSIdeal``

Nominally responsible for

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
be retireved using the ``EnzoReconstructor::construct_reconstructor`` static
factory method.

Public Interface
----------------
The main interface function provided by this class is:

reconstruct_interface (Block *block, Grouping &prim_group,
				      Grouping &priml_group,
				      Grouping &primr_group, int dim,
				      EnzoEquationOfState *eos,
				      int stale_depth)

How to extend
-------------
To add a new reconstructor, subclass ``EnzoReconstructor`` and provide
definitions for the virtual methods. Additionally, modify the factory method
``EnzoReconstructor::construct_reconstructor`` to return pointers to
instances of the new class when the appropriate name is passed as an argument.

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
and requires that we specify the name of the solver ``solver``. We note that
the names of the integrable quantites should match the names specified in
``FIELD_TABLE`` that are at the start of specified in
src/Enzo/enzo_EnzoCenteredFieldRegistry.hpp

To actually use the Riemann solver, the virtual solve method should be called.
The signature for this method is

.. code-block:: c++

   void solve (Block *block, Grouping &priml_group, Grouping &primr_group, 
               std::string pressure_name_l, std::string pressure_name_r,
               Grouping &flux_group, int dim, EnzoEquationOfState *eos,
               int stale_depth);

In this function, the ``priml_group`` and ``primr_group`` arguments are
references to ``Grouping`` objects that have groups named for each item in the
``integrable_group`` and ``passive_groups`` arguments originally passed to the
factory method. Within each group, the ``Grouping`` objects should contain the
field names holding the left and right reconstructed values that represent a
quantity. For ``SCALAR`` quantites held in integrable_group, a group should only have 1 field while for ``VECTOR`` quantites, a group should have 3 fields (labelled for the x, y, and z components.

Explain more




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
advectable, integrable quantities and the members corresponding to
the registered integrable quantities assigned values that
correspond to indices in an array. See ? for a more detailed description
of the struct. This struct is used inplace of the more tradional
approach of defining global enums or macros to map quantity names to
array indices. Note that we when we tried to use a built-in hash tables to
perform this task, had a significant performance cost.

Below, a brief overview of the ``EnzoRiemannImpl::solve`` control flow is
provided. Basically the function loops over all cell interfaces, along a given
dimension, where the flux should be computed. At each location, the following
sequence of operations are performed:

  1. Retrieve the left and right primitives at the given location from the
     fields and stores them in arrays of ``enzo_float`` elements called
     ``wl`` and ``wr``. The elements are ordered based on a preconfigured
     instance of ``EnzoAdvectionFieldLUT`` called ``lut_``.
  2. The left and right pressure values are retrieved from the temporary
     fields holding the values that were precomputed from the reconstructed
     quantities (presumably using a concrete subclass of
     ``EnzoEquationOfState``). The values are stored in ``pressure_l`` and
     ``pressure_r``.
  3. The conserved forms of the left and right reconstructed primitives and
     stored in the arrays called ``Ul`` and ``Ur``. Note that the primitives
     that are always in conserved form (e.g. density or magnetic field) are
     simply copied over. The elements of ``Ul`` / ``Ur`` maintain the same
     ordering as those of ``wl`` / ``wr`` (e.g. the index for a given component
     of the velocity in ``wl`` / ``wr`` is the index for the same component of
     the momentum in ``Ul`` / ``Ur``).
  4. The standard left and right hydro/mhd fluxes are computed using
     the above quantities and stored in ``Fl`` and ``Fr``; the elements are
     again ordered by ``lut_``.
  5. In principle, non-standard fluxes are then computed and stored in ``Fl``
     and ``Fr`` (this might include quantities like cosmic ray energy density
     and flux density OR internal energy for the dual energy formalism)
  6. These quantities are all passed to the static public
     ``calc_riemann_fluxes`` method provided by ``ImplStruct``. This method
     then directly updates the fields provided to hold each Riemann Flux.

A separate method is provided to compute the fluxes for the passively advected
quantities.
     
*Currently EnzoRiemannImpl has only been tested and known to work for 3D
 problems. Additionally, no solvers are currently implemented that explicitly
 support barotropic equations of state, but all of the machinery is in place
 to support them.*

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

*why its factored out* (to provide easy reuse)
