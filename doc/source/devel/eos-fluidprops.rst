************************
Fluid Properties and EOS
************************

*[This page is under development]*

Fluid properties are centralized in the :cpp:class:`!EnzoPhysicsFluidProps`.
The user documentation is described :ref:`here <user-fluid_props>`.

As an aside, this section uses terms like "integration" and "primitive" quantities or "stale depth" that are defined in the :ref:`Hydro/MHD Infrastructure section <HydroMHDInfrastructure-page>`.

===========
Fluid props
===========

The primary goal of the :cpp:class:`!EnzoPhysicsFluidProps` class is to hold general data about fluid properties (that may affect multiple parts of the codebase) in a central location and provide methods for accessing this information.
An important guiding principle is that the class is immutable: once it's initialized it should **not** change.

With that said, it does provide some other miscellaneous useful methods.

Primary methods
---------------

All of the primary methods return references to immutable objects that can be used to query information related to the fluid properties.

.. cpp:function:: const EnzoDualEnergyConfig& \
                  EnzoPhysicsFluidProps::dual_energy_config() \
                  const noexcept

  Access a constant reference to the contained dual energy configuration object.

.. cpp:function:: const EnzoFluidFloorConfig& \
                  EnzoPhysicsFluidProps::fluid_floor_config() const noexcept

  Access a constant reference of the object that encodes the fluid floor properties.

.. cpp:function:: const EnzoEOSVariant& EnzoPhysicsFluidProps::eos_variant() \
                  const noexcept

  Access a constant reference to the contained :cpp:class:`!EnzoEOSVariant` class.
  It holds an object that represents the caloric EOS used by the simulation.
  See the :ref:`EOS-Developer-Guide` section for more details about the design of the EOS functionality (and how to use it).

Misc useful methods
-------------------
The functionality described in this subsection are only defined as methods of :cpp:class:`!EnzoPhysicsFluidProps` for a lack of better place to define them.
In the future, it might make sense to move them around.

It's important that none of these functions actually mutate the contents of :cpp:class:`!EnzoPhysicsFluidProps`.
As in the Hydro/MHD Infrastructure section, many of these function operations act on the contents of :cpp:class:`!EnzoEFltArrayMap` rather than directly on :cpp:class:`!Block` objects for additional flexibility.

.. cpp:function:: void EnzoPhysicsFluidProps::primitive_from_integration \
                  (const EnzoEFltArrayMap &integration_map, \
                  EnzoEFltArrayMap &primitive_map, \
                  int stale_depth, \
                  const std::vector<std::string> &passive_list, \
                  bool ignore_grackle = false) const

  This method is responsible for computing the primitive quantities (to be held in ``primitive_map``) from the integration quantities (stored in ``integration_map``).
  Non-passive scalar quantities appearing in both ``integration_map`` and ``primitive_map`` are simply deepcopied and passive scalar quantities are converted from conserved-form to specific form.
  If :class:`!EnzoPhysicsFluidProps` holds a non-barotropic EOS, this method also computes pressure (by calling :cpp:func:`EnzoEquationOfState::pressure_from_integration`).

.. cpp:function:: void EnzoPhysicsFluidProps::pressure_from_integration \
                  (const EnzoEFltArrayMap &integration_map, \
                  const CelloArray<enzo_float, 3> &pressure, \
                  int stale_depth, bool ignore_grackle = false) const

  This method computes the pressure from the integration quantities (stored in ``integration_map``) and stores the result in ``pressure``.
  This wraps the :cpp:class:`!EnzoComputePressure` object whose default behavior is to use the Grackle-supplied routine for computing pressure when the simulation is configured to use :cpp:class:`!EnzoMethodGrackle`.
  The ``ignore_grackle`` parameter can be used to avoid using that routine (the parameter is meaningless if the Grackle routine would not otherwise get used).
  This parameter's primary purpose is to provide the option to suppress the effects of molecular hydrogen on the adiabatic index (when Grackle is configured with ``primordial_chemistry > 1``).

.. cpp:function:: void EnzoPhysicsFluidProps::apply_floor_to_energy_and_sync \
                  (EnzoEFltArrayMap &integration_map, const int stale_depth) \
                  const

   This method applies the pressure floor to the ``"total_energy"`` array specified in ``integration_map``.
   If using the dual-energy formalism the floor is also applied to the ``"internal_energy"`` (also specified in ``integration_map``) and synchronizes the ``"internal_energy"`` with the ``"total_energy"``.

   If :cpp:class:`!EnzoPhysicsFluidProps` holds a barotropic EOS, this method should do nothing.

   .. note::
      In the future, it may make sense to directly pass the pressure floor. 


    .. _EOS-Developer-Guide:

===
EOS
===

Overviews
---------
This section documents how different caloric/isothermal equations of state are supported in Enzo-E.
For reference, these govern the relationship between density, pressure, and internal (or thermal) energy.

Currently, Enzo-E supports relatively few equations of state.
Over time, Hydro codes have a tendency to add support for multiple different types of equations of state.
It’s therefore important to have a solid strategy in place early on to support multiple equations of state.
Unlike other simulation codes (e.g., Athena++) that partially configure physics-features (like the choice of EOS) at compile-time, Enzo-E tends to takes the approach of compiling all physics at once. Thus, Enzo-E needs to support the selection of the EOS at runtime.

The obvious strategy (and the original approach that we took) is to use inheritance with virtual methods.
However, virtual methods are not well-suited for being used within compute kernels (i.e., in the body of a for-loop).
Issues arise because: (i) there is overhead associated with virtual method calls and (ii) there are problems with invoking the virtual methods on GPUs.
While we can get around this to some degree by designing virtual methods to be called outside of the for-loop, there will always be cases where EOS details must be known within a for-loop (e.g., in a Riemann Solver).

For this reason, we choose a different approach for achieving polymorphism, which involves using a `tagged union <https://en.m.wikipedia.org/wiki/Tagged_union>`_.
The idea is that we represent each type of EOS as a stand-alone class and the type-safe union holds an instance of one of those classes (unlike a typical union, the type-safe union explicitly prohibits unsafe access of any union member other than the one that is currently stored)
This is an approach popular in functional programming, in modern languages (e.g., Rust), and that has received support in C++17.
While this approach still incurs some overhead analogous to that of a virtual function, it provides much greater flexibility in choosing where/when we pay this overhead.
For example, we can choose to pay this cost just before a for-loop.

High-Level Design
-----------------
As mentioned above, :cpp:class:`!EnzoEOSVariant`, is simply a container that holds an EOS object.
An EOS object is an instance of one of a few (unrelated, standalone) classes like :cpp:class:`!EnzoEOSIdeal` or :cpp:class:`!EnzoEOSIsothermal`.
that acts as a typesafe union which always holds an instance of one of these classes.

When the :cpp:class:`!EnzoPhysicsFluidProps` physics object is constructed, it creates an instance of the EOS class and holds it internally within an instance of :cpp:class:`!EnzoEOSVariant`.
Throughout the remainder of the simulation, the :cpp:class:`!EnzoPhysicsFluidProps` physics object prevents mutation of the :cpp:class:`!EnzoEOSVariant` instance that it owns (and consequently to the contained EOS object).
Users can access this object with :cpp:func:`EnzoPhysicsFluidProps::eos_variant`.

:cpp:class:`!EnzoEOSVariant` implements a type-safe union that is modelled after the :cpp:class:`!std::variant` template class introduced in C++17.
This class was designed in a way that the internals can easily be replaced with :cpp:class:`!std::variant` when Enzo-E eventually transitions to using C++17.

.. note::

  The choice to model :cpp:class:`!EnzoEOSVariant` after :cpp:class:`!std::variant` leads to slightly more complicated code than is strictly necessary.


    .. _EOSClassDescription-section:

EOS Classes
-----------

These classes are supposed to be lightweight struct/classes that encapsulate an equation of state.
It's also important that these objects are cheap to copy.
They are all entirely defined in header files to facilitate inlining.

We currently expect an EOS class, ``EOSClass``, to provide the following methods:

  * ``constexpr static const char* EOSClass::name() noexcept`` This returns the name of the EOS class (it should match the name a user would specify in a parameter file)

  * ``constexpr static bool EOSClass::is_barotropic() noexcept`` This should return true when the pressure field is just a function of density (e.g., in an isothermal gas).

  * ``std::string EOSClass::debug_string() const noexcept`` This should return a string (for debugging purposes) that represents the internal state of the EOS object.

Other methods supported by an EOS may include calculation of sound speed, fast magnetosonic speed, internal energy, etc.
Essentially all (non-static) methods of an EOS-object are declared as ``const`` (i.e. there's no reason for them to mutate internal state).

One of the perks of using tagged unions is that different types of EOS objects don't NEED to implement the same methods.
For example, it doesn't make much sense for an isothermal eos to support methods that compute the thermal energy.

Currently, two EOS classes exist: :cpp:class:`!EnzoEOSIdeal` and :cpp:class:`!EnzoEOSIsothermal`.
The :cpp:class:`!EnzoEOSIdeal` class implements methods that, given the density and pressure, will compute the following quantities:

  * specific internal energy
  * internal energy density
  * sound speed
  * fast magnetosonic speed (this requires magnetic field values to
    also be specified)

 At the time of writing this section, :cpp:class:`!EnzoEOSIsothermal` is mostly just a placeholder that is used alongside the PPML method (it's not actually used within the PPML method, but it indicates the choice of EOS when other methods are used alongside PPML).

.. note::

  At this time, temperature-related stuff is handled entirely outside of the EOS.
  The rationale for this choice is that this functionality is somewhat unrelated to hydro-solvers (but this is something that can be revisited in the future).

  Grackle has also **NOT** been integrated with the EOS solver at this time.
  (this may need to be revisited in the future).

.. note::

  Currently, to ensure that they are lightweight, all of the EOS classes are "aggregates", which means that they are classes with:

    1. no user-provided or explicit constructors

    2. no private or protected data members (attributes)

    3. no default member initializers (this can be relaxed in C++ 14)

    4. no base classes or virtual methods

  Invariants that might be enforced in a constuctor are instead enforced by a factory method (e.g. :cpp:func:`!EnzoEOSIdealt::construct`). 

  In reality, it would probably simplify the code quite a bit, without sacrificing much/any performance, if we just required that the class was trivially copyable (that's possible without it being an aggregate)

Using ``EnzoEOSVariant`` (accessing stored EOS object)
------------------------------------------------------

The :cpp:class:`!EnzoEOSVariant` class is a type-safe union that ALWAYS holds an instance of one of the types representing an EOS.
EOS objects instances are lightweight structs.

When discussing how to use :cpp:class:`!EnzoEOSVariant`, it is most instructive to describe different operations with examples (rather than providing a detailed API).

Retrieving the EOS Object
~~~~~~~~~~~~~~~~~~~~~~~~~

Let's first imagine we want to write some code that assumes that Enzo-E is configured with an ideal EOS and requires knowledge of the adiabatic index, ``gamma``.
If Enzo-E is configured with a different type of EOS the codebase should terminate with an error.

The following snippet shows a verbose approach for accomplishing this:

.. code-block:: c++

  void my_func(/* args... */) {

    // 1. retrieve pointer to the PE's EnzoPhysicsFluidProps instance
    const EnzoPhysicsFluidProps* fluid_props = enzo::fluid_props();

    // 2. fetch a const reference to the EnzoEOSVariant instance held within
    //    the object pointed to by fluid_props
    const EnzoEOSVariant& eos_variant = fluid_props->eos_variant();

    // 3. fetch a const reference to the eos within eos_variant, while
    //    enforcing the assumption that it's an EnzoEOSIdeal instance
    const EnzoEOSIdeal& eos = eos_variant.get<EnzoEOSIdeal>();

    // fetch the value of gamma
    enzo_float gamma = eos.get_gamma();

    // do work with gamma...
  }

Now, let's break this down in slightly more detail.

  1. :cpp:func:`!enzo::fluid_props()` returns a pointer to the instance of the :cpp:class:`!EnzoPhysicsFluidProps` that is configured for the Processing Element (PE).
     This pointer can't be a ``nullptr`` (if it is, the function will abort with an error).

  2. fetch a const reference to the :cpp:class:`!EnzoEOSVariant` instance held within the object pointed to by ``fluid_props``

  3. fetch a const reference to the eos within ``eos_variant`` if it currently holds an :cpp:class:`!EnzoEOSIdeal`.
     In other cases, the program aborts with an error.

We can write a much more concise form of the above function:

.. code-block:: c++

  void my_func(/* args... */) {

    // the program aborts with an error if Enzo-E was not configured with an
    // ideal eos
    const EnzoEOSIdeal& eos = enzo::fluid_props()->eos_variant().get<EnzoEOSIdeal>();
    enzo_float gamma = eos.get_gamma();
    // do work with gamma...
  }

In both of these snippets we make use of the method:

.. cpp:function:: template<typename T> \
                  const T& EnzoEOSVariant::get() const

  Accessor method that returns a reference to the contained EOS object if ``this`` currently holds the EOS object of type ``T``.
  Otherwise, the program aborts with an error message.
  A non-``const``-qualified version of this method also exists.

  This is a counterpart of the ``std::get`` template function.

Retrieving the EOS Object with Detailed Error Message
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now let's consider a variation on the last case.
In this situation let's imagine that we want to write a more detailed error message in the case where it is executed and Enzo-E is not configured with an ideal EOS:

.. code-block:: c++

  void my_func(/* args... */) {

    const EnzoEOSIdeal* eos
      = enzo::fluid_props()->eos_variant().get_if<EnzoEOSIdeal>();

    if (eos == nullptr) {
      ERROR("my_func",
            "my_func only works when Enzo-E is configured with an ideal EOS");
    }
    enzo_float gamma = eos->get_gamma();
    // do work with gamma...
  }

This snippet makes use of

.. cpp:function:: template<typename T> \
                  const T* EnzoEOSVariant::get_if() const

  Accessor method that returns a pointer to the contained EOS object, if ``this`` curently holds the EOS object of type ``T``.
  Otherwise, a ``nullptr`` is returned.
  The program aborts with an error if ``T`` is a type that :cpp:class:`!EnzoEOSVariant` is incapable of holding.
  A non-``const``-qualified version of this method also exists.

  This is a counterpart of the ``std::get_if`` template function.

Alternatively we could also accomplish the above by writing:

.. code-block:: c++

  void my_func(/* args... */) {

    const EnzoEOSVariant& eos_variant = enzo::fluid_props()->eos_variant();
    if (!eos_variant.holds_alternative<EnzoEOSIdeal>()) {
      ERROR("my_func",
            "my_func only works when Enzo-E is configured with an ideal EOS");
    }
    enzo_float gamma = eos_variant().get<EnzoEOSIdeal>().get_gamma();
    // do work with gamma...
  }

This last snappet employs the following method:

.. cpp:function:: template<typename T> \
                  T* EnzoEOSVariant::holds_alternative() const

  Returns whether ``this`` currently holds the alternative EOS type, ``T``.
  The program aborts with an error, if ``T`` is a type that :cpp:class:`!EnzoEOSVariant` is incapable of holding.

  This acts as a backport for one of C++17's ``std::holds_alternative``

Using ``EnzoEOSVariant`` (General semantics)
--------------------------------------------

The :cpp:class:`!EnzoEOSVariant` class has semantics just like :cpp:class:`!std::variant` (albeit, slightly more limited).

For example, :cpp:class:`!EnzoEOSVariant` is never empty.
If you call:

.. code-block:: c++

  EnzoEOSVariant my_eos_variant;

Then the variable ``my_eos_variant`` holds a default-constructed instance of :cpp:class:`!EnzoEOSVariant`.
At the time of writing this documentation this object will contain an instance of an :cpp:class:`!EnzoEOSIsothermal`, but that is an implementation detail that may change over time.

Like :cpp:class:`!std::variant`, :cpp:class:`!EnzoEOSVariant` also has value-like semantics.
This means that any time you perform a copy on an instance of :cpp:class:`!EnzoEOSVariant` it's a deepcopy.

.. code-block:: c++

  const EnzoEOSVariant& eos_variant = enzo::fluid_props()->eos_variant();

  // make a copy of eos_variant
  EnzoEOSVariant my_eos_variant = eos_variant;

  
Any mutations to the contents of ``my_eos_variant`` will not affect the contents of ``enzo::fluid_props()->eos_variant()``.
Examples might include:

  * changing the type of object stored within ``my_eos_variant`` (if it initially stores an instance of :cpp:class:`!EnzoEOSIsothermal`, we could replace it with an instance of :cpp:class:`!EnzoEOSIdeal`)

  * mutating the attributes of a stored object within ``my_eos_variant`` (one could imagine mutating the value of gamma stored within a :cpp:class:`!EnzoEOSIdeal` instance)

**As an aside**, the API of :cpp:class:`!EnzoPhysicsFluidProps` is designed so that the user can't accidentally mutate the PE's EOS object (you can only mutate copies of that object).

    .. _eos-timestep-example:

Using ``EnzoEOSVariant`` (A concrete example)
---------------------------------------------

Many hydro methods need to determine the maximum timestep that they allow.
In the process, they may need to compute:

.. math::

  C_0 \times \min\left(\frac{\Delta x}{|c_{s,ijk}+ v_{x, ijk}|},
                       \frac{\Delta y}{|c_{s,ijk}+ v_{y, ijk}|},
                       \frac{\Delta z}{|c_{s,ijk}+ v_{z, ijk}|}\right)

where :math:`C_0` is the courant factor (a constant between 0 and 1) and :math:`\Delta x,\, \Delta y,\, \Delta z` specify cell widths.

The following code snippet shows a somewhat simplified example of how you might perform this calculation.
This snippet will abort with an error if each Processing Element's global :cpp:class:`!EnzoPhysicsFluidProps` instance was configured to hold anything other than an ideal EOS.


.. code-block:: c++

   double timestep(CelloArray<const enzo_float, 3> density,
                   CelloArray<const enzo_float, 3> velocity_x,
                   CelloArray<const enzo_float, 3> velocity_y,
                   CelloArray<const enzo_float, 3> velocity_z,
                   CelloArray<const enzo_float, 3> pressure,
                   double dx, double dy, double dz,
                   double courant_factor)
   {
     // the program aborts with an error if Enzo-E was not configured with an
     // ideal EOS (as an aside, we are technically making a copy of the EOS
     // here - that should be totally fine sinces it's just a memcpy)
     const EnzoEOSIdeal eos = enzo::fluid_props()->eos_variant().get<EnzoEOSIdeal>();

     const int mx = density.shape(2);
     const int my = density.shape(1);
     const int mz = density.shape(0);

     double dt = std::numeric_limits<double>::max();
     for (int iz = 0; iz < mz; iz++) {
       for (int iy = 0; iy < my; iy++) {
         for (int ix = 0; ix < mx; ix++) {

           double cs = (double) eos.sound_speed(density(iz,iy,ix),
                                                pressure(iz,iy,ix));
           double abs_vx = std::fabs((double) velocity_x(iz,iy,ix));
           double abs_vy = std::fabs((double) velocity_y(iz,iy,ix));
           double abs_vz = std::fabs((double) velocity_z(iz,iy,ix));
           double tmp = std::min(std::min(dx/(abs_vx + cs),
                                          dy/(abs_vy + cs)),
                                 dz/(abs_vz + cs));
           dt = std::min(dt, tmp);
         }
       }
     }

   return courant_factor * dt;
  }


Using ``EnzoEOSVariant`` (visitor pattern)
------------------------------------------

The :cpp:func:`EOSVariant::visit` method can be used to dispatch code based on the type of the EOS that is stored within the EOSVariant. 
This method effectively implements the `visitor design pattern <https://en.m.wikipedia.org/wiki/Visitor_pattern>`_.
While this is generally most helpful when you have a collection of objects, it may be helpful in simplifying some code in Enzo-E.

The method that is used to accomplish this is defined below, but it's most useful to consider example cases


.. cpp:function:: template<class Visitor> \
                  EnzoEOSVariant::visit(Visitor&& vis) const noexcept

  invokes the callable visitor, ``vis``, by passing the EOS instance held by ``this``.
  The visitor must accept any of the EOS variants passed as an argument, by value, and return an output with a consistent type for all of them.

  This acts like a very crude port of :cpp:func:`!std::visit` from C++17.

Query whether the EOS is barotropic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's consider an example where we want to query whether the EOS is barotropic.
We will make use of the ``is_barotropic`` method that is defined for all EOS classes.

Now the obvious way to write this is:

.. code-block:: c++

   bool is_barotropic_eos(const EOSVariant& eos_variant) {

     if (eos_variant.holds_alternative<EnzoEOSIdeal>()) {
       return eos_variant.get<EnzoEOSIdeal>().is_barotropic();
     } else if (eos_variant.holds_alternative<EnzoEOSIsothermal>()) {
       return eos_variant.get<EnzoEOSIsothermal>().is_barotropic();
     } else {
       ERROR("is_barotropic_eos", "eos_variant holds an unknown eos");
     }
   }

The code snippet shown above is fine, but it could get tedious to have to modify that code every time that we introduce a new type of EOS.
Instead we can write the following snippet, which accomplishes the same thing (but without the caveat):

.. code-block:: c++

   struct IsBarotropicVisitor {
     template <typename T>
     bool operator()(T eos) const { return T::is_barotropic(); }
   };

   bool is_barotropic_eos(const EOSVariant& eos_variant) {
     return eos_variant.visit(IsBarotropicVisitor());
   }

This second snippet is still a little verbose.
It can further simplify in C++14 to:

.. code-block:: c++

   bool is_barotropic_eos(const EOSVariant& eos_variant) {
     return eos_variant_.visit([](auto eos) { return eos.is_barotropic(); });
   }

One could imagine that this example generalizes to any case where all EOS classes provide a common interface (e.g., calling the ``name`` method or the ``debug_string`` method).

More Sophisticated cases
~~~~~~~~~~~~~~~~~~~~~~~~
One could also apply the visitor design pattern in more sophisticated cases, like our :ref:`timestep-example <eos-timestep-example>`.

.. note::

   It’s a little unclear how well this visitor design pattern works with compute kernels.
   At the end of the day, it may make sense to just drop the ``visit`` method.
   (The method's complexity may not be worthwhile)


How to extend this machinery
----------------------------

When you introduce a new EOS class, you need to do three things:

1. You need to update a small subsection of the declaration of :cpp:class:`!EnzoEOSVariant` where the names of the EOS classes are listed.

2. You need to update the :cpp:func:`!pup` routine implemented in the source file for :cpp:class:`!EnzoEOSVariant`.

3. You need to update the :cpp:func:`!EnzoConfig::read_physics_fluid_props_` method to allow the user to specify a new type of EOS.
