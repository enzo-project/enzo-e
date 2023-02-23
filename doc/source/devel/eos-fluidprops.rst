************************
Fluid Properties and EOS
************************

*[This page is under development]*

Fluid properties are centralized in the :cpp:class:`!EnzoPhysicsFluidProps`.
The user documentation is described :ref:`here <user-fluid_props>`.

As an aside, this section uses terms like "integration" and "primitive" quantities or "stale depth" that are defined in the Hydro/MHD Infrastructure section.

===========
fluid props
===========

The primary goal of the :cpp:class:`!EnzoPhysicsFluidProps` class is to hold general data about fluid properties (that may affect multiple parts of the codebase) in a central location and provide methods for accessing this information.
An important guiding principle is that the class is immutable: once it's initialized it should **not** change.

With that said, it does provide some other miscellaneous useful methods.

primary methods
---------------
*[To be filled in]*

misc useful methods
-------------------
The functionallity described in this subsection are only defined as methods of :cpp:class:`!EnzoPhysicsFluidProps` for a lack of better place to define them.
In the future, it might make sense to move them around.

It's important that none of these functions actually mutate the contents of :cpp:class:`!EnzoPhysicsFluidProps`.
As in the Hydro/MHD Infrastructure section, many of these functions operations act on the contents of :cpp:class:`!EnzoEFltArrayMap` rather than directly on :cpp:class:`!Block` objects for additional flexibility.

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

===
EOS
===
*[To be filled in]*
