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
