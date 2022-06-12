.. _about-physics-groups:

****************************
Enzo-E Physics Config-Groups
****************************

*[ This page is under development ]*

This section decribes the Physics subgroups that can be specified in
the ``Physics`` Group of the parameter file. Each subgroup directly
maps to a different type of ``Physics`` object in the codebase. These
objects hold information that needs to be accessible across different
Enzo-E methods.

``"cosmology"``
===============

Specifies cosmological parameters

``"fluid_props"``
=================

Specifies parameters related to the fluid properties. These parameters
are largely into a few subcategories:

  1. ``dual_energy``: these parameters govern whether the simulation
     uses the dual-energy formalism and (if applicable) which
     formulation is used.

  2. ``floors``: these parameters specify floors for used for a small
     assortment of quantities.

  3. ``eos``: specifies parameters related to the (caloric) equation of
     state.

.. _using-fluid_props-de:

dual-energy formalism
---------------------

While the dual-energy formalism is inherently tied to the user's
choice of hydrodynamics solver, the choice to use the dual-energy
formalism can have implications for other methods and for initial
conditions. This is the primary reason why the dual-energy
formalism properties are tracked as part of the ``fluid_props``
physics object.

.. list-table:: Method ``fluid_props:dual_energy`` parameters
   :widths: 10 5 1 30
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``"type"``
     - `string`
     - ``"disabled"``
     - ``Specifies handling of magnetic fields (or lack thereof)``
   * - ``"eta"``
     - `list(float)`
     - `-`
     - ``Specifies handling of magnetic fields (or lack thereof)``



Two separate formulations of the dual-energy formalism are nominally
supported.

The ``"ppm"`` method uses the original formulation described in
`Bryan et al (1995)
<https://ui.adsabs.harvard.edu/abs/1995CoPhC..89..149B>`_

ADD MORE
