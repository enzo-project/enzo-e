.. _about-physics-groups:

****************************
Enzo-E Physics Config-Groups
****************************

*[ This page is under development ]*

This section decribes the Physics subgroups that can be specified in
the ``Physics`` Group of the parameter file. Each subgroup directly
maps to a different type of ``Physics`` object in the codebase. These
objects hold information that needs to be accessible across different
Enzo-E methods and/or initializers.

In a parameter file, a user currently needs to explicitly list the
names of all of the physics objects that they are configuring within
the :par:param:`Physics:list` parameter.

.. note::

   Currently, the ``"fluid_props"`` and ``"gravity"`` don't need to be
   :par:param:`Physics:list` for Enzo-E to parse and make use of
   parameters in the respective Physics subgroups. This choice is made
   purely for the sake of backwards compatability (so that parameter
   files designed for earlier versions of the code will continue to
   work properly).

   With that said, this is mostly just an implementation detail (that
   is subject to change, especially as deprecated parameters are
   eventually removed). At this time, Users should still explicitly
   list these physics groups within the :par:param:`Physics:list`
   parameter if they want to use/configure them.

``"cosmology"``
===============

Specifies cosmological parameters

.. _user-fluid_props:

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

A typical Godunov scheme might carry the specific total energy, :math:`E`, as one of its primary variables (the ``"total_energy"`` field).
To compute the pressure (e.g. for reconstruction), the scheme needs to first compute the specific internal energy, :math:`e=p/((\gamma - 1) \rho)`, by subtracting the non-thermal energy components from :math:`E`.
Numerical problems arise under extreme conditions when :math:`E_{\rm non-thermal}/e` is extremely large (e.g. from some combination of high Mach number or magnetic energy) since :math:`e` is the difference between 2 large numbers, :math:`E-E_{\rm non-thermal}`.

The dual-energy formalism is a technique used to avoid these numerical issues.
Under this technique we track an additional field called ``"internal_energy"`` that we frequently synchronize with with the ``"total_energy"`` field.
Under these extreme conditions, then use the ``"internal_energy"`` field to provide extra precision for :math:`e`.

In some sense the the dual-energy formalism can be considered an implementation detail of the hydro/MHD solver.
However, the choice whether to use the dual-energy formalism has important implications for other methods and for initializers.
This is the primary reason why the dual-energy formalism properties are tracked as part of the ``fluid_props`` physics object.

The relevant parameters are listed below:

.. list-table:: Method ``fluid_props:dual_energy`` parameters
   :widths: 10 5 1 30
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - :par:param:`~Physics:fluid_props:dual_energy:type`
     - `string`
     - ``"disabled"``
     - choice of formalism: ``"disabled"``, ``"bryan95"``, ``"modern"``
   * - :par:param:`~Physics:fluid_props:dual_energy:eta`
     - `list(float)`
     - `-`
     - Interpretation (and defaults) depend on value of :par:param:`~Physics:fluid_props:dual_energy:type`

When :par:param:`~Physics:fluid_props:dual_energy:type` is ``"disabled"``, the simulation runs without the dual-energy formalism
This is the default configuration and in this case, the :par:param:`~Physics:fluid_props:dual_energy:eta` parameter should **NOT** be specified.

The other two choices for :par:param:`~Physics:fluid_props:dual_energy:type` refer to two variants of the dual-energy formalism that we describe below.

``"bryan95"`` variant
~~~~~~~~~~~~~~~~~~~~~
This is the original formulation of the dual-energy formalism described in
`Bryan et al (1995)
<https://ui.adsabs.harvard.edu/abs/1995CoPhC..89..149B>`_.
It is used by the ``"ppm"`` solver and it is parameterized by two values: :math:`\eta_1\, \&\, \eta_2`.

:par:param:`Physics:fluid_props:dual_energy:eta` expects a list of 2 entries: :math:`[\eta_1, \eta_2]` (users are **NOT** permitted to provide a single entry).
When this parameter isn't specified, it defaults to ``[0.001, 0.1]``.

``"modern"`` variant
~~~~~~~~~~~~~~~~~~~~
This implementation is used by the ``"mhd_vlct"`` solver and it more closely resembles the implementation employed in Enzo's Runge–Kutta integrator than the  ``"bryan95"`` variant.
This variant is parameterized by a single value: :math:`\eta`, and thus :par:param:`Physics:fluid_props:dual_energy:eta` should only provide a single entry.

There are 3 primary differences from the ``"bryan95"`` variant:

  1. the ``"internal_energy"`` field is always used to compute pressure. Under the ``bryan95`` variant, pressure could be computed from ``"total_energy"`` or ``"internal_energy"`` (the decision was independent of synchronization).
  2. pressure and ``"internal_energy"`` are not separately reconstructed. Instead, just the pressure is reconstructed. The ``"internal_energy"`` is computed at the left and right interfaces from the reconstructed quantities.
  3. Synchronization of the total and internal energies is a local operation that doesn't require knowledge of cell neighbors. The ``"bryan95"`` variant requires knowledge of the immediate neighbors (each synchronization incremented the stale depth - so 3 extra ghost zones would have been required for the ``"mhd_vlct"`` solver).

For clarity, the conditions for synchronization are provided below. The specific ``internal_energy``, :math:`e`, is set to :math:`e'= E - (v^2 + B^2/\rho)/2` (where :math:`E` is the specific ``total_energy``) when the following conditions are met:

  * :math:`c_s'^2 > \eta v^2`, where :math:`c_s'^2=\gamma(\gamma - 1) e'`.
  * :math:`c_s'^2 > \eta B^2/\rho` (this is always satisfied in hydro mode)
  * :math:`e' > e /2`

If the above conditions are not met, then ``total_energy`` is set to :math:`e + (v^2 + B^2/\rho)/2` in MHD mode (in hydro mode, it's set to :math:`e + v^2/2`).

*Note: in the future, the behavior described in difference 2, may change
to achieve better compatibility with Grackle.*

AMR Inconsistency
~~~~~~~~~~~~~~~~~

A minor inconsistency is present in our simulation when we use either dual-energy formalism in AMR simulations.
Each formulation of the dual-energy formalism requires the hydro-solver to add a source term to the ``"internal_energy"`` field.
This source term involves the velocity on each cell face (just the component normal to the face).
In practice, we compute the velocity in the Riemann Solver.

When neighboring cells compute the source term, they should theoretically use consistent values for the velocity component on the shared face.
While this is true for neighbors with the same refinement level, it’s not true for neighbors with different refinement levels.
In principle, a coarse cell should actually compute the source term from some kind of average of the velocity component on the faces of the neighboring finer cells (this could hypothetically be achieved with the equivalent of a "flux-correction").

We don’t expect this inconsistency to be important given that the dual-energy formalism only affects highly supersonic flows and isn’t conservative anyway.

EOS
---

The ``"fluid_props:eos"`` subgroup holds parameters that configure the nominal (caloric) equation of state.
These parameters primarily affect the Hydro/MHD methods and the ``"grackle"`` method.
It also affects the calculation of pressure and temperature fields.

At this time, all simulations are assumed to have an ideal EOS and the only configurable parameter is provided below.
In the future, further EOS customization will be supported in this section

.. list-table:: Method ``fluid_props:eos`` parameters
   :widths: 10 5 1 30
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - :par:param:`~Physics:fluid_props:eos:gamma`
     - `float`
     - ``5.0/3.0``
     - Adiabatic index (a.k.a. the ratio of specific heats)

See :ref:`using-grackle-gamma-with-HD` for further discussion about
how the equation of state is handled when
:par:param:`Method:grackle:primordial_chemistry` exceeds ``1`` (under
these conditions Grackle models a spatially varying adiabatic index).

.. _using-fluid_props-floors:

Floors
------

The ``"fluid_props:floors"`` subsection is used for specifying the floors of different fluid quantities. A list of the quantities whose floors can be configured are provided below.

.. list-table:: Method ``fluid_props:floors`` parameters
   :widths: 10 5 1 30
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - :par:param:`~Physics:fluid_props:floors:density`
     - `float`
     - `-`
     - The floor to apply to the ``"density"`` field.
   * - :par:param:`~Physics:fluid_props:floors:pressure`
     - `float`
     - `-`
     - The floor to apply to the ``"pressure"`` field.
   * - :par:param:`~Physics:fluid_props:floors:temperature`
     - `float`
     - `-`
     - The floor to apply to the ``"temperature"`` field.
   * - :par:param:`~Physics:fluid_props:floors:metallicity`
     - `float`
     - `-`
     - This multiplied by the ``"density"`` field and ``enzo_constants::metallicity_solar`` gives the floor for the ``"metal_density"`` field.

See :ref:`using-methods` for discussions of the floors that are actually used by a given method.
Be mindful that unlike the other parameters, the :par:param:`~Physics:fluid_props:floors:metallicity` parameter doesn't directly specify the floor for a fluid field (the actual floor depends on other quantities).

.. note::

   The ``"pressure"`` and ``"temperature"`` fields can be written to disk as derived quantities (if the fields are specified in the "derived" grouping).
   In these cases, these quantities are computed using ``EnzoComputePressure`` and ``EnzoComputeTemperature``, respectively.
   You may want to check these classes to see if/when the floors get applied.

``"gravity"``
=============

Specifies the gravitational constant. In the future, additional
gravity-related parameters could be introduced.

.. list-table:: Physics ``gravity`` parameters
   :widths: 10 5 1 30
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - :par:param:`~Physics:gravity:grav_const_codeU`
     - `float`
     - `-`
     - The gravitational constant specified in code units. When not
       specified, it's automatically computed from the real-world
       reference value :math:`G\approx 6.67\times 10^{-8}\, {\rm
       cm}^3\, {\rm g}^{-1} {\rm s}^{-2}` (see codebase for exact
       value).

In most cases, users should not specify
:par:param:`~Physics:gravity:grav_const_codeU` at all (so that the
appropriate default value is used). This parameter mostly exists to
help simplify some test problems in non-cosmological simulations.

Users are **NOT** allowed to specify
:par:param:`~Physics:gravity:grav_const_codeU` parameter in
cosmological simulations. This is because cosmological code-units are
defined such that :math:`4\pi G\bar{\rho}` has the value ``1.0``,
where :math:`\bar{\rho}` is the mean physical matter density of the
universe.

We generally advise users to include ``"gravity"`` within
:par:param:`Physics:list` whenever they use any method involving
gravity, event if they aren't explicitly assigning a value to
:par:param:`Physics:gravity:grav_const_codeU`.

.. note::

   At the time of writing, if the user specifies both ``"cosmology"``
   and ``"gravity"`` within :par:param:`Physics:list`, it's important
   that ``"cosmology"`` comes first. In the future, we can hopefully
   relax these requirements.
