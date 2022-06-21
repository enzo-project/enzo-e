************
Enzo-E Units
************

In non-cosmological simulations, the user is free to specify length, time, and either density
or mass units. This is done by setting values for ``Units:length``, ``Units:time``, and
``Units:density`` or ``Units:mass``, which correspond to the unit length / time / density / mass
in cgs units. If running with gravity, the user must set a value for
``Method:gravity:grav_const`` which is consistent with their choice of units.

In cosmological simulations, the units are defined by physical constants and cosmological
parameters.

The length unit is specified by ``Physics:cosmology:comoving_box_size``, which gives the length
unit in terms of comoving :math:`Mpc/h`.

The density unit is defined so that the comoving mean matter density of the universe is 1, where
the mean comoving matter density is given by :math:`\frac{3 H_0^2 \Omega_m}{8 \pi G}`.

The time unit is defined so that :math:`\frac{3}{2} H_0^2 \Omega_m (1+z_i)^3 = 1`, where
:math:`z_i` is the initial redshift of the simulation.

For cosmological simulations, the value set for ``Method:gravity:grav_const`` is ignored.




