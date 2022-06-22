************
Enzo-E Units
************

In non-cosmological simulations, the user is free to specify length, time, and either density
or mass units (only one can be set).
This is done by setting values for ``Units:length``, ``Units:time``, and
either ``Units:density`` or ``Units:mass``, which correspond to the unit
length / time / density / mass in cgs units. If running with gravity, and if the user wants to use
the standard value for the gravitational constant, the user must set a
value for ``Method:gravity:grav_const`` which is consistent with their choice of units; i.e.,
its value must be :math:`G_{cgs}\times M \times T^2 \times L^{-3}`, or equivalently,
:math:`G_{cgs}\times D \times T^2`, where :math:`M, D, T, L` are the mass, density, time, and length
units, and :math:`G_{cgs}` is the value of the gravitational constant in cgs units.

In cosmological simulations, the code ignores any specified units and instead operates in a
coordinate system which is comoving with the universal expansion, defining the length, time,
velocity, and density units as given below (length and density units depend on time / redshift.)

The length unit is specified by ``Physics:cosmology:comoving_box_size``, which gives the length
unit in terms of comoving :math:`Mpc/h`.

The density unit is defined so that the comoving mean matter density of the universe is 1, where
the mean comoving matter density is given by :math:`\frac{3 H_0^2 \Omega_m}{8 \pi G}`.

The time unit is defined so that :math:`\frac{3}{2} H_0^2 \Omega_m (1+z_i)^3 = 1`, where
:math:`z_i` is the initial redshift of the simulation. This is the free-fall time at
:math:`z = z_i`, which has the effect of simplifying Poisson's equation.

The velocity unit is defined as :math:`\frac{1+z_i}{1+z} L / T`, where :math:`L` is the length
unit, and :math:`T` is the time unit.

For cosmological simulations, the value set for ``Method:gravity:grav_const`` is ignored.

In all simulations, the ``"temperature"`` field always has units of Kelvin.
