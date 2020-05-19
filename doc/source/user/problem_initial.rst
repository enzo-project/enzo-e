*************************
Enzo-E initial conditions
*************************

*[ This page is under development ]*

Initial conditions define the initial setup of a problem.  They
include field values and particle values at the start of the
simulation.  They can be declared in the parameter file using Cello's
"value" initializer, or specific problem initializers can be invoked,
such as "implosion_2d" for the "Implosion test", or "sedov_array_3d" for
a 3D array of Sedov blast waves.

``"cloud"``
   Initialize a spherical cloud embedded in a hot wind.

``"collapse"``
   Initialize a spherical collapse test.

``"file"``
   Initialize from input HDF5 file.  Not implemented yet.

``"grackle_test"``
   Initialize for a grackle 2.0 chemistry and cooling library test
   (depreciated).

``"implosion_2d"``
   Initialize an "implosion" test problem.

``"inclined_wave``
   Initialize an inclined wave test problem. (Primarily used for
   testing the VL+CT MHD solver).
  
``"pm"``
   Initialize ``"dark"`` matter particles in either a regular uniform
   array with one particle per cell, or randomly following the ``"density"``
   field distribution.

``"ppml_test"``
   Initialize fields for the PPML solver for a high-pressure shpere in
   an anisotropic magnetic field.

``"sedov"``
   Calls either ``"sedov_array_2d"`` or ``"sedov_array_3d"``
   initializer depending on the problem rank.

``"sedov_array_2d"``

   Initialize a regular 2D array of Sedov blast problems.  Used for
   parallel-scaling studies without load balancing.

``"sedov_array_3d"``

   Initialize a regular 3D array of Sedov blast problems.  Used for
   parallel-scaling studies without load balancing.

``"sedov_random"`` [ Thomas Bolden ]

   Initialize a regular 3D array of Sedov blast problems.  Used for
   parallel-scaling studies with load balancing.

``"shock_tube"``
   Initialize an axis-aligned shock tube test problem (Primarily used for
   testing the VL+CT MHD solver).

``"soup"``
   
   Similar to the ``"sedov"`` problem, but with letters instead of spheres.

``"trace"``

   Initialize ``"trace"`` tracer particles in either a regular uniform
   array with one particle per cell, or randomly following the
   ``"density"`` field distribution.

``"turbulence"``
   
   Initialize fields for driving turbulence, including ``"driving_[xyz]"``
   fields.

``"value"``

   Initialize fields using expressions directly from the parameter
   file.

   .. warning::

      For technical reasons, ``"value"`` does not work reliably for
      multi-node problems.

``"vlct_bfield"``

   Initialize the cell-centered magnetic fields for use by the VL + CT
   method.  This initialization can be performed from expressions
   specified in the parameter file for each component of the magnetic
   vector potential (in this case, the face-centered magnetic fields
   are also initialized) OR from the face-centered magnetic fields
   that have already been initialized by a separate
   initializer. Additionally, it provides the ability to update
   partially initialized ``"total_energy"`` fields with the specific
   magnetic energy computed from the newly computed cell-centered
   bfields and pre-initialized ``"density"`` fields.
