.. include:: ../roles.incl

************************
Enzo-E Coding Guidelines
************************

**[ This page is under development:** *last updated 2019-03-26* **]**


This page describes current coding guidelines for Enzo-E development.  It is
a working document and always open to suggestions for changes.  

-----------------------

============
Naming files
============

All Enzo-E-related files should be in the ``src/Enzo`` subdirectory of
the Enzo-E source code directory.  File names have the form
*component*\ ``_``\ *class*\ ``.[hc]pp``, where *component* is the
high-level code component (for Enzo-E code it is always "``enzo``"),
*class* is the class-name (see below for Enzo-E class naming
guidelines), and extension ``.hpp`` for header (declaration) files and
``.cpp`` for source (definition) files.

+---------------------+---------------------------------------------+
| File type           |  Enzo-E file name                           |
+=====================+=============================================+
| Physics method      |  ``enzo_EnzoMethod``\ *Name*\ ``.[hc]pp``   |
+---------------------+---------------------------------------------+
| Linear solver       |  ``enzo_EnzoSolver``\ *Name*\ ``.[hc]pp``   |
+---------------------+---------------------------------------------+
| Initial conditions  | ``enzo_EnzoInitial``\ *Name*\ ``.[hc]pp``   |
+---------------------+---------------------------------------------+
| Boundary conditions | ``enzo_EnzoBoundary``\ *Name*\ ``.[hc]pp``  |
+---------------------+---------------------------------------------+
| Interpolation       | ``enzo_EnzoProlong``\ *Name*\ ``.[hc]pp``   |
+---------------------+---------------------------------------------+
| Restriction         | ``enzo_EnzoRestrict``\ *Name*\ ``.[hc]pp``  |
+---------------------+---------------------------------------------+

==============
Naming classes
==============

Enzo-E classes all begin with ``Enzo``, followed by the (capitalized)
general type indicating the type of Cello class ("``Method``",
"``Solver``", "``Initial``" for initial conditions, etc.), followed by
the specific name for the class type.  For example, the Conjugate
Gradient (CG) linear solver class is named ``EnzoSolverCg``.  Note
that all code implementing an Enzo-E class should (generally) live in
the files ``enzo_``\ *class-name*\ ``.hpp`` and ``enzo_``\ *class-name*\ ``.cpp``.

+---------------------+---------------------------+------------------+
| Class type          |  Enzo-E class name        | Cello base class |
+=====================+===========================+==================+
| Physics method      |  ``EnzoMethod``\ *Name*   | ``Method``       |
+---------------------+---------------------------+------------------+
| Linear solver       |  ``EnzoSolver``\ *Name*   | ``solver``       |
+---------------------+---------------------------+------------------+
| Initial conditions  | ``EnzoInitial``\ *Name*   | ``Initial``      |
+---------------------+---------------------------+------------------+
| Boundary conditions | ``EnzoBoundary``\ *Name*  | ``Boundary``     |
+---------------------+---------------------------+------------------+
| Interpolation       | ``EnzoProlong``\ *Name*   | ``Prolong``      | 
+---------------------+---------------------------+------------------+
| Restriction         | ``EnzoRestrict``\ *Name*  | ``Restrict``     |
+---------------------+---------------------------+------------------+

====================
Naming class methods
====================

Methods (functions associated with a specific class) are generally
named beginning with a lower-case letter, and underscores for spacing.
Public methods end in an underscore ``_``.

+---------------------------------+-------------------------+
| Method type                     | Method name             |
+=================================+=========================+
| public methods                  | ``a_public_thing()``    |
+---------------------------------+-------------------------+
| private methods                 | ``a_private_thing_()``  |
+---------------------------------+-------------------------+
| Charm++ entry methods           | ``p_blah()``            |
+---------------------------------+-------------------------+
| Charm++ reduction entry methods | ``r_reduce()``          |
+---------------------------------+-------------------------+

Note that Charm entry methods have very different behavior than
regular C++ methods--they are (usually) asynchronous, return
immediately, and are called using a Charm++ proxy to a class typically
residing on a different processing element in a different memory
space.  So it's important to name entry methods in a way that they are
obviously entry methods!  See the Charm++ manual for more details.


====================
Accessing Field data
====================

The (current) approach to accessing Cello Field block arrays is shown
below.  The table summarizes the suggested names for Field-related
variables used to access the array elements, and the following code
initializes active zones (non-ghost zones) of the ``"density"`` field
to be equal to 0.0.

A (preferred) approach would be to use a multi-dimensional array
template, such as ``MultiArray`` in the boost library.  Direct support
for ``MultiArray`` will at some point be directly supported in Cello,
though for the time-being the application developer can wrap the raw
array pointers provided by Cello using the multidimensional array
template of their choice.  Using array templates will improve code
clarity and safety, though may result in some performance degredation.

Summary of field attributes
---------------------------

+------------------------+-----------------+
| Field-related variable | suggested names |
+========================+=================+
| Array dimensions       | ``mx,my,mz``    |
+------------------------+-----------------+
| Active region size     | ``nx,ny,nz``    |
+------------------------+-----------------+
| Ghost zone depth       | ``gx,gy,gz``    |
+------------------------+-----------------+
| Loop variables         | ``ix,iy,iz``    |
+------------------------+-----------------+

Sample code for clearing active density zones
---------------------------------------------

|    ``Field field = cello::field();``
|    ``id = field.field_id("density");``
|    ``field.dimensions (id, &mx, &my, &mz);``
|    ``field.ghost_depth(id, &gx, &gy, &gz);``
|
|    ``enzo_float * d = (enzo_float *) field.values(id);``
|
|    ``for (int iz=gz; iz<mz-gz; iz++) {``
|       ``for (int iy=gy; iy<my-gy; iy++) {``
|          ``for (int ix=gx; ix<my-gx; ix++) {``
|             ``int i = ix + mx*(iy + my*iz);``
|             ``d[i] = 0.0;``
|          ``}``
|       ``}``
|    ``}``

 
===================
Accessing Particles
===================
