.. include:: ../roles.incl

************************
Enzo-E Coding Guidelines
************************

**[ This page is under development:** *last updated 2022-02-24* **]**


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

There are 2 approaches for accessing Cello Field block arrays.

  1. The traditional approach is to directly manipulate the raw pointers managed by Cello.
     These can be accessed with the ``Field::values`` method.

  2. The preferred approach is to use the ``CelloView`` multidimensional array templates that wrap the raw pointers managed by Cello.
     Using ``CelloView``\s improve code clarity and safety at the cost of very marginal losses in performance (according to benchmarks).
     ``Field::view`` is used to help facillitate this approach.
     For more details about ``CelloView``\s, see :ref:`using-CelloView`

Below, we provide a table that summarizes the suggested names for Field-related variables used to access the array elements.
We additionally provide brief descriptions of relevant ``Field`` methods.
Sample code is also included for both approaches.
In each case, the sample code initializes active zones (non-ghost zones) of the ``"density"`` field to be equal to ``0.0``.

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

Relevant ``Field`` methods
--------------------------

The traditional approach for accessing fields relies on:

.. code-block:: c++

   char * Field::values (int id_field, int index_history=0) throw ();

This method returns the pointer to the data for the corresponding field (specified by ``id_field``). If no field can be found, this returns a ``nullptr``.
The inclusion of ghost zones is determined by whether or not they've been allocated.
While invoking this method as part of a ``Method`` object's implementation, it's fairly safe to assume that the returned array will indeed include the ghost zones.

The ``index_history`` argument is used to optionally specify the generation of field values that you want to load (higher values correspond to older generations of values).
The default value, ``0``, will give you the current generation of values (this is usually what you want).

The preferred approach for accessing fields relies upon:

.. code-block:: c++

   template<class T>
   CelloView<T,3> Field::view (int id_field,
                                ghost_choice choice = ghost_choice::include,
                                int index_history=0) throw ();

This returns a `CelloView` that acts as a view of the specified field.
The meaning of the ``id_field`` and ``index_history`` arguments are unchanged from ``Field::values``.
Unlike ``Field::values``, when an invalid ``id_field`` is specified, the program will abort with an error.

The template parameter ``T`` specifies the expected datatype of the field.
If the field does not have the expected datatype, the program will abort with an explanatory error message.
While implementing a ``Method`` object in the ``Enzo`` layer, this parameter is frequently ``enzo_float``.

By default, the returned ``CelloView`` **always** include ghost zones.
This behavior is controlled by the ``choice`` argument, which can be passed the following values:

  * ``ghost_choice::include``: the returned view always includes ghost zones (the program aborts with an error message if ghost zones aren't allocated). This is the default value.
  * ``ghost_choice::exclude``: the returned view always excludes ghost zones.
  * ``ghost_choice::permit``: the returned view includes ghost zones if they are allocated (replicating the behavior of ``Field::values``).

Overloads are also provided for ``Field::values`` and ``Field::view`` that:
  * provide read-only access to field arrays from a ``const Field`` instance
  * let you replace the first argument with a string holding the field's name

Sample code for clearing active density zones
---------------------------------------------

Traditional Approach
~~~~~~~~~~~~~~~~~~~~

.. code-block:: c++

   Field field = cello::field();
   int id = field.field_id("density");
   int mx, my, mz;
   field.dimensions (id, &mx, &my, &mz);
   int gx, gy, gz;
   field.ghost_depth(id, &gx, &gy, &gz);

   enzo_float * d = (enzo_float *) field.values(id);

   for (int iz=gz; iz<mz-gz; iz++) {
     for (int iy=gy; iy<my-gy; iy++) {
       for (int ix=gx; ix<my-gx; ix++) {
         int i = ix + mx*(iy + my*iz);
         d[i] = 0.0;
       }
     }
   }

Preferred Approach
~~~~~~~~~~~~~~~~~~

A shorter version of the following snippet is also possible where we pass ``ghost_choice::exclude`` as the second argument to ``field.view`` to entirely exclude the ghost zone from the array.

.. code-block:: c++

   Field field = cello::field();
   int id = field.field_id("density");
   int gx, gy, gz;
   field.ghost_depth(id, &gx, &gy, &gz);

   CelloView<enzo_float,3> d = field.view<enzo_float>(id);

   // we can get the array shape directly from the array (a minor design quirk
   // may make the arguments for the shape method seem a little unintuitive)
   int mz = d.shape(0);
   int my = d.shape(1);
   int mx = d.shape(2);

   for (int iz=gz; iz<mz-gz; iz++) {
     for (int iy=gy; iy<my-gy; iy++) {
       for (int ix=gx; ix<my-gx; ix++) {
         d(iz, iy, ix) = 0.0;
       }
     }
   }

More code examples that use this functionallity can be found in the implementation of the ``EnzoInitialCloud`` and ``EnzoInitialBCenter`` classes

===================
Accessing Particles
===================

===================
Preprocessor Macros
===================

In case you add a new preprocessor macro or definition that is being used in a ``.ci`` file,
you have to add it to the ``CHARM_PREPROC_DEFS`` list in the main ``CMakeLists.txt`` so that
it is being used/forwarded when processing the ``.ci`` files.

For example, for Grackle we set
``set(CHARM_PREPROC_DEFS ${CHARM_PREPROC_DEFS} "-DCONFIG_USE_GRACKLE ")``, which appends to the
existing list of defintions alredy stored in the ``CHARM_PREPROC_DEFS`` variable.
