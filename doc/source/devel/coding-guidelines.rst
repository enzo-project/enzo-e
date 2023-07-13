.. include:: ../roles.incl

************************
Enzo-E Coding Guidelines
************************

**[ This page is under development:** *last updated 2023-07-13* **]**


This page describes current coding guidelines for Enzo-E development.  It is
a working document and always open to suggestions for changes.  

-----------------------

============
Naming files
============

All Enzo-E-related files should be placed in an appropriate
subdirectory of ``src/Enzo`` (the Enzo-E source code directory).
File paths have the form
``src/Enzo/``\ *component*\ ``/`` *class*\ ``.[hc]pp``, where

* *component* is the organizing directory from the Enzo-Layer. You can
  investigate 
* *class* is the class-name (see below for Enzo-E class naming
  guidelines)
* and the ``.hpp`` extension is for header (declaration) files while
  ``.cpp`` is for source (definition) files.

Some examples are shown below for classes that nicely fall into the
Cello class-hierarchy:

.. list-table:: Sample Enzo-E file-naming
   :widths: 10 30 18
   :header-rows: 1

   * - File type
     - Enzo-E file name
     - Other comments :superscript:`1`
   * - Physics method
     - ``src/Enzo/``\ *component*\ ``/EnzoMethod``\ *Name*\ ``.[hc]pp``
     - these are split out over a number of different components
   * - Linear solver
     - ``src/Enzo/``\ *component*\ ``/EnzoSolver``\ *Name*\ ``.[hc]pp``
     - currently, all solvers are located in the ``gravity/solvers`` (since
       they are currently just used by the Gravity Solver)
   * - Initial conditions
     - ``src/Enzo/``\ *component*\ ``/EnzoInitial``\ *Name*\ ``.[hc]pp``
     - currently, these are mostly located in ``initial`` or ``tests``
       (the latter is for initializers that are just used by
       test-problems)
   * - Boundary conditions
     - ``src/Enzo/``\ *component*\ ``/EnzoBoundary``\ *Name*\ ``.[hc]pp``
     - currently, there is one instance of this class within ``enzo-core``
   * - Interpolation
     - ``src/Enzo/``\ *component*\ ``/EnzoProlong``\ *Name*\ ``.[hc]pp``
     - currently, these are all located inside of ``mesh``
   * - Restriction
     - ``src/Enzo/``\ *component*\ ``/EnzoRestrict``\ *Name*\ ``.[hc]pp``
     - currently, these are all located inside of ``mesh``

:superscript:`1` The comments in this table may become somewhat outdated
over time. Your best bet is to always look at the actual file structure in
the source code

All Cello-related files should be in the ``src/Cello`` subdirectory
(it is organized has a flat directory structure).  File names should
have the form *component*\ ``_``\ *class*\ ``.[hc]pp``, where
*component* is the high-level code component (e.g. "``problem``",
"``io``", "``mesh``"), *class* is the class-name, and extension
``.hpp`` for header (declaration) files and ``.cpp`` for source
(definition) files.


==============
Naming classes
==============

Enzo-E classes all begin with ``Enzo``.

Subclasses of the Cello-Hierarcy
--------------------------------

For functionality to be executed within a simulation, that functionality must be executed by a subclass of a Cello class. For that reason, most classes in the Enzo-layer are subclasses of a Cello class.
Such subclasses all begin with ``Enzo``, followed by the (capitalized)
general type indicating the type of Cello class ("``Method``",
"``Solver``", "``Initial``" for initial conditions, etc.), followed by
the specific name for the class type.  For example, the Conjugate
Gradient (CG) linear solver class is named ``EnzoSolverCg``.  Note
that all code implementing an Enzo-E class should (generally) live in
the \ *class-name*\ ``.hpp`` and \ *class-name*\ ``.cpp`` files (more
discussion about file naming can be found in the previous section)

+---------------------+---------------------------+------------------+
| Class type          |  Enzo-E class name        | Cello base class |
+=====================+===========================+==================+
| Physics method      |  ``EnzoMethod``\ *Name*   | ``Method``       |
+---------------------+---------------------------+------------------+
| Linear solver       |  ``EnzoSolver``\ *Name*   | ``Solver``       |
+---------------------+---------------------------+------------------+
| Initial conditions  | ``EnzoInitial``\ *Name*   | ``Initial``      |
+---------------------+---------------------------+------------------+
| Boundary conditions | ``EnzoBoundary``\ *Name*  | ``Boundary``     |
+---------------------+---------------------------+------------------+
| Interpolation       | ``EnzoProlong``\ *Name*   | ``Prolong``      | 
+---------------------+---------------------------+------------------+
| Restriction         | ``EnzoRestrict``\ *Name*  | ``Restrict``     |
+---------------------+---------------------------+------------------+

Classes outside of the Cello-Hierarcy
-------------------------------------

There's a temptation to try to squeeze everything into a subclass
extending the Cello-class hierarchy. For that reason, it's worth
emphasizing that developers should feel free to implement classes that
are not part of this hierarchy. In these cases, the convention is
generally for the class name to start with ``Enzo`` and for the rest
of the name to avoid causing confusion with names in the above table.

This is especially useful to consider when implementing complex
functionality that involve a lot of code. There are certain cases
where this can this can lead to massive classes (the implementation is
well over 1000 lines). In such cases, it can be hard to keep track of
everything that is going on (especially if you are not the original
developer of that code). In these cases, there are sometimes
opportunities where you can define a separate helper class that does
a good job abstracting some subset of this functionality.

The :cpp:class:`!EnzoMethodM1Closure` class provides examples of where
this is done. This class requires a data table to be read in from an
external file. Rather than implementing the functionality to read the
table and access the table data (after it has been read into memory)
within the :cpp:class:`!EnzoMethodM1Closure` class, this functionality
is part of the :cpp:class:`!M1Tables` class. Other examples
arise in :cpp:class:`!EnzoMethodMHDVlct` and
:cpp:class:`!EnzoInitialInclinedWave`.

.. note::

    Developers are encouraged to put classes into their own files. (Of
    course, exceptions can be made - especially if the class is really
    tiny).

    However, there are a handful of examples in the codebase where
    this is not the case. Sometimes the header and source file
    dedicated to a class inheriting from a Cello-class will contain
    declarations/definitions of other classes (that are used to help
    implement the primary class).  Such examples arose for historical
    reasons: the enzo-layer was previously organized with a flat
    directory structure, which made it hard to identify groups of
    files that implemented functionality that was used together.

    In some of these cases, the names of the secondary classes are not
    prefixed with ``Enzo``. This should be avoided going forward
    (especially if the class is defined/implemented in its own
    file/header pair).

    The primary reason for the ``Enzo`` prefix is to make it easier to
    differentiate between classes from the ``Cello`` and ``Enzo``
    layers at a glance.

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

========================
General Coding Practices
========================

We recommend some of the following coding practices:

* Prefer Composition over Inheritance

   * In C++, inheritance is the primary way to implement runtime
     polymorphism. [#f1]_ It is most useful to implement an interface
     as an abstract base class (e.g. :cpp:class:`Method` or
     :cpp:class:`!Initial`). Then all derived classes can be used
     interchangably (after they've been constructed). This is a
     **great** uses of inheritance (that are used throughout Cello).
     The Google C++ Style Guide calls this "Interface Inheritance".

   * Inheritance can also be used in other cases as a mechansim for
     reusing code.

       * For example, one might want to make slight changes to the
         behavior of a base class by subclassing it.
       * Unfortunately, this practice easily/commonly produces code
         where the control flow is less explicit, which makes the code
         harder to reason about.If the base class has a chain of
         methods, where some can be selectively overwritten, there is
         a lot more to think about at any given time.
       * An alternative approach for code reuse in these circumstances
         involves composition. Essentially, you compose your objects
         out of self-contained (possibly reusable) components. These
         components may be stored as temporary variables in a function
         call or as attributes of a class. This typically produces
         much more explicit code. (Note: Interface inheritance may be
         used to make it possible to switch between different
         components).
       * A simple Google search of "composition over inheritance" will
         produce lots of additional discussion about this topic. The
         wikipedia article about this topic can be found `here
         <https://en.wikipedia.org/wiki/Composition_over_inheritance>`_

* Where possible, use ``const`` to denote that a function doesn't mutate an argument or a member function doesn't mutate the members of a class. The concept of immutability helps make code much easier to think about (especially when you aren't the original author). The `C++ core guidelines <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#con-constants-and-immutability>`_ talk a little about why this can be useful.

* In the arguments of a function or the arguments of a member-function, prefer to pass C++ references instead of regular C pointers. There are a few exceptions to this rule:
  * You are passing in an array of data (as a pointer)
  * You explicitly want to allow the argument to be a ``nullptr``.

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


.. rubric:: Footnotes

.. [#f1] There are other approaches for acheiving runtime polymorphism
         in c++. The main alternative used in the codebase is the idea
         of a tagged-union (this is used inside of the
         :cpp:class:`!Parameters` class and the
         :cpp:class:`ViewCollec` class template. It is also leveraged
         by the :ref:`EOS functionality <EOS-Developer-Guide>`). This
         can provide certain speed advantages over inheritance, but
         the tradeoff is that you need to declare the closed set of
         all allowed types as part of the tagged union. In
         contrast, inheritance lets you freely introduce new types in
         other sections of the code.
