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

Subclasses of the Cello-Hierarchy
---------------------------------

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

Classes outside of the Cello-Hierarchy
--------------------------------------

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

=============================
Headers and File Organization
=============================

The Cello and Enzo layers of the codebase are both organized into subcomponents, but there are slightly different guidelines dictating the file organization.

Cello
-----

The files are all organized in a flat directory structure. Every file
is prefixed by the name of the component that it belongs to. Individual
headers are not intended to be used. Instead, each component defines 2
standard header files.

1. ``_``\ *component*\ ``.hpp``: This is a private header that is used
   to aggregate the headers for all functions/classes defined as part
   of the component.

2. *component*\ ``.hpp``: This is the public header. It contains
   includes the ``_``\ *component*\ ``.hpp``. It also contains all of
   the necessary include statements for the headers from other
   components that are necessary for the declarations/definitions
   present in the current headers (usually it includes the
   ``_``-prefixed header).

Because the directory structure is flat all include-directives just
specify filenames (there are no paths).

Enzo
----

Over the last several years, the Enzo layer has grown significantly (it
is now comparable in size to the Cello layer), and it is likely to continue
growing. Due to the large (and increasing) size of the Enzo layer, we
take some steps to improve build-times.  We are in the midst of
rolling-out an updated policy.


Traditional Approach
~~~~~~~~~~~~~~~~~~~~

Historically, the Enzo layer was structured just like one of the
components in the Cello Layer in a flat structure. All individual
header files were aggregated inside of the private
``src/Enzo/_enzo.hpp`` header and there was a public header called
``src/Enzo/enzo.hpp``. The files were also originally organized in a
flat directory structure, but that is no longer the case.
Additionally, just about every source file started with:

.. code-block:: c++

   #include "cello.hpp"
   #include "enzo.hpp"

While this approach has been **very** successful and there’s nothing wrong with it *per se*, it does trigger a full rebuild of the entire Enzo layer any time any header file changes.

New Approach
~~~~~~~~~~~~

Under our new approach, the contents of the Enzo layer are now organized into subdirectories, which corresponds to a subcomponents.
Each subcomponent has an aggregate header named after the component (e.g. the ``Enzo/mesh`` subcomponent should have an associated header called ``Enzo/mesh/mesh.hpp``).

.. note::

   Nested subdirectories are handled on a case-by-case basis.  In some
   cases, each nested subdirectory is treated as a separate
   subcomponent. In other cases, nested subdirectories are only used
   for organization-purposes.

**Unlike with the Cello layer,** there is currently no distinction between a private and public header: the aggregate header **is** the public header.
Each public header file should be self-contained.
In other words, the header should compile on its own, without requiring it to be included alongside other header files (or requiring a particular order of include-directives).

To be self-contained, each public header needs to have the necessary include directives so that all of the other aggregated headers within the public header have access to the necessary symbols.

* Within the enzo-layer, all include-directives should specify paths to relative to the ``src`` directory.
  Thus the include directives may look like:

  .. code-block:: c++

     #include "Cello/cello.hpp"
     #include "Cello/mesh.hpp"
     #include "Enzo/enzo.hpp"

  This ensures that there is no ambiguity if there are similarly named subcomponents in the Enzo- and Cello-layers.

* A public header in a given subcomponent should generally only include the public headers from other subcomponents.
  Because including the public header from another component exposes symbols beyond what is strictly necessary, we recommend adding a comment next to the include directive that specifies the name of the symbol (e.g. a class, type, function)that is required from the header.
  By doing this, future developers can more easily remove unnecessary ``#include`` directives if a particular dependency is no longer necessary OR is moved

* It is important to be mindful that include-directives introduce transitive dependencies, and any time a header file changes, all ``.cpp`` files that depend on that header (whether directly or transitively) needs to be recompiled.

* Do your best to ensure that the public header files only include what is necessary and avoid unnecessary include statements (to keep compile-times shorter).
  With that said, it's better to explicitly write out include-directives to other headers, even if that header would transitively be included by some other unrelated header.


.. note::

   In the long-term, it may make sense to make each individual header self-contained, which is the recommended strategy by the `Google C++ Style Guidelines <https://google.github.io/styleguide/cppguide.html#Header_Files>`_

   Ideally, each header would include just the headers defining other symbols (classes/types/functions) that it needs.
   Under this approach, inclusion order would become less problematic and it would further speed up incremental compilation.
   Transitioning to this kind of approach all at once would be intractable since there are currently over 125 header files in the Enzo layer.
   However, this alternative approach is definitely worth exploring and could be implemented on a subcomponent-by-subcomponent basis in the future.

What is actually necssary to include in a header file?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is useful to have a brief discussion about what is actually necessary to include in a header file.

As noted in the `Mozilla Style-Guide <MozillaStyleHeaders_>`_, a full definition of a type is required for the type to be used

* as a Base class,
* as a member/local variable
* in a function declaration, where an instance of the type is passed as an argument, by-value, OR is returned from the function, by value.
* with ``delete`` or ``new``
* as a template argument (in certain cases)
* if it refers to a non-scoped enum type (in all cases)
* when refering to the values of a scoped enum

The `Mozilla Style-Guide <MozillaStyleHeaders_>`_ also notes that a forward declarations of a type will suffice when the type is used

* for declaring a member/local variable that holds a reference or pointer to that type.
* in a function declaration that accepts reference or pointer to that type as an argument or returns a reference or pointer to that type.
* in the definition of a type-alias.

.. note::

   In general, there is a difference of opinion about whether forward declarations should be used.
   For example, the `Mozilla Style-Guide <MozillaStyleHeaders_>`_ and the `LLVM Style Guide <https://llvm.org/docs/CodingStandards.html#include-as-little-as-possible>`_ generally encourage usage of forward declarations.
   In contrast, the `Google Style Guide <https://google.github.io/styleguide/cppguide.html#Forward_Declarations>`_ discourages this practice.

   We don't currently take a strong position on this.
   In the codebase's current state, there are definitely cases where using forward declarations is **VERY** useful (it may be helpful to leave a comment explaining why its useful/necessary in a particular case).

   If we do shift to making **all** header-files self-contained, it may make the codebase more readable if we avoided forward declarations in the future.
   But, we can cross that bridge, when we get to it.

Finally, it's worth noting that there is a large temptation to implement functions or a class's member-functions (aka methods) inside of the header where they are declared.
Sometimes this is unavoidable - it's necessary when defining templates and sometimes it's necessary for performance-purposes (to facilitate inlining of a function called within a tight for-loop).
However, in a lot of cases (especially when implementing a ``virtual`` method), this can and should be avoided.
This practice has a tendency to introduce additional include-directives into a header file that could otherwise just be located within a source file.

For example, a handful of functions in Enzo-E make use of functionality defined in the ``<algorithm>``, ``<random>``, and ``<sstream>`` headers without needing to pass around types defined in these headers between functions.
In these cases, by implementing such functions in ``.cpp`` files, we can directly include these headers in the ``.cpp`` source files and avoid including them in the header files (this can actually save a lot of time during compiling since standard library headers can be large).

As we finish transitioning the Enzo layer to using separate subcomponents, additional opportunities will arise for including headers in source files rather than inside of headers.
For example, most times when you access instances of :cpp:class:`!EnzoPhysicsCosmology`, :cpp:class:`!EnzoPhysicsFluidProps`, or :cpp:class:`!GrackleChemistryData` in the method of a class, ``MyClass``, the header defining ``MyClass``, doesn't actually require knowledge about the definitions of these other classes.
Often times, a method ``MyClass`` will simply use :cpp:expr:`!enzo::cosmology()`, :cpp:expr:`!enzo::fluid_props()`, or :cpp:expr:`!enzo::grackle_chemistry()` to retrieve an instance of one of these classes.
Then the method will query a piece of information stored in the retrieved instance and it will never touch the instance again.


.. _MozillaStyleHeaders: https://firefox-source-docs.mozilla.org/code-quality/coding-style/coding_style_cpp.html#header-files



Questions and Answers about introducing new files to the Enzo Layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**How do I add a new file to the Enzo Layer?**

Usually you will introduce a header (``.hpp``) and source (``.cpp``) file at the same time.

1. Find an appropriate subdirectory to put the new file in.
2. *Identify the relevant aggregate header file.*
   First, check for a header file that shares the name of the current subdirectory (with a ``.hpp`` suffix).
   If that doesn't exist, check for the aggregate header file in the parent directory.
   Repeat the process until you find subdirectory has an aggregate header file (you should only really need to go up 1 level).
3. Once you find the appropriate aggregate header file, add an include statement to your header file (if there aren't include-statements to other header files in the same subdirectory, you're probably in the wrong place).
4. *Identify the relevant* ``CMakeLists.txt`` *file for other files in the current subdirectory.*
   First, check for this file in the subdirectory where you have placed your new files.
   If it doesn't exist, check the parent directory.
   Repeat the process until you find it.
5. Add the paths to your new header and source file to the existing list of other files that are located in the same subdirectory as your new files.
   (If no such list exists, but you see something about ``GLOB``, you may not need to do anything).

**How do I delete a file from the Enzo Layer?**

This is pretty self-explanatory.
Just make sure to delete entries in the aggregate header file and (if applicable) the ``CMakeLists.txt`` file that specify the path(s) to your deleted file.

**How do I move files between subdirectories in the Enzo Layer??**

This is also pretty easy.
Just make sure to delete the old paths from the aggregate header and the ``CMakeLists.txt`` files where they were originally listed.
And, make sure to add the new paths to the appropriate aggregate header and the ``CMakeLists.txt`` files.

Enzo Header Guards
------------------

To avoid issues with a header being included multiple times (i.e. if it is a transitive dependency of another header), we make use of ``#define`` header guards in every header. In general, the header guard symbol looks something like ``ENZO_<PATH>_<FILE>_HPP``, where ``<PATH>`` is replaced by the path to the file and ``<FILE>`` is the filename (excluding the ``.hpp`` suffix).

.. note::

   A number of headers have guard symbols that are unchanged from before the Enzo subdirectory was reorganized.
   As an example, ``EnzoMethodPpm`` was previously defined in ``src/Enzo/enzo_EnzoMethodPpm.hpp`` and had a header guard symbol called ``ENZO_ENZO_METHOD_PPM_HPP`` (the underscore was used as a separator between Camel-Case separated words and in place of the period``.

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
         harder to reason about. If the base class has a chain of
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

* When declaring enumerations, prefer to use a `scoped enum <https://en.cppreference.com/w/cpp/language/enum>`_.
  The concept of a scoped enums was introduced in C++11. 
  An example of a scoped enum from the codebase is the :cpp:enum:`!ghost_choice` enumeration.
  The declaration of this type looks something like the following snippet:

  .. code-block:: c++

     enum class ghost_choice {exclude, include, permit};

  In the above snippet, if the ``class`` keyword is replaced by the ``struct`` keyword, the result is completely equivalent.
  If neither the ``class`` nor ``struct`` keyword were present, then ``ghost_choice`` would be an unscoped enum (or a C-style enum).

  There are three main differences to be aware of when using a scoped enum:

    1. The enumerators must be specified as :cpp:enumerator:`!ghost_choice::exclude`, :cpp:enumerator:`!ghost_choice::include`, and :cpp:enumerator:`!ghost_choice::permit`.
       This is much more explicit (and consequently better) than the alternative.
       If :cpp:enum:`!ghost_choice` were instead defined as an unscoped enum, :cpp:enumerator:`!exclude`, :cpp:enumerator:`!include`, and :cpp:enumerator:`!permit` would be introduced as symbols in the global namespace. [#f2]_

    2. Because :cpp:enum:`!ghost_choice` is a scoped enum, ``ghost_choice`` is directly recognized as a type of a variable. If it were an unscoped enum, you must 
       If it were an unscoped enum, you could instead declare a variable of type ``enum ghost_choice``.

    3. Perhaps most importantly, scoped enums have better type safety.
       Specifically, integer values cannot be implicitly converted to an enum.
       Any conversions from an integer value requires an explicit cast.
       This generally leads to more explicit code.
       Furthermore, it lets the compiler identify cases where the order of an integer argument and a scoped enum argument are accidentally permuted.

* **ERROR REPORTING:** We don't use C++ exceptions in the codebase. In general, when an error arises, we generally abort the program with an informative error message. You can signal that an error occured by using the ``ERROR`` family of macros or you can use the ``ASSERT`` family of macros to have the program conditionally abort if some condition is not satisfied.

  * When you are implementing new functionality, you are encouraged to liberally use the :c:macro:`!ERROR` and :c:macro:`!ASSERT` families of macros to ensure that Enzo-E loudly fails and aborts when the functionality is used in unexpected ways.
    When running a really expensive simulation, a user should generally prefer that a simulation loudly fails.
    The extreme alternative case is for the simulation to run to completion while silently having problems, which likely invalidates the results (and these problems may not be detected until MUCH later).
    Furthermore, it's easy enough for a user to comment out an error message that they wish to ignore.

  * We briefly describe the arguments of :c:macro:`!ERROR`, :c:macro:`!ERROR1`, :c:macro:`!ERROR2`, ..., :c:macro:`!ERROR8`, :c:macro:`!ASSERT`, :c:macro:`!ASSERT1`, :c:macro:`!ASSERT2`, ..., :c:macro:`!ASSERT8` down below:

    * The first argument is always the name of the function where the macro is being invoked (to assist debugging in the future).

    * The second argument is a c-string that provides the error message.
      Printf formatting specifiers can be used within the error message, but the number of specifiers must match the integer at the end of the macro (e.g. :c:macro:`!ERROR2` or :c:macro:`!ASSERT2` expects 2 printf specifiers while :c:macro:`!ERROR` or :c:macro:`!ASSERT2` expects none).

    * The next arguments specifiy the variables used by the formatting specifier (if there are any).

    * The :c:macro:`!ASSERT` macro-family expects 1 last argument: the boolean condition dictating whether the program aborts.

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

.. [#f2] You can mimic the scoping behavior to some degree with an unscoped enum if you place the definition of the unscoped inside of a class or struct.
