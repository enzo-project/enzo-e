***********************
Writing Physics Classes
***********************

*[ This page is under development ]*

This page tries to provide an overview for how to write :cpp:class:`!Physics` classes and documents some quirks about existing classes.

All physics classes are subclasses of the :cpp:class:`!Physics` class.
Unlike other classes in the Enzo-E layer that descend from the Cello class-hierarchy (:cpp:class:`!Method`, :cpp:class:`!Initial`, :cpp:class:`!Prolong`, :cpp:class:`Refine`, etc.), the Cello-layer doesn't really interact much with instances of the :cpp:class:`!Physics` classes - beyond storing them.

In practice, :cpp:class:`!Physics` classes are commonly used to store problem-specific configuration information that needs to be accessed by multiple different :cpp:class:`!Method` classes and/or initializers.
Some functions that make use of this information are also sometimes introduced to these classes.

.. _how-to-store-global-data:

When to write a ``Physics`` class
=================================

Other ways to store data
~~~~~~~~~~~~~~~~~~~~~~~~

To understand when it may be useful to write a :cpp:class:`!Physics` class, it's first useful to discuss some of the ways one could access/store cross-cutting configuration data:

1. store this information in a global variable (DON'T DO THIS)

     - From a coding-style perspective, this is almost always the **wrong** way to store such information.

     - Moreover, Charm++ `does not support global variables <https://charm.readthedocs.io/en/latest/charm%2B%2B/manual.html#read-only-data>`_ unless they are properly declared in the .CI (note: the code will still compile, but it should be avoided all the same).

2. access information stored in :cpp:class:`!EnzoConfig`, which is accessible through :cpp:expr:`enzo::config()` (TRY NOT TO DO THIS)

     - this class already stores a lot of values parsed from parameter files.

     - while this approach can be convenient in simple cases, it generally leads to brittle code that is hard to refactor (challenges can come up if you want to alter the way that different options are stored in the parameter file).

     - This approach has been used a fair amount in the past, but we are actively moving away from it (and discourage this approach).

3. Cross-cutting configuration information is commonly only relevant when you are using a particular :cpp:class:`!Method` class.
   In this case, you can store the configuration information within that :cpp:class:`!Method` class and define public instance-methods on that particular class that accesses the information.

     - This approach is strongly prefered over the preceding approaches 1 and 2.
       Refactoring is generally easier in this approach.
       Moreover, this can be easier than defining a :cpp:class:`!Physics` class.

     - For this approach, it's important to understand how to access an instance of a particular kind of :cpp:class:`!Method` at an arbitrary point in the code (after all :cpp:class:`!Method` classes have been constructed).
       One can use :cpp:expr:`enzo::problem()->method("<name>")` to return a pointer to the instance of the :cpp:class:`!Method` class for which :cpp:func:`Method::name()` returns ``"<name>"``.
       If no such instance can be found, the expression returns a ``nullptr``.
       You then need to cast that pointer to the appropriate :cpp:class:`!Method` subclass before you access the information.

     - At the time of writing, this approach is commonly used to store information encoded within the :cpp:class:`!EnzoMethodGrackle` class.
       To access such information, one could write

       ..  code-block:: c++

           const EnzoMethodGrackle *ptr = static_cast<const EnzoMethodGrackle*>
             (enzo::problem()->method("grackle"));
           if (ptr != nullptr) {
             // maybe do stuff with ptr->try_get_chemistry() ...
           }

       In practice, some convenience functions have been written to help with these sorts of operations like :cpp:expr:`enzo::grackle_method()` or :cpp:expr:`enzo::grackle_chemistry()`

     - In principle, one could do something analogous involving subclasses of :cpp:class:`!Initial`, but that could potentially introduce problems during a simulation restart.

Storing information in a ``Physics`` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It's often most useful to encode configuration-information within a :cpp:class:`!Physics` class when there isn't an obvious single :cpp:class:`!Method` class where it should be stored.

A particular scenario where this is relevant is when separate (somewhat-interchangable) :cpp:class:`!Method` classes implement different algorithms to model the same set of physics.
For example, consider the storage/access of the dual-energy formalism configuration.
Since this is mostly relevant in the context of a hydro-solver it may make sense to store this information in the :cpp:class:`!Method` class that encapsulates a hydro-solver.
However, because Enzo-E has :cpp:class:`!Method` classes that implement different hydro-solvers (that use the dual-energy formalism), we instead encode this information in a :cpp:class:`!Physics` class.

There are also scenarios where some configuration information isn't really associated with any singular :cpp:class:`!Method`.
For example the Equation-Of-State is important to a number of different methods.
Another example includes the Gravitational Constant - this is important in self-gravity solvers and external-potential solvers, which are implemented in different :cpp:class:`!Method` classes.

General Tips
============

The general advice is to implement a :cpp:class:`!Physics` class so that it is immutable (after construction the instance's state doesn't change).

This makes the behavior of :cpp:class:`!Physics` classes much easier to reason about because a single PE (processing element)

  - only has one instance of a given :cpp:class:`!Physics` class
  
  - AND is responsible for evolving one or more instances of :cpp:class:`EnzoBlock`.


Quirky Implementations
======================

:cpp:class:`!EnzoPhysicsCosmology` currently tracks some mutable state (e.g. the current scale-factor, the current rate of expansions, current redshift).
This is just something to be mindful of.

It's worth noting that the initialization of :cpp:class:`!EnzoPhysicsFluidProps` and :cpp:class:`!EnzoPhysicsGravity` are a little quirky.
These objects are ALWAYS initialized, regardless of whether a user specifies the names of these objects in the :par:param:`Physics:list` configuration-file parameter.
This choice was made for the sake of maintaining backwards compatability with older versions of parameter-files that were created before these classes were invented (since they encode some information that was previously stored elsewhere).

.. note::

   This means that all simulations have an instance of :cpp:class:`EnzoPhysicsGravity` (regardless of whether or not gravity is actually in use).
