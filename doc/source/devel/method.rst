***************
Writing Methods
***************

*[ This page is under development ]*

All method classes are subclasses of the ``Method`` class. The virtual
methods that concrete classes override are shown below:

.. doxygenfunction:: Method::compute

.. doxygenfunction:: Method::name

.. doxygenfunction:: Method::timestep

.. doxygenfunction:: Method::compute_resume


.. note::
   This page is very incomplete. Among other things, we have not discussed
   ``pup`` routines.


Standard properties tracked in base class
=========================================

All ``Method`` classes provide some standardized properties that are managed through the base class.

In some of the following cases, we will talk about how the parameter get's specified for a hypothetical method called ``"my_method"`` (in this hypothetical scenario, subclass's implementation of :cpp:func:`~Method::name` would return ``"my_method"``).

Courant Number
~~~~~~~~~~~~~~

All ``Method`` classes have an associated courant condition.
For a method named ``"my_method"``, the courant value is specified via the parameter called ``Method:my_method:courant``.


This parameter is automatically parsed by machinery in the ``Cello`` layer and the machinery will update the ``Method`` objects with this value right after the constructor is called.
The value of this parameter can be accessed in a ``Method`` subclass with the following function:

.. doxygenfunction:: Method::courant

This parameter is usually accessed in the subclass's implementation of :cpp:func:`~Method::timestep`.

At this time, developers should avoid parsing and tracking the courant value separately within the subclass.

.. note::
   This should not be confused with the ``Method:courant`` parameter.
   This parameter specifies a global courant factor that should never be touched by a ``Method`` subclass.
   Instead, this parameter is entirely handled by the rest of the ``Cello`` infrastructure.

Scheduling
~~~~~~~~~~

All ``Method`` classes support the ability to be scheduled.
For a method named ``"my_method"``, the schedule is specified via a subgroup called ``schedule``.
The rules for specifying a schedule are fairly standard and are described elsewhere in the documentation.

The initialization and usage of an associated schedule are all handled by external ``Cello`` machinery.
A ``Method`` subclass should never need to interact with it (in fact, interacting with it improperly could cause problems).

Refresh Machinery
~~~~~~~~~~~~~~~~~

``Cello`` provides some standardized machinery for specifying requirements related to the fields and particles that need to be refreshed.
The refresh operations are automatically handled by the Cello machinery prior to calls to the ``method`` class.

Configuration of this machinery is typically handled in the constructor of a ``Method`` subclass.

*[ This section is incomplete ]*

