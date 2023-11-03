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

========
Overview
========

When writing a ``Method`` class, it's useful to understand how it is used by ``Cello/Enzo-E``.

Recall that when you launch ``Enzo-E``, you specify how many PEs (processing elements) should be used on the command line.
During startup, a list of ``Method`` objects are constructed on each PE, based on the parameter file (this list is managed by the PE's ``Problem`` instance).
It's also useful to remember simulation data (e.g. fields and particles) are associated with ``Block`` objects.
Throughout the simulation each PE is responsible for evolving a local set of 1 or more ``Block`` objects (note that load-balancing can theoretically migrate ``Block`` objects).

Before each compute cycle, ``Cello/Enzo-E`` determines the current timestep.
For each local ``Block``, a PE invokes the ``timestep`` for each of its ``Method`` instances to determine constraints on the next timestep.
Some other considerations (e.g. user-specified scheduling of operations/stopping at certain simulation times) may alter the duration of the timestep.
Finally, a reduction is performed to pool together the constraints from all blocks.

During the compute cycle, each PE executes a control flow similar to the following code snippet.
Be mindful, that the following snippet doesn't actually exist anywhere in the codebase.

..  code-block:: c++

  void call_compute_on_all_methods(std::vector<Method*> &method_l,
                                   std::vector<Block*> &local_block_l){
    std::size_t num_methods = method_l.size();
    std::size_t num_local_blocks = local_block_l.size();

    for (std::size_t method_ind = 0; method_ind < num_methods; num_methods++){
      Method* cur_method = method_l[method_ind];

      for (std::size_t j = 0; j < num_local_blocks; j++){
        /* do some fancy stuff related to refreshing fields with data from 
           neighboring blocks */

        // call compute on the method
        cur_method->compute(local_block_l[j])
      }

      /* apply a synchronization barrier to make sure that all blocks across
       * all processes have finished completing the current Method.
       *
       * (essentially, wait for block->compute_done() to be called on every
       * block...)
       */
    }

  }
  
Pitfalls
========
The previous section should have made it clear that a given ``Method`` instance generally has its ``timestep`` and ``compute`` method invoked on one or more ``Block`` per cycle.
Consequently, problems can arise if you mutate the attributes of a ``Method`` instance based on data from a given ``Block`` instance.

A good rule-of-thumb for new developers is that you should generally avoid mutating attributes ``Method`` object outside of the constructor.
If you need to associate data with a given ``Block``, you should consider using one of the specialized data interfaces that exist for:

  * field data (managed by ``Field``)
  * particle data (managed by ``Particle``)
  * scalar data (managed by ``Scalar``)

An advantage of using these interfaces is that the associated data will be appropriately migrated if a ``Block`` migrates between PEs.
In certain cases one might alternatively add an attribute to ``EnzoBlock``, but that's generally discouraged if it can be avoided (the ``Scalar`` interface is usually a better choice).

As an aside, there may be times where it makes sense to violate this guideline (e.g. to facillitate optimizations).

Standard properties tracked in base class
=========================================

All ``Method`` classes provide some standardized properties that are managed through the base class.

In some of the following cases, we will talk about how the parameter gets specified for a hypothetical method called ``"my_method"`` (in this hypothetical scenario, subclass's implementation of :cpp:func:`~Method::name` would return ``"my_method"``).

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
