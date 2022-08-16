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

.. warning::

   This page is very incomplete. Among other things, we have not discussed
   ``pup`` routines.

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
  
Gotchas
=======
The previous section should have made it clear that a given ``Method`` instance may have it's ``timestep`` and ``compute`` method to one or more ``Block`` per cycle.
Consequently, problems if you mutate the attributes of ``Method`` object based on data from a ``Block``.

A good rule-of-thumb for new developers is that you should generally avoid mutating attributes ``Method`` object outside of the constructor.
Instead you should consider using one of the specialized data interfaces.
This includes field data (managed by ``Field``), particle data (managed by ``Particle``), or scalar data (managed by ``Scalar``).
The advantage of doing this is that the associated data will be appropriately migrated if a ``Block`` migrates between PEs.
In certain cases you can also add an attribute to ``EnzoBlock``, but that's generally discouraged (the ``Scalar`` interface is usually a better choice).

Naturally, this isn't a hard-and-fast rule.
There are times where you may need to violate it (e.g. to facillitate optimizations).
