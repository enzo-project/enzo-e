.. include:: ../roles.incl

**********************
Inference Array Design
**********************

.. toctree::

============
Requirements
============


========
Approach
========

======
Design
======

.. image:: level-array.png
           :width: 800

The purpose of the EnzoMethodInferenceArray method is to create arrays
of block field data to pass to an external machine learning inference
method, invoke the inference method, then update associated block data
based on the inference method results.

Some characteristics include:

   1. arrays are about 64^3
   2. arrays are associated with a specific non-negative refinement level
   3. arrays are only created where needed

Data must be copied from block fields to arrays, using either linear
restriction or prolongation, and back again. Inference arrays are
created based on some simple criteria, such as density threshold,
possibly coupled with a restriction on the block's minimum refinement
level.

Some assumptions we make include:

   1. inference array edges are aligned with coarse-grid cells (but
      not necessarily coarse-level blocks)
   2. inference arrays may overlap each other if "ghost zones" are used
   3. A given block may overlap multiple inference arrays
   4. A given inference array may overlap multiple blocks

Since inference arrays are associated with a specific refinement
level, we will henceforth use the term "level array" to refer to the
array. The level array can be implemented as a (sparse) 3D Charm++
chare array.

Some issues include the following:

   1. multiple level array "create" requests may be received from overlapping blocks
   2. level arrays should be distributed across the machine to prevent
      memory overflow
   3. level arrays won't a priori know in which levels overlapping
      blocks live
   4. synchronization between blocks and inference arrays will be needed
   5. level arrays will tend to be clustered, which can cause
      load imbalances

Phases of the algorithm include the following:

   1. Blocks apply criteria to determine which level arrays to create
   2. Create level arrays
   3. Level arrays request data from overlapping blocks
   4. Blocks pack and send data to level arrays
   5. Level arrays receive and unpack block data
   6. Level arrays call ML inference model
   7. Level arrays send results to overlapping blocks
   8. Leaf Blocks process recevied data

These phases are described in more detail below:

Phase 1. Blocks apply criteria to determine which level arrays to create
========================================================================

Control enters the Method at the Block level, in which all blocks call
``Method::compute()``. Blocks in negative levels will not be involved
in the Method, so can immediately call ``compute_done()``.

Initially we use the "Refine" capability used for mesh adaptation as
the criteria whether to create a given inference array. This criteria
may be combined with additional criteria, in particular a minimum
refinement level for the block.  Note that this does `not` mean all
blocks overlapping a given level array will satisfy that minimum
refinement level, since not all overlapping blocks need to satisfy the
"create" criteria.

After a leaf block applies the criteria, if any cells satisfy the
criteria, the associated overlapped level arrays are tagged for
creation. While the easiest way to do this would be to send a message
directly to the level array element, that element will by definition
not exist until it is created. We instead send the request to the
ancestor root (level-0) block, which uses an array mask to keep track
of which overlapping level arrays need to be created.

Synchronization is required to determine when this phase is complete.
Since root blocks span the domain, we require all leaf blocks to
update their root block ancestor (that is, all leaf blocks send their
"create inference array" result, whether false or true). The root
block will update a volume-based counter with each receive, and all
root blocks will in turn synchronize with the root process.

.. code-block:: C++

   EnzoMethodInference::compute(Block * block)
   {
     // Negative level blocks can exit immediately
     if (block->level() < 0) block->compute_done();

     // Leaf blocks apply criteria and forward to base level
     if (block->is_leaf()) {
       bool create_array = apply_criteria_(block);
       // Result is generated on base level blocks
       send_create_array_(block,create_array);
     }
   }

   void send_create_array_(Block * block, bool create_array)
   {
     if (block->level() > level_base_) {
       // If level > base level, forward to parent
       block_proxy[index_parent]->p_method_inference_send_create_array(create_array);
     } else if (block->level() == level_base_) {
       // If level == base level, notify level arrays
       level_array_mask[*] = false;
       for (cell in cell_create_list) {
           // first determine level arrays to create using mask...
           for (inference array overlapping cell) {
             level_array_mask[infrence_array] = true;
           }
       }
       for (level_array) {
         // ...then create level arrays
         level_array_proxy[level_array] = CkNew();
       }
     }
   }

Phase 2. Create level arrays
============================

Control after phase 1 is in the root-level Simulation object, which
needs to call ``doneInserting()`` to finalize the insertion of the
level array chare elements. We may need some extra synchronization to
ensure all elements are actually created, since Block B calling
``LevelArray::ckNew()``, then calling
``EnzoSimulation::p_done_creating()``, does not ensure that
``ckNew()`` will actually be called before ``p_done_creating()`` is
called.  Synchronization could be handled by the LevelArray
constructor.

After ``doneInserting()`` is called, we call the LevelArray elements
to request data from the blocks in the next phase.

Phase 3. Level arrays request data from overlapping blocks
==========================================================

Control begins in ``LevelArray::p_request_data()``. A level array computes
its overlapping root blocks, and calls
``Block::p_request_level_array_data()`` to request block data. The root
blocks will forward requests to any overlapping child blocks. Note not
all child blocks may overlap the level array.

Alternatively, we could have all leaf blocks in levels L >= K restrict
their data to their level-K ancestors. This could be an effective way
to "prefetch" much of the data required to fill the inference arrays,
leaving only interpolating from overlapping coarse blocks, if any.

Phase 4. Blocks pack and send data to level arrays
==================================================

Control begins in ``Block::p_request_level_array_data()`` in the leaf
blocks. These will pack all required field data and send it to the
LevelArray. Field data on a Block in level L > K will be restricted
(L-K) times before sending to the LevelArray. Blocks in level L < K
will send full data, which will be prolonged (K-L) times at the
receiving end. Note a given Block may overlap multiple Level Array
elements.

If the "prefetch" approach in the previous phase is used, only blocks
in levels L <= K are involved.

Phase 5. Level arrays receive and unpack block data
===================================================

Level Arrays will receive data from overlapping blocks via
``LevelArray::p_receive_block_data()``.  Synchronization will be via a
volume-based counter or by counting updated overlapping cells.  When
all field arrays are updated, the next phase is invoked.

Phase 6. Level arrays call ML inference
=======================================

After a Level Array has initialized all of its arrays from the
received Block data, it calls the external inference method, after
which the next phase is immediately invoked (no synchronization is
required).

Phase 7. Level arrays send results to overlapping blocks
========================================================

After the level array has been processed, data may need to be sent to
the overlapping Blocks. This is done by using essentially the reverse
process of receiving the data. Data will be sent to overlapping root
blocks, which are forwarded to the leaf blocks.

Synchronization is done by the receiving blocks counting the number of
receives until that count reaches the number of overlapping blocks.

Phase 8. Leaf Blocks process recevied data
==========================================

Lastly, the leaf blocks self-update based on received data.

To finalize the Method, all (non-negative level) Blocks must call
``Block::compute_done()``. This is done starting at the leaf-level,
with each leaf block calling ``compute_done()`` after completing this
phase, then calling ``p_method_level_array_done()`` on its
parent. Non-leaf blocks will count receives until ``num_children()``
are reached, call ``compute_done()``, and call
``p_method_level_array_done()`` on its parent, if any.

