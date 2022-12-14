.. include:: ../roles.incl

.. _Inference Design:

**********************
Inference Array Design
**********************

.. toctree::


=======
Purpose
=======


.. image:: array.png
           :width: 346
           :alt: Image showing an example of the placement of inference arrays in a small cosmology simulation

.. image:: infer.png
           :width: 346
           :alt: Image showing an example of the bubbles output by the inference method

This page describes the design of the ``EnzoMethodInference`` Method.
The purpose of the EnzoMethodInferenceArray method is to create a collection of
regular arrays ("inference arrays"), each containing a subset of block field
data to pass to an external machine learning inference method. After
the inference method is invoked, the overlapping blocks are then
provided with the pertinent output of the inference method, such as
the locations of bubbles where star formation is expected to occur. A
mock-up of inference arrays and generated bubbles is shown in the
above figures.

Some characteristics of inference arrays include:

   1. inference array sizes are typically about 64^3
   2. all inference arrays have the same resolution
   3. all inference arrays have the same size
   4. inference arrays are only created where needed

Inference array locations are determined by some relatively simple
(local) criteria, such as density threshold, possibly coupled with a
restriction on the block's minimum refinement level. Field data are
then copied from overlapping block fields to the inference arrays,
using either linear restriction or prolongation.

Some assumptions we make include:

   1. A given block may overlap multiple inference arrays
   2. A given inference array may overlap multiple blocks
   3. Inference arrays are aligned with blocks in some specific AMR level, ``level_array``
   4. Inference array resolution matches that of blocks in some (finer) AMR level, ``level_infer``.

Since inference arrays are associated with blocks in a specific
refinement level, we use the term "level array" to refer to the sparse array
of inference-arrays. The level array is implemented as a sparse 3D
Charm++ chare array, where each element of the chare array is a
collection of inference arrays for a specific "zone".

Some issues that should be addressed in the design include the following:

   1. multiple level array "create" requests may be received from overlapping blocks, but can only be created once
   2. inference arrays tend to be clustered, so level array elements should be distributed across compute nodes to reduce compute and memory load imbalances
   3. level arrays won't a priori know in which levels overlapping
      blocks live
   4. synchronization between blocks and inference arrays will be needed

======
Phases
======

Below is a UML sequence diagram illustrating the different phases in
``EnzoMethodInference``.  The left blue columns represent inference
arrays, the red right columns represent blocks in different refinement
levels, and the center yellow column represents the root-node
Simulation object, used for synchronization and counting.

.. image:: level-array.png
           :width: 800

Phases of the algorithm include the following:

   1. **evaluate**: Blocks apply criteria to determine where to create inference arrays
   2. **allocate**: The "level array" chare array of inference arrays is created
   3. **populate**: Inference arrays request and receive field data from overlapping blocks
   4. **apply inference**: Inference arrays apply the external ML inference method
   5. **update blocks**: Inference arrays send results back to overlapping blocks

These phases are described in more detail below:

Phase 1. Evaluate
=================

In the "Evaluate" phase, Blocks apply criteria to determine where to
create inference arrays.  Control enters the Method at the Block
level, in which all blocks call ``Method::compute()``. Blocks in
negative levels will not be involved in the Method, so can immediately
call ``compute_done()`` to exit the method.

The criteria currently implemented is whether the point-wise density
is greater than the block-local average by some specified threshold.
To improve performance, this is applied only on "sufficiently fine"
level blocks, specified by "level_base" (level_base=2 is
typical). Inference arrays are guaranteed not to overlap leaf nodes in
levels coarser than level_base. Conversely, all blocks in level =
level_base that overlap inference arrays are guaranteed to exist.

After a leaf block applies the criteria ``apply_criteria()``, if any
cells satisfy the criteria, the associated overlapped level array
elements are tagged for creation. Note there may be multiple such
elements, based on whether the block is coarser or finer than
"level_array" (the level at which blocks and inference arrays coincide
in size). If there are multiple overlapping inference arrays for a
block, a logical "mask" array is used for keeping track of which
inference arrays to create. If only one inference array overlaps a
block, the mask size is 1.

These masks are merged toward the coarser level_base level, using the
``p_merge_masks()`` entry method called on block parents. At each
step, the child masks are merged in their parent block using
logical-OR. When level_base is reached (level 2 in the figure), each
block in the level_base level will have a mask specifying where each
inference array needs to be created. At this step, the level array elements
are created using ``p_create_level_array()``.

The reduction operation continues with counting the number of
inference arrays created, using ``p_count_arrays()``. This continues down to the
root level blocks, which send the accumulated counts to the root Simulation
object. After all root-level block counts have been received, the
Simulation object will contain the total number of inference arrays to
be created.

Phase 2. Allocate
=================

The count of number of inference arrays to create, determined in the
previous phase, is used to determine when all level array elements
have been created.  (As a technicality, the count is set to one more
than the count to prevent the algorithm from hanging if no level array
elements need to be created, which is possibile). As level array
elements are created, the constructor notifies the root Simulation
object, which decrements the counter.  When zero is reached, all level
array elements are guaranteed to have been created, and the Simulation
object can then finalize the chare array by calling the Charm++
"doneInserting()" method, and proceed to the next phase.

@@@@@@@@@@@@@

Phase 3. Populate
=================

Control begins in ``LevelArray::p_request_data()``. A level array computes
its overlapping root blocks, and calls
``Block::p_request_level_array_data()`` to request block data. The root
blocks will forward requests to any overlapping child blocks. Note not
all child blocks may overlap the level array.

Alternatively, we could have all leaf blocks in levels L >= K restrict
their data to their level-K ancestors. This could be an effective way
to "prefetch" much of the data required to fill the inference arrays,
leaving only interpolating from overlapping coarse blocks, if any.

Phase 4. Apply inference
========================

Control begins in ``Block::p_request_level_array_data()`` in the leaf
blocks. These will pack all required field data and send it to the
LevelArray. Field data on a Block in level L > K will be restricted
(L-K) times before sending to the LevelArray. Blocks in level L < K
will send full data, which will be prolonged (K-L) times at the
receiving end. Note a given Block may overlap multiple Level Array
elements.

If the "prefetch" approach in the previous phase is used, only blocks
in levels L <= K are involved.

Phase 5. Update blocks
======================

Level Arrays will receive data from overlapping blocks via
``LevelArray::p_receive_block_data()``.  Synchronization will be via a
volume-based counter or by counting updated overlapping cells.  When
all field arrays are updated, the next phase is invoked.

After a Level Array has initialized all of its arrays from the
received Block data, it calls the external inference method, after
which the next phase is immediately invoked (no synchronization is
required).

After the level array has been processed, data may need to be sent to
the overlapping Blocks. This is done by using essentially the reverse
process of receiving the data. Data will be sent to overlapping root
blocks, which are forwarded to the leaf blocks.

Synchronization is done by the receiving blocks counting the number of
receives until that count reaches the number of overlapping blocks.

Lastly, the leaf blocks self-update based on received data.

To finalize the Method, all (non-negative level) Blocks must call
``Block::compute_done()``. This is done starting at the leaf-level,
with each leaf block calling ``compute_done()`` after completing this
phase, then calling ``p_method_level_array_done()`` on its
parent. Non-leaf blocks will count receives until ``num_children()``
are reached, call ``compute_done()``, and call
``p_method_level_array_done()`` on its parent, if any.

