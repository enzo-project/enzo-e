************
Adapt Design
************

In the *adapt phase*, blocks may refine or coarsen to adapt to the
evolving resolution requirements of a simulation.  The main
complication is enforcing the "level-jump" condition, which prohibits
adjacent blocks from being in non-consecutive mesh refinement levels.


Maintaining the level-jump conditions may require refining blocks that
would not otherwise be refined, or may require not coarsening blocks
that would otherwise be coarsened.  The process of refining blocks in
a mesh hierarchy solely to maintain the level-jump condition across
block faces is called *balancing* the mesh (not to be confused with
*dynamic load balancing*)

Figure 1. illustrates the steps used in the adapt phase.  Suppose we
begin with the mesh hierarchy at the left, which contains seven
blocks: three in a coarse level and four in the next finer level.  The
first step involves applying local refinement criteria to each block;
in this particular example, the center-most fine block is tagged for
refinement, here indicated by a "+" in the left-most image.

.. figure:: adapt-balance.png

   **Figure 1.** Refining a block (left) may trigger further
   refinements (center) to maintain the level-jump condition.

If we were to only refine this block, however, level jumps would be
introduced across the faces marked by red lines in the center image
(here we optionally include corners as "faces").  The final mesh after
completing the balancing step is shown on the right.

We note that blocks marked for refinement solely to maintain the
level-jump condition may themselves trigger further refinement in
neighboring blocks.  While cascades repeat multiple steps, blocks in
each successive step are in a coarser level, so cascades are always
guaranteed to terminate.  However, cascades still complicate
parallelizing the algorithm, since a given block may not immediately
know whether it needs to refine (or not coarsen) or not, so
determining when the balancing step of the adapt phase is actually
complete is non-trivial.

===================================
Revised adapt algorithm description
===================================

In this section we describe a revised algorithm for the adapt phase in
Enzo-E/Cello. This algorithm was developed by Phil Miller, and was
originally presented in his University of Illinois at Urbana-Champaign
Ph.D. Dissertation, `Reducing synchronization in distributed parallel
programs <\http://hdl.handle.net/2142/95305>`_ (2016). The previous
parallel algorithm relied on Charm++'s support for **"quiescence
detection"**, which is defined as *"the state in which no processor is
executing an entry point, no messages are awaiting processing, and
there are no messages in-flight"* (see `The Charm++ Parallel
Programming System
<https://charm.readthedocs.io/en/latest/charm++/manual.html#quiescence-detection>`_)
Getting the algorithm to work "correctly" required considerable effort
and debugging, and even after several years of development on Enzo-E /
Cello, users still occasionally ran into issues of level-jumps in the
resulting mesh hierarchy, which catastrophically stops simulations
with an error message.

Miller's algorithm avoids using quiescence detection in favor of a
more direct approach.  First, as with the previous algorithm, each
block evaluates its local adapt criteria to determine whether it needs
to refine, stay in the same level, or can coarsen.  Next, bounds on
mesh levels are estimated for each block and communicated with
neighbors.  Bounds for a block may be adjusted as updated bounds
arrive from neighboring blocks.  When a block's minimum and maximum
levels match, the block's next level is decided, and it calls a global
reduction.  All leaf blocks are guaranteed to reach this state, which
can be easily proven by induction on the mesh level starting with the
finest level.

Before presenting the algorithm, we define the following notation:

* :math:`B_i` *block i*
* :math:`B_j` *a block adjacent to block i*
* :math:`L_i^{k}` *the level of Block i in cycle k*
* :math:`\hat{L}_i^{k+1}` *block i's desired next level as locally-evaluated from refinement criteria*
* :math:`\underline{L}_i^{k+1} \leq L_i^{k+1} \leq \bar{L}_i^{k+1}`: *current level bounds which are dynamically updated*
* :math:`L_i^{k+1}` *the next level which is decided when* :math:`\underline{L}_i^{k+1} = \bar{L}_i^{k+1}`
 
We can now write the two main conditions that we use to initialize and
update the level bounds:

* :math:`|L_i^k - L_i^{k+1}| \le 1` *the (temporal) level-jump condition: a block can refine or coarsen at most once per adapt cycle*
* :math:`|L_i^{k} - L_j^{k}| \le 1` *the  (spacial) level-jump condition: refinement levels of adjacent blocks can differ by at most one*

Level bounds are initialized to be :math:`\underline{L}_i^{k+1}
\leftarrow \hat{L}_i^{k+1}` and :math:`\bar{L}_i^{k+1} \leftarrow
L_i^{k} + 1`. That is, the minimum level is initially the level
determined by the local refinement criteria, and the maximum level is
initially one level of refinement more than the current level.

The balancing step of the algorithm proceeds by alternately sending a
block's level bounds to its neighbors, and updating the block's level
bounds given received updated bounds from its neighbors. Bounds are updated
according to the following:

:math:`\underline{L}_i^{k+1} \leftarrow \max ( \underline{L}_i^{k+1}, \max_j (\underline{L}_j^{k+1} - 1))`

:math:`\bar{L}_i^{k+1} \leftarrow \max ( \underline{L}_i^{k+1}, \max_j(\bar{L}_j^{k+1} - 1))`

The minimum bound, which acts like the current desired level in the
previous algorithm, is updated if any neighbor's minimum bound is
greater than one plus the block's current minimum bound.  The maximum
bound, which is used to determine when the balancing algorithm
terminates, is defined as the maximum of the minimum bound, and the
maximum of all neighboring maximum bounds minus one. Note that in
general the maximum bound can only be defined after all neighboring
blocks have been heard from.

