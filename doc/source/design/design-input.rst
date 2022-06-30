.. include:: ../roles.incl

*************************
Checkpoint/Restart Design
*************************

.. toctree::

============
Requirements
============

Three code functional requirements of I/O in Cello are:

  1. writing data dumps for subsequent reading by external
     analysis/visualization applications
  2. writing checkpoint files, and
  3. reading checkpoint files to restart a previously run simulation

(While writing image files such as "png" files is also included in the
I/O component of Cello, here we focus on HDF5 files containing
actual block data.)

Additionally, writing and reading disk files must be scalable to the
largest simulations runnable on the largest HPC platforms available,
which necessarily include the largest parallel file systems available.

This scalable I/O approach has been implemented for
checkpoint/restart, and will be adapted for use with data dumps in the
near future.

-----------------------------
Outline of implementation
-----------------------------

The implementation includes modifications to how blocks are mapped
to files, what data are written to the files, and how file I/O is
parallelized.

Ordering
--------

The approach involves a generalization of the existing
``MethodOutput`` method, but enables load-balancing of data between
disk files through the use of block `orderings` to define how blocks
are mapped to files. Currently, the ordering used in ``MethodOutput``,
which is implicit and embedded in the code, is based on a regular
partitioning of root-level blocks together with their descendents. The
updated implementation factors out this ordering into an ``Ordering``
class, provides a Morton space-filling curve ordering, and allows
enables defining other orderings, such as Hilbert curves

File content
------------

The content of the data files must be augmented to include all state
data required to recreate a previously saved AMR block array on
restart. Some information such as block connectivity are generated as
blocks are inserted into the mesh hierarchy. Other information such as
method or solver parameters are not stored, but are taken from
the parameter file. This allows for "tweaking" of parameters on
restart, for example to adjust refinement criteria or solver
convergence criteria.

Control flow
------------

Control flow is handled by separate ``IoWriter`` or ``IoReader`` chare
arrays, where each element is associated with a single HDF5
file. Advantages over previous approaches are better load-balancing of
I/O operations, and decoupling of I/O operations from the Block chare
array. For Enzo-E checkpoint/restart data in particular,
``IoEnzoReader`` and ``IoEnzoWriter`` chare arrays are used.

======
Design
======


Components of the new I/O approach include

  1. Control management

     * ``control_restart.cpp``

        - ``Main::r_restart_enter()``
        - ``Main::p_restart_done()``
        - ``Main::restart_exit()``

  2. New Classes

     * ``EnzoMethodCheck``
     * ``IoEnzoReader``
        - ``IoEnzoReader::IoEnzoReader()``
     * ``IoEnzoWriter``
        - ``IoEnzoWriter::IoEnzoWriter()``
     * ``IoReader``
        - ``IoReader::IoReader()``
     * ``IoWriter``
        - ``IoWriter::IoWriter()``
     * ``MethodOrderMorton``

----------
Algorithms
----------

Output: checkpoint
------------------

.. image:: io-output.png
           :width: 800

Input: restart
--------------

The UML sequence diagram below shows how the Simulation group,
IoReader chare array, and Block chare array interoperate to read data
from a checkpoint directory. Time runs vertically starting from the top,
and the three Charm++ group/arrays are arranged into three columns.

.. image:: io-read.png

startup
-------

In the "startup" section.

Level 0
-------

In the level-0 or root-level section,

Level k
-------

In the level-k section,

cleanup
-------

In the cleanup section,

-------
Classes
-------

EnzoMethodInput

===========
Data format
===========

Data for a given checkpoint dump are stored in a single checkpoint
directory, specified in the user's parameter file using the
``Method:check:dir`` parameter.

The number of data files in the directory is specified using the
``Method:check:num_files`` parameter. A rule-of-thumb is to use the
same number of files as (physical) nodes in the simulation.

Data files are named ``block_data-`` `x` ``.h5``, where 0 <= x <
``num_files``. The format of data files is given in the next section.


Each data file has an associated `block-list` text file named
``block_data-`` `x` ``.block_list``. The block-list file contains a
list of all block names in the associated data file, together with each
block's mesh refinement level. There is one block listed per line, and
the block name and level are separated by a space.

A ``check.file_list`` text file is also included, which includes the
number of data files, and a list of the file prefixes ``block_data-`` `x`.

Note all blocks are included in the files, not just leaf-blocks, and
including blocks in "negative" refinement levels.

------------------
Data file contents
------------------

The HDF5 data files are used to store all block state data, as well as
some global data.

Simulation attributes
---------------------

Metadata for the simulation are stored in the top-level "/" group.
These include the following:

* `cycle`: Cycle of the simulation dump.
* `dt`: Current global time-step.
* `time`: Current time in code units.
* `rank`: Dimensionality of the problem.
* `lower`: Lower extents of the simulation domain.
* `upper`: Upper extents of the simulation domain.
* `max_level`: Maximum refinement level.

Block attributes
----------------

Block attributes and data are stored in HDF5 groups with the same name
as the block, e.g. "B00:0_00:0_00:0".

Block attribute data include the following:

* `cycle`: Cycle of this block.
* `dt`: Current block time-step.
* `time`: Current time of this block.
* `lower`: Lower extents of the block.
* `upper`: Upper extents of the block.
* `index`: Index of the block, specified using three 32-bit integers.
* `adapt_buffer`: Encoding of the block's neighbor configuration.
* `num_field_data`: currently unused.
* `array`: Indices identifying the octree containing the block in the "array-of-octrees".
* `enzo_CellWidth`: Corresponds to the EnzoBlock ``CellWidth`` parameter.
* `enzo_GridDimension`: Corresponds to the EnzoBlock ``GridDimension`` parameter.
* `enzo_GridEndIndex`: Corresponds to the EnzoBlock ``GridEndIndex`` parameter.
* `enzo_GridLeftEdge`: Corresponds to the EnzoBlock ``GridLeftEdge`` parameter.
* `enzo_GridStartIndex`: Corresponds to the EnzoBlock ``GridStartIndex`` parameter.
* `enzo_dt`: Corresponds to the EnzoBlock ``dt`` parameter.
* `enzo_redshift`: Corresponds to the EnzoBlock ``redshift`` parameter.

Block data
----------

Block data are stored as HDF5 datasets.

Fields are currently stored as
arrays of size ``(mx,my,mz)``, where ``mx``, ``my``, and ``mz`` are
the dimensions of the field data `including` ghost data. (Note that
future checkpoint versions may only include non-ghost data to reduce
disk space.) Dataset names are field names with ``"field_`` prepended,
for example ``"field_density"``.

Particles are stored as one-dimensional HDF5 datasets, one dataset per
attribute per particle type. Datasets are named using ``"particle"`` +
`particle-type` + `particle attribute`, delimited by underscores. For
example, ``"particle_dark_vx"`` for the x-velocity particle attribute
``"vx"`` values of the ``"dark"`` type particles in the block.  The
length of the arrays equals the number of that type of particle in the
block.

