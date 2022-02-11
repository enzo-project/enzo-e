.. include:: ../roles.incl

******************
Scalable IO Design
******************
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
I/O component of Cello, here we will focus on HDF5 files containing
actual block data.)

Additionally, writing and reading disk files must be scalable to the
largest simulations runnable on the largest HPC platforms available,
which necessarily include the largest parallel file systems available.

----------------------
Current implementation
----------------------

In Cello, "scalable I/O" is currently implemented using
``MethodOutput``. This is close to being suitable for requirement 1,
but has limitations in scalability and flexibility; in particular,
disk files are not load balanced. Checkpoint and restart (requirements
2 and 3) are currently implemented using Charm++, which has some
limitations in scalability (the file count can be excessive), and in
flexibility (simulations cannot be restarted with modified
parameters).

-----------------------------
Outline of new implementation
-----------------------------

Ordering
--------

Our new approach will involve a generalization of the ``MethodOutput``
method, but will enable load-balanced disk file through the use of
block `orderings`. Currently, the ordering used in ``MethodOutput``,
which is implicit and embedded in the code, is based on a regular
partitioning of root-level blocks together with their descendents. Our
new approach will factor out this ordering, and allow flexibility in
defining other orderings, such as space-filling curves like Hilbert or
Morton.

AMR support
-----------

Our new approach must also be extended to allow reading files that
contain full AMR hierarchies; currently, while Cello can write AMR
data, it is limited to reading only unigrid data files in the MUSIC
format. AMR support will be used for both data and checkpoint files.

File content
------------

Lastly, the content of the data files must be augmented to include all
data required to recreate a previously saved AMR block array. This is
required for both checkpoint/restart. Current data dump files do not
include all required data, such as block connectivity, scalar
data, and other Block attributes used internally.

Parallelization
---------------

Another difference between the proposed and current designs are how
individual tasks in an I/O operation are scheduled and synchronized.
In the current approach, a subset of root blocks are defined to be
"writers", and they control requesting and writing data from other
blocks within their scope. In the new approach, a separate ``IoWrite``
or ``IoRead`` chare array will be used instead, with different writers
or readers used for different file types, e.g. ``IoReadCheck`` or
``IoWriteData``.  Advantages are better load-balancing of I/O
operations and decoupling I/O operations from the Block chare array,
and allow multiple file operations to run concurrently. For EnzoBlock
data, ``IoEnzoReader`` and ``IoEnzoWriter`` chare arrays are used.

======
Design
======

Components of the new I/O approach include

  1. Control management
     InitialInput
     MethodOutput
  2. Control definition
  3. Block-to-file mapping
  4. File type definition
  5. File operation implementation

----------
Interfaces
----------

InitialRestart class
--------------------

.. glossary::

   ``MethodCheckpoint::MethodCheckpoint()``
   *Create a new MethodCheckpoint object*

----


MethodOutput class
------------------

.. glossary::

   ``MethodCheckpoint::MethodCheckpoint()``
   *Create a new MethodCheckpoint object*

----

MethodOutput class
------------------

----

.. glossary::

   ``MethodCheckpoint::MethodCheckpoint()``
   *Create a new MethodCheckpoint object*

=============
Communication
=============

=======
Testing
=======
