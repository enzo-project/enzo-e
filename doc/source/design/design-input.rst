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

The implementation will include modifications to how blocks are mapped
to files, what data are written to the files, and how file I/O is
parallelized.

Ordering
--------

Our new approach will involve a generalization of the ``MethodOutput``
method, but will enable load-balanced disk file through the use of
block `orderings` to define how blocks are mapped to files. Currently,
the ordering used in ``MethodOutput``, which is implicit and embedded
in the code, is based on a regular partitioning of root-level blocks
together with their descendents. The updated implementation will
factor out this ordering into an ``Ordering`` class, and allow
flexibility in defining other orderings, such as space-filling curves
like Hilbert or Morton.

File content
------------

The content of the data files must be augmented to include all data
required to recreate a previously saved AMR block array on restart.
Current data dump files do not include all required data, such as
block connectivity, scalar data, and other Block attributes used
internally.

Control flow
------------

Another difference between the proposed and current designs are how
individual tasks in an I/O operation are scheduled and synchronized.
In the current approach, a subset of root blocks is defined to be
"writers" or "readers", and they control requesting and writing data
from other blocks within their assigned scope. In the new approach, a
separate ``IoWriter`` or ``IoReader`` chare array will be used
instead.  Advantages are better load-balancing of I/O operations, and
decoupling of I/O operations from the Block chare array, which could
improve load balancing. For Enzo-E checkpoint/restart data,
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
     * ``Ordering``
     * ``OrderingHilbert``
     * ``OrderingRootBlocks``

----------
Algorithms
----------

Input algorithm
---------------

..
   .. image:: io-input.png

..
   .. image:: io-read.png

.. uml::
   @startuml

   box "Simulation array" #Pink
   participant simulation_0 as s0
   participant simulation_k as sk
   end box

   box "IoReader array" #LightBlue
   participant reader_0 as r0
   participant reader_1 as r1
   end box

   box "Block array" #LightGreen
   participant block_0 as b0
   participant "block_k-1" as bk
   participant "block_k" as bkp1
   end box

   == startup ==

   b0 -> s0 : p_restart_enter()
   s0 -> r0 ** : ckNew()
   s0 -> r1 **
   s0 <-> sk : p_set_io_reader()
   hnote over s0,sk : r_restart_start()

   == level 0 ==

   s0 -> r0 : p_init_root()
   s0 -> r1
   r0 ->o b0 : p_restart_set_data()
   r1 ->o b0
   b0 o-> r0 : p_block_ready()
   b0 o-> r1
   hnote over r0,r1 : sync
   r0 -> s0 : r_restart_next_level()
   r1 -> s0
   hnote over s0 : sync

   == level k ==
   loop for k=1 to L
   s0 -> bk **
   s0 -> r0 : p_create_level(k)
   s0 -> r1
   r0 ->o bk : p_restart_refine()
   r1 ->o bk
   bk o-> bkp1 ** : insert
   bkp1 o-> r0 : p_block_created()
   bkp1 o-> r1

   hnote over r0,r1 : sync
   r0 -> s0 : r_restart_level_created()
   r1 -> s0
   hnote over s0 : sync
   s0 -> r0 : p_init_level(k)
   s0 -> r1
   r0 ->o bkp1 : p_restart_set_data()
   r1 ->o bkp1
   bkp1 o-> r0 : p_block_ready()
   bkp1 o-> r1
   hnote over r0,r1 : sync
   r0 -> s0 : r_restart_next_level()
   r1 -> s0
   hnote over s0 : sync
   end
   == cleanup ==
   s0 -> s0 : doneInserting()
   s0 -> r0 : delete
   destroy r0
   s0 -> r1
   destroy r1
   s0 ->o b0 : p_restart_done()
   s0 ->o bk
   s0 ->o bkp1

   @enduml

// Begin reading restart data and create the mesh hierarchy of
// EnzoBlocks. Replaces functionality of control_adapt.

.. code-block:: C++

    entry void main::r_restart_enter(std::string file_hierarchy)
    {
       // open hierarchy file
       // read hierarchy file
       // create block_array
       // create IoEnzoReader array
       // initialize sync_file(num_io_reader)
       for (i_f = files in restart) {
          io_reader[i_f].insert(file_name[i_f]);
       }
       // close hierarcy file
    }

.. code-block:: C++

    // Read in a Block file and create all blocks in the file
    IoEnzoReader::IoEnzoReader(file_block)
    {
       // open block file
       // read block file
       // initialize sync_block(num_blocks)
       for (i_b = loop over blocks) {
          // read Block data
          // create Block, initialize, notify caller when done
          enzo_factory.create_block(block_data, io_reader, i_f);
       }
       // close file
    }

.. code-block:: C++

    // Initialize a block then notify caller when done
    Block::Block (block_data, io_reader, i_f)
    {
       // initialize Block using block_data
       io_reader[i_f].p_block_created(index);
    }

.. code-block:: C++

    // After all blocks created, notify main when done
    entry IoReader::p_block_created(Block index)
    {
       // Count blocks created
       if (sync_block.done()) {
          // After last block created, exit restart phase
          main::p_restart_done();
       }
    }

.. code-block:: C++

    // After all files have been read, proceed with the simulation
    entry main::p_restart_done()
    {
       if (sync_file.done()) {
          restart_exit();
       }
    }

    void main::restart_exit()
    {
       // exit restart phase
    }

Output algorithm
----------------

.. image:: io-output.png

.. code-block:: C++

    // Begin writing a checkpoint file
    entry void EnzoMethodCheck::apply(Block)
    {
       contribute (main::r_check_enter(std::string file_hierarchy));
    }

.. code-block:: C++

    entry void main::r_check_enter(std::string file_hierarchy)
    {
       // create hierarchy file
       // write hierarchy file
       // create IoEnzoWriter array
       // initialize sync_writer(num_io_writer)
       for (i_f = files in checkpoint) {
          io_writer[i_f].insert(file_name[i_f]);
       }
       // [wait for IoWriters to all be created before proceeding]
    }

.. code-block:: C++

   IoEnzoWriter::IoEnzoWriter()
   {
      // open block file
      // notify main that this IoWriter has been created
      main.p_writer_created();
   }

.. code-block:: C++

   entry void main::p_writer_created()
   {
      // after all writers check-in, get first blocks
      if (sync_writer.done()) {
         // initialize sync_file
         // ask blocks to self-identify as the first block in a file
         block_array.p_write_first();
      }
   }

.. code-block:: C++

    entry void EnzoBlock::p_write_first(io_writer, ordering)
    {
       // if this block is first in the partitioned ordering,
       // send data to assigned writer i_w
       if (ordering.is_start(index, count)) {
          write_next(io_writer,i_w);
       }
    }

    void EnzoBlock::write_next(io_writer, i_w)
    {
      // pack data
      // determine next Block in the ordering, else signal if last
      // send data to assigned io writer to output to file
      io_writer[i_w].p_write(data_buffer, index_next, is_last)
    }

.. code-block:: C++

    entry void IoEnzoWrite::p_write(data_buffer,index_next, is_last);
    {
       // unpack block data
       // write block data to file
       // request next block if any, else signal main we're done
       if (is_last) {
          // close file
          main.p_check_done();
       } else {
          block_array[index_next].p_write_next();
       }
    }

.. code-block:: C++

    entry void EnzoBlock::p_write_next(io_writer, i_w)
    {
       write_next(io_writer, i_w);
    }

.. code-block:: C++

    // After all files have been written, proceed with the simulation
    entry void main::p_check_done()
    {
       if (sync_file.done()) {
          check_exit();
       }
    }

    void main::check_exit() {
    {
       // exit checkpoint phase
    }

-------
Classes
-------

EnzoMethodInput

=============
Communication
=============

=======
Testing
=======
