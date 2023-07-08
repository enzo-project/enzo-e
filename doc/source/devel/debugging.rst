******************
Tips for Debugging
******************

This page contains a collection of useful tips for debugging.

General Advice
==============

In general, if you find simple calls to ``CkPrintf`` are inadequate, you should consider whether the GDB debugger can help.
A common scenario where the GDB debugger is useful is when you encounter a segmentation fault (this can occur when you dereference a ``NULL`` pointer).
Essentially, the debugger runs the program until it crashes and then it will let you inspect variables in the stack frames just before the segmentation fault occurred.
It can also be very useful to inspect variables in stack frames when the program aborted early during a call to the ``ERROR`` or ``ASSERT`` macros.

Debugging a program with GDB is easiest when debugging a problem: replicated in a version of Enzo-E/Cello that was built on top of a netlrts-based build of charm++ on a local machine with xterm windows.
In this scenario, you simply need to append ``++debug-no-pause`` to the arguments of the ``charmrun`` launcher (e.g. near the arguments specifying the number of nodes).
This will cause charm++ to open an xterm window for each node and within each window the program is launched under GDB.


.. note::
   To make optimal use of GDB, you should be sure that you compiled Enzo-E/Cello with debugger symbols.
   Usually, setting ``-DCMAKE_BUILD_TYPE=RELWITHDEBINFO`` during the build process is adequate, however some compiler optimizations will prevent you from querying the values of some local variables.
   This can be addressed by recompiling with ``-DCMAKE_BUILD_TYPE=DEBUG``.

   You generally don't need to worry about having symbols for other libraries (like charm++).

.. note::
   TODO: add instructions for using GDB on a remote cluster using an MPI-based build of charm++

   TODO: it may be useful to provide a link to documentation describing how to use GDB to investigate stack frames

   TODO: it may be useful to talk a little about coredumps. Specifically, how do we enable them and how do we use gdb to inspect them?


Iterating on a bug captured by a unit test
==========================================

When addressing a bug in in the Cello Layer that can be captured by a unit test, it can be a lot faster to simply rebuild the file containing the unit test each time you make a modification, rather than rebuilding all targets.

To be concrete, imagine we were fixing a bug in a subclass of :cpp:class:`Schedule`, that was reproduced in the ``test_schedule`` unit-test.
Each time we want to check whether a change to a source/header file fixes the bug, we can rebuild ``test_schedule`` by executing ``make test_schedule`` (or ``ninja make_schedule``, depending on how your build is set up) from your build directory.
In comparison, calling ``make`` (or ``ninja``) would rebuild all binaries and can take a lot more time.

Dumping data to disk
====================

For certain classes of problems, it can be useful to dump data (like field data) to disk (exactly when the problem arises) and investigate the values at a later time

We provide 2 functions for doing this:

.. doxygenfunction:: disk_utils::dump_view_to_hdf5

.. doxygenfunction:: disk_utils::dump_views_to_hdf5
