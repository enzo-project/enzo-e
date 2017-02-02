======
Design
======

The software design of Enzo-P / ``Cello`` is based on object-oriented
programming (OOP), a proven software design and implementation
paradigm that helps reduce software complexity and improves code
maintainability.

.. image:: components.png

The ``Cello`` framework is functionally decomposed into a collection
of software components.  The ``Enzo-P`` application is implemented on
top of Cello, which is comprised of the components within the central
boxed region.  ``Cello`` is in turn implemented using the Charm++
parallel programming system.  Components in ``Cello`` are organized
into multiple layers, with software dependencies directed primarily
downwards; cross-cutting software concerns are also implemented as
components.  Each component is composed of one or more inter-operating
C++ classes.

Enzo-P
======

``Enzo-P`` is the astrophysics and cosmology application built on top
of the ``Cello`` scalable AMR framework.  ``Enzo-P`` interfaces with
``Cello`` primarily through C++ class inheritance of ``Cello`` classes
in the ``Simulation``, ``Problem``, ``Method``, and ``Mesh``
components.  Due to the separation of concerns between ``Enzo-P`` and
Cello, development of ``Enzo-P`` requires no knowledge or awareness of
parallel programming, and development of ``Cello`` requires no
knowledge of specific numerical methods.  This allows relatively
independent software development of ``Enzo-P`` by physical science and
numerical methods experts, and of ``Cello`` by computer science and
parallel programming experts.

High-level components
=====================

Top level components of ``Cello`` include ``Simulation``, ``Problem``,
and ``Method``.  A ``Simulation`` defines and manages a computational
astrophysics ``Problem``, which defines physics, method, and data
structure ``Parameters``, initial conditions, and boundary conditions.
The ``Simulation`` class, which is implemented as a ``Charm`` *process
group*, initializes and begins parallel execution of the simulation.
``Method`` classes are used to implement individual numerical methods
that compute on ``Field`` (current) data on a block, as well as
``Particle`` and ``Ray`` data (proposed).

Data structure components
=========================

Data structure components include the existing ``Mesh`` and ``Field``
components, and the proposed ``Particle`` and ``Ray`` components, for
implementing the distributed computational data containers. The
``Mesh`` component includes classes for representing and operating on
an adaptive mesh hierarchy, implemented as a fully-distributed
array-of-octrees.  Octree nodes are associated with Blocks, which are
containers for the ``Field`` and other data, and implemented as a
Charm++ *chare array*.


Middle-level components
=======================

Middle-level components include ``Control``, ``Parameters`` and ``Io``
``Control`` handles the time stepping of ``Method`` s to advance the
problem forward in time, as well as sequencing adaptive mesh
refinement data structure operations, including remeshing, scheduling
dynamic load balancing (which will be delegated to Charm++), and
refreshing ghost zones on Block boundaries.  The ``Parameters``
component serves to read, store, and provide access to parameters
defined in an input configuration file.  To improve usability over
Enzo, configuration files are more structured, and support
floating-point and logical expressions to greatly simplify
initializing problems with complex initial conditions.  The ``Io``
component serves as a layer to coordinate the disk output of data
structure components, such as ``Simulation`` ``Hierarchy`` and
``Field`` data.  It calls the ``Disk`` component to handle actual file
operations.

Hardware-interface components
=============================

The lower-level hardware-interface components include ``Disk``,
``Memory`` and ``Parallel``.  The ``Disk`` component implements basic
disk operations, isolating the specific file format from the
higher-level ``Io`` component.  ``Disk`` currently supports HDF5, and
we propose to support the Adaptable IO System (ADIOS) in the future to
enhance transfer of data to and from other HPC software
components. The ``Memory`` component controls dynamic memory
allocation and management.  Currently ``Memory`` handles allocating
and monitoring heap memory usage; proposed functionality includes
allocating, deallocating, and transferring data between main memory,
hardware accelerator (GPU) memory, and many-core coprocessors
(e.g.~the Intel Xeon Phi).  As with the ``Disk`` component, this
serves to isolate lower-level details from higher-level components.
The ``Parallel`` component currently supplies basic access to core
rank and core count, and is being depreciated.

Interface components
====================

Interface components include ``Monitor`` (current) and ``Portal``
(proposed).  The ``Monitor`` component controls the user-readable
summary of progress to stdout, and the proposed ``Portal`` component
will control the interaction of ``Enzo-P`` with external applications
running concurrently, such as inline analysis or real-time
visualization.  One particular such analysis and visualization
application is yt, which we will use to help drive the design and
development of the ``Portal`` component.

Cross-cutting components
========================

Some ``Cello`` components can in principle be called from any software
layer---these include ``Performance`` and ``Error``. The
``Performance`` component dynamically collects performance data for
the running ``Enzo-P`` simulation, and provides a holistic summary of
performance data to the user, as well as to software components that
can adapt to optimize desired performance metrics.  Current metrics
measured include memory usage (via the ``Memory`` component), and
computation amount and memory access amount (via the Performance
Application Programming Interface (PAPI).  Future support will include
metrics for monitoring parallel communication, dynamic load balancing,
and disk usage.  The ``Error`` component will be used to detect,
evaluate, and decide what to do about software errors; higher-level
error detection and recovery will be handled by Charm++, which
supports both simple checkpoint to disk, as well as double in-memory
checkpoint with automatic restart.

