.. Enzo-P/Cello documentation master file, created by
   sphinx-quickstart on Thu Sep  1 17:54:50 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |date| date::

Enzo-P/Cello Documentation
============================

**Cello** is a highly scalable, fully-distributed array-of-octree
parallel adaptive mesh refinement (AMR) framework, and **Enzo-P** is an
Eulerian hydrodynamics and MHD astrophysics application that is built
using Cello.  Enzo-P is a branch of the `Enzo
<http://enzo-project.org/>`_ parallel astrophysics and cosmology
application.  Enzo-P / Cello is currently funded by
the National Science Foundation (NSF) grant SI2-SSE-1440709,
with previous funding through NSF grants PHY-1104819 and AST-0808184.

Two fundamental differences between Enzo-P and Enzo are their AMR
design and parallelization.  Cello implements *array of octree* AMR,
which has demonstrated scalability to date through 256K floating-point
cores of the `NSF <http://www.nsf.gov/>`_ `Blue Waters supercomputer
<http://www.ncsa.illinois.edu/enabling/bluewaters>`_ at the `National
Center for Supercomputing Applications
<http://www.ncsa.illinois.edu/>`_.  Unlike Enzo, which is parallelized
using MPI, Enzo-P/Cello is parallelized using `Charm++
<http://charm.cs.uiuc.edu/research/charm/>`_, an externally-developed
OOP parallel programming system targeting the implementation of Exascale
applications.

Enzo-P currently has two hyperbolic solvers: `PPM
<http://adsabs.harvard.edu/abs/1995CoPhC..89..149B>`_, an enhanced
piecewise parabolic method solver that was migrated to Enzo-P from the
Enzo code base, and `PPML <http://arxiv.org/abs/0905.2960>`_, an ideal
compressible MHD solver originally implemented in serial Fortran.
More recently, physics and infrastructure capabilities have been
developed for particle methods, including an implementation of Enzo's
CIC particle-mesh gravity solver.  Currently we are collaborating with
`Prof. Daniel Reynolds <http://faculty.smu.edu/reynolds/>`_ on
developing and implementing a highly scalable multigrid-based linear
solver.

.. toctree::
   :hidden:

   getting_started
   parameters-list
   parameters-example
   parameters-file
   using      
   devel
   principles
   design

Getting started
---------------

The `getting started`_ section covers everything you need to know to download Enzo-P / Cello and its dependent software, configure, port, build, and run an example test problem.

.. _getting started: getting_started.html

Parameters
----------

The `parameters list`_ section is a reference page for all Enzo-P and
Cello parameters, which are used to write parameters files defining
simulations to run.

.. _parameters list: parameters-list.html

How parameters are organized within a parameter file is described in
`parameter files`_, which covers parameter groups, subgroups,
parameters, and data types recognized.

.. _parameter files: parameters-file.html

An example parameter file is described in its entirety in `parameter
file example`_.

.. _parameter file example: parameters-example.html

Using Enzo-P
------------

The new `Using Enzo-P`_ section will describe in detail how to use
Enzo-P, including what methods, fields, particle types etc. are
available.  While only an outline exists now, we are working on
getting this section completed soon.

.. _Using Enzo-P: using.html

Developing with Cello
---------------------

The new `Developing with Cello`_ section will describe how to add new
functionality to an application such as Enzo-P built on Cello,
including adding new computational methods, initial conditions, and
refinement criteria.  As with the `Using Enzo-P`_ section, this
documentation is under active development, and is mostly an outline at
this point.

.. _Developing with Cello: devel.html


PDF User and Developer Guide
----------------------------

An Enzo-P / Cello User and Developer Guide is available from the link
below.  Warning, it's big (0.2GB), and currently written as a
presentation.  Parts are also a little outdated, though I will be
transferring content to the above online "Using Enzo-P" and
"Developing with Cello" sections over the coming weeks, and will
update things along the way.

   :download:`Using and Developing Enzo-P/Cello <./enzo-p-cello.pdf>`

James Bordner
jobordner@ucsd.edu


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

[Modified 2014-07-24]
