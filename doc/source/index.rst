.. Enzo-E/Cello documentation master file, created by
   sphinx-quickstart on Thu Sep  1 17:54:50 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |date| date::

Enzo-E/Cello Documentation
============================

**Cello** is a highly scalable, fully-distributed array-of-octree
parallel adaptive mesh refinement (AMR) framework, and **Enzo-E** is a
scalable branch of the `ENZO <http://enzo-project.org/>`_ parallel
astrophysics and cosmology application that has been ported to use
Cello.  Enzo-E / Cello is currently funded by the National Science
Foundation (NSF) grant OAC-1835402, with previous funding through
NSF grants PHY-1104819, AST-0808184, and SI2-SSE-1440709.

Two fundamental differences between Enzo-E and ENZO are their AMR
design and code parallelization.  Cello implements *array of octree* AMR,
which has demonstrated scalability to date through 256K floating-point
cores of the `NSF <http://www.nsf.gov/>`_ `Blue Waters supercomputer
<http://www.ncsa.illinois.edu/enabling/bluewaters>`_ at the `National
Center for Supercomputing Applications
<http://www.ncsa.illinois.edu/>`_.  Unlike ENZO, which is parallelized
using MPI, Enzo-E/Cello is parallelized using `Charm++
<http://charm.cs.uiuc.edu/research/charm/>`_, an OOP parallel
programming system, targeting the development of Exascale software
applications, and actively developed at the Parallel Programming
Laboratory at the University of Illinois, Urbana-Champaign.

Enzo-E currently has two hyperbolic solvers: `PPM
<http://adsabs.harvard.edu/abs/1995CoPhC..89..149B>`_, an enhanced
piecewise parabolic method solver that was migrated to Enzo-E from the
ENZO code base, and `PPML <http://arxiv.org/abs/0905.2960>`_, an ideal
compressible MHD solver originally implemented in serial Fortran.
More recently, physics and infrastructure capabilities have been
developed for particle methods, including an implementation of ENZO's
CIC particle-mesh gravity solver, and cosmological expansion with
comoving coordinates.  Currently we are collaborating with
`Prof. Daniel Reynolds <http://faculty.smu.edu/reynolds/>`_ on
developing and implementing a highly scalable multigrid-based linear
solver.

.. toctree::
   :maxdepth: 1
	   
   tutorial/index
   user/index      
   reference/index
   devel/index
   design/index
   project/index

Getting started
---------------

The `Getting started using Enzo-E`_ section should cover everything
one needs to know to download Enzo-E / Cello and its dependent
software, configure, build, and run a small example test problem.

.. _Getting started using Enzo-E: getting_started.html

Parameters
----------

The `Enzo-E / Cello parameter reference`_ section is a reference page
for all Enzo-E and Cello parameters, which are used to write
parameters files defining simulations to run.

.. _Enzo-E / Cello parameter reference: parameters-list.html

How parameters are organized within a parameter file is described in
`Parameter files`_, which covers parameter groups, subgroups,
parameters, and data types recognized.

.. _Parameter files: parameters-file.html

An example parameter file is described in detail in `Parameter
file example`_.

.. _Parameter file example: parameters-example.html

Using and developing Enzo-E / Cello
-----------------------------------

The new `Using Enzo-E`_ section will describe in detail how to use
Enzo-E, including what methods, fields, particle types etc. are
available.  While only an outline exists now, we are working on
getting this section completed soon.

.. _Using Enzo-E: using.html

The new `Developing with Cello`_ section will describe how to add new
functionality to an application such as Enzo-E built on Cello,
including adding new computational methods, initial conditions, and
refinement criteria.  As with the `Using Enzo-E`_ section, this
documentation is under active development, and is mostly an outline at
this point.

.. _Developing with Cello: devel.html


PDF User and Developer Guide (depreciated)
------------------------------------------

An Enzo-E / Cello User and Developer Guide is available from the link
below.  Warning, it's big (0.2GB), and currently written as a
presentation, and some sections are somewhat outdated.  This document,
while currently still useful, is being phased out in favor of the
above online content.

   :download:`Using and Developing Enzo-E/Cello <http://client65-90.sdsc.edu/cello-doc/_downloads/4e523174654f2009e024ac3aa718a5d1/enzo-p-cello.pdf>`

Presentations given at the 2018 Enzo Days Workshop are also available
below.  They are more up-to-date, but are also still rather large.

   :download:`Enzo-E / Cello Status and What's New? <http://client65-90.sdsc.edu/cello-doc/_downloads/11ba173e00a03e5688c2d351f2476881/1805-1-status.pdf>`

   :download:`Introduction to Enzo-E <http://client65-90.sdsc.edu/cello-doc/_downloads/014af2086ca3c5899855903ea7d79d78/1805-2-intro.pdf>`
	     
   :download:`First Steps with Enzo-E <http://client65-90.sdsc.edu/cello-doc/_downloads/0bb42a270c3e64905ffce6517df0e971/1805-3-using.pdf>`

James Bordner
jobordner@ucsd.edu


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

