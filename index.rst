.. Enzo-P/Cello documentation master file, created by
   sphinx-quickstart on Thu Sep  1 17:54:50 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |date| date::

Enzo-P/Cello Documentation
============================

**Cello** is a highly scalable, fully-distributed array-of-octree
parallel adaptive mesh refinement (AMR) framework, and **Enzo-P** is a
scalable branch of the `ENZO <http://enzo-project.org/>`_ parallel
astrophysics and cosmology application that has been ported to use
Cello.  Enzo-P / Cello is currently funded by the National Science
Foundation (NSF) grant SI2-SSE-1440709, with previous funding through
NSF grants PHY-1104819 and AST-0808184.

Two fundamental differences between Enzo-P and ENZO are their AMR
design and code parallelization.  Cello implements *array of octree* AMR,
which has demonstrated scalability to date through 256K floating-point
cores of the `NSF <http://www.nsf.gov/>`_ `Blue Waters supercomputer
<http://www.ncsa.illinois.edu/enabling/bluewaters>`_ at the `National
Center for Supercomputing Applications
<http://www.ncsa.illinois.edu/>`_.  Unlike ENZO, which is parallelized
using MPI, Enzo-P/Cello is parallelized using `Charm++
<http://charm.cs.uiuc.edu/research/charm/>`_, an OOP parallel
programming system, targeting the development of Exascale software
applications, and actively developed at the Parallel Programming
Laboratory at the University of Illinois, Urbana-Champaign.

Enzo-P currently has two hyperbolic solvers: `PPM
<http://adsabs.harvard.edu/abs/1995CoPhC..89..149B>`_, an enhanced
piecewise parabolic method solver that was migrated to Enzo-P from the
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
	   
   getting_started
   parameters-list
   parameters-file
   parameters-example
   using      
   devel
   principles
   design

Getting started
---------------

The `Getting started using Enzo-P`_ section should cover everything
one needs to know to download Enzo-P / Cello and its dependent
software, configure, build, and run a small example test problem.

.. _Getting started using Enzo-P: getting_started.html

Parameters
----------

The `Enzo-P / Cello parameter reference`_ section is a reference page
for all Enzo-P and Cello parameters, which are used to write
parameters files defining simulations to run.

.. _Enzo-P / Cello parameter reference: parameters-list.html

How parameters are organized within a parameter file is described in
`Parameter files`_, which covers parameter groups, subgroups,
parameters, and data types recognized.

.. _Parameter files: parameters-file.html

An example parameter file is described in detail in `Parameter
file example`_.

.. _Parameter file example: parameters-example.html

Using and developing Enzo-P / Cello
-----------------------------------

The new `Using Enzo-P`_ section will describe in detail how to use
Enzo-P, including what methods, fields, particle types etc. are
available.  While only an outline exists now, we are working on
getting this section completed soon.

.. _Using Enzo-P: using.html

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
presentation, and some sections are somewhat outdated.  This document,
while currently still useful, is being phased out in favor of the
above online content.

   :download:`Using and Developing Enzo-P/Cello <./enzo-p-cello.pdf>`

James Bordner
jobordner@ucsd.edu


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

[Modified 2017-11-30]
