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
   :maxdepth: 1

   principles
   parameters-list
   getting_started
   parameters-example
   parameters-file
   design

An Enzo-P / Cello User and Developer Guide is available from the link below.  (Warning, it's big (0.2GB).  And mostly outdated, though I will be updating it over the coming weeks.)

   :download:`Using and Developing Enzo-P/Cello <./enzo-p-cello.pdf>`

James Bordner
jobordner@ucsd.edu


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

[Modified 2014-07-24]
