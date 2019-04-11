.. Enzo-E/Cello documentation master file, created by
   sphinx-quickstart on Thu Sep  1 17:54:50 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |date| date::

Enzo-E / Cello Documentation
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
	   
   index-tutorial
   index-use
   index-develop
   index-reference	      
James Bordner
jobordner@ucsd.edu


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

