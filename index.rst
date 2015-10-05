.. Enzo-P/Cello documentation master file, created by
   sphinx-quickstart on Thu Sep  1 17:54:50 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |date| date::

Enzo-P/Cello Documentation
============================

**Cello** is a highly scalable, fully-distributed forest-of-octree
adaptive mesh refinement (AMR) framework, and **Enzo-P** is an
Eularian hydrodynamics and MHD astrophysics application that is built
using Cello.  Enzo-P is a branch of the `Enzo
<http://enzo-project.org/>`_ parallel astrophysics and cosmology
application.  Enzo-P and Cello is currently funded by
the National Science Foundation (NSF) grant SI2-SSE-1440709.
Previous funding was through NSF grants PHY-1104819 and AST-0808184.

Two fundamental differences between Enzo-P and Enzo are its AMR design
and parallelization.  Cello implements *forest of octree* AMR, which
has demonstrated scalability to date through 32K floating-point cores
of the `NSF <http://www.nsf.gov/>`_ `Blue
Waters supercomputer
<http://www.ncsa.illinois.edu/enabling/bluewaters>`_ at NCSA.  Unlike Enzo,
which is parallelized using MPI, Enzo-P/Cello is parallelized using `Charm++
<http://charm.cs.uiuc.edu/research/charm/>`_, an externally-developed
OOP parallel language targeting the implementation of Exascale
applications.

Enzo-P currently has two hyperbolic solvers: `PPM
<http://adsabs.harvard.edu/abs/1995CoPhC..89..149B>`_, an enhanced
piecewise parabolic method solver that was migrated to Enzo-P from the
Enzo code base, and `PPML <http://arxiv.org/abs/0905.2960>`_, an ideal
compressible MHD solver originally implemented in serial Fortran.
Additional physics and infrastructure capabilities are currently being
developed, including self-gravity and support for particle methods.

.. toctree::
   :maxdepth: 1

   principles
   getting_started
   parameters-example
   parameters-file
   parameters-list
   design

**More complete documentation in the form of a large-scale PDF
presentation is currently being developed, and will be available here
about 2015-10-01.**

James Bordner
jobordner@ucsd.edu


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

[Modified 2014-07-24]
