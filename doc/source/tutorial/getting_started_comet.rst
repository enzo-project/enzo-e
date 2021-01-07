.. _Comet:

Getting Started: SDSC Comet
===========================

This page is for getting started using Enzo-E / Cello on the "Comet"
supercomputer at the San Diego Supercomputer Center.  For more general
help in getting Enzo-E compiled and running, see the `Getting started
using Enzo-E`_ page; for general information about using the Comet
supercomputer, see the `Comet User Guide`_.

.. _Getting started using Enzo-E: getting_started.html

.. _Comet User Guide: https://www.sdsc.edu/support/user_guides/comet.html

Install ``Charm++`` on Comet
----------------------------

``Charm++`` can be `downloaded here <http://charm.cs.illinois.edu/software>`_.

To build Charm++ on Comet, try this:

   ``./build charm++ verbs-linux-x86_64 -j8  --with-production``

Charm++ also provides a helpful interactive script ``smart-build.pl``
that can be used to configure and install ``Charm++``.
   
Configuring on Comet
--------------------


---------------------
Environment variables
---------------------

Several environment variables should be initialized before compiling
and running Enzo-E / Cello on Comet.

First, as in the general "getting started" page, the ``CELLO_ARCH``,
``CELLO_PREC``, and ``CHARM_HOME`` environment variables must be
initialized.  ``CELLO_PREC`` can be set to ``single`` or ``double`` as
usual, and ``CELLO_ARCH`` should be set to ``comet_gnu``.  This tells
the build system to use the ``config/comet_gnu.py`` configuration
file, which is used for Comet.

.. code-block:: bash

   $ export CELLO_PREC=single
   $ export CELLO_ARCH=comet_gnu
   $ export CHARM_HOME=$HOME/Charm/charm.6.10

-------
Modules
-------

On Comet, as on many other supercomputers, `Modules`_ are used to aid
in initializing the shell environment.  Modules useful for Enzo-E /
Cello are ``python``, ``boost``, ``gnu``, and ``hdf5-serial``.  To
load these modules, you can use the following:

.. code-block:: bash

   $ module load python boost gnu hdf5-serial
   
.. _Modules: https://www.sdsc.edu/support/user_guides/comet.html#modules

You may want to add this line to your ``$HOME/.bashrc`` bash
startup file to avoid having to load them each time you log in.

You should now be able to compile Enzo-E using ``make``, as described
in the ``Building`` section of the `Getting started using Enzo-E`_
page.

--------------
Running Enzo-E
--------------

Enzo-E can be run on Comet using Comet's batch system, or in interactive
mode.  See the `Running Jobs on Regular Compute Nodes`_ section of
the Comet User Manual for how to run batch or interactive jobs.

.. _Running Jobs on Regular Compute Nodes: https://www.sdsc.edu/support/user_guides/comet.html#running

To run Enzo-E in either interactive mode or in a batch script,
use the following, modifying the number of processors (``8``) and
parameter file (``input/method_ppm-8.in``) accordingly.

..  code-block:: bash
		 
    $ ./charmrun +p8 ++mpiexec bin/enzo-p input/method_ppm-8.in


This assumes that the charmrun command is in your path--if it is not,
then you will need to include the path name as well.  Also note the
addition of "++mpiexec".

----

2020-04-10: Updated with corrections from Joshua Smith.
