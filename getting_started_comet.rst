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

First, as in the general "getting started" page, the ``CELLO_ARCH``
and ``CELLO_PREC`` environment variables must be initialized.
``CELLO_PREC`` can be set to ``single`` or ``double`` as usual, and
``CELLO_ARCH`` shoulde be set to ``gordon_gnu``.  This tells the build
system to use the ``config/gordon_gnu.py`` configuration file, which
is used for Comet as well as Gordon.

.. code-block:: bash

   $ export CELLO_PREC=single
   $ export CELLO_ARCH=gordon_gnu

Secondly, the ``HDF5HOME`` environment variable should be initialized
to tell the build system where the HDF5 library is installed.  This
can be set to ``/home/ux452912/public``, which is a publicly readable
page in my home directory.  Or you can compile HDF5 yourself in your
own home directory and set ``HDF5HOME`` accordingly.

.. code-block:: bash

   $ export HDF5HOME=/home/ux452912/public

-------
Modules
-------

On Comet, as on many other supercomputers, `Modules`_ are used to aid
in initializing the shell environment.  Modules useful for Enzo-E / Cello
are ``python``, ``boost``, and ``gnu``.  To load these modules, use the
following:

.. code-block:: bash

   $ module load python boost gnu
   
.. _Modules: https://www.sdsc.edu/support/user_guides/comet.html#modules

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

After compiling Enzo-E, an extra step beyond those in the general
`Getting started using Enzo-E`_ page are required to run Enzo-E on
Comet, due to the HDF5 library.

Above, we set the ``HDF5HOME`` environment variable to
`/home/ux452912/public`` when compiling Enzo-E.  To run Enzo-E, it
needs the ``LD_LIBRARY_PATH`` environment variable to include this
path (+ '/lib').  However, Charm++ does not automatically pass on the
``LD_LIBRARY_PATH`` to the application, so we have to do that
explicitly using a short run script.  First, create a short script
``run-enzoe.sh`` containing the following:

.. code-block:: csh

    #!/bin/csh
    setenv LD_LIBRARY_PATH /home/ux452912/public/lib:$LD_LIBRARY_PATH
    $*

Now, to run Enzo-E in either interactive mode or in a batch script,
use the following, modifying the number of processors (``8``) and
parameter file (``input/method_ppm-8.in``) accordingly.

..  code-block:: bash
		 
    $ ./charmrun ++runscript ./run-enzoe.sh +p8 bin/enzo-p input/method_ppm-8.in


This assumes that the charmrun command is in your path--if it is not,
then you will need to include the path name as well.
