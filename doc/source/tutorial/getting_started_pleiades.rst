.. _Pleiades:

Getting Started: NASA Pleiades/Electra/Aitken
=============================================

Latest updated/verification: 05 May 2021

This page is for getting started using Enzo-E / Cello on the Pleiades,
Electra and Aitken supercomputers at NASA's HECC.  For more general
help in getting Enzo-E compiled and running, see the `Getting started
using Enzo-E`_ page; for general information about using the NASA
supercomputers, see the `NASA HECC User Guide`_.

.. _Getting started using Enzo-E: getting_started.html

.. _NASA HECC User Guide: https://www.nas.nasa.gov/hecc/support/kb/


Software environment
--------------------

Note, this guide currently covers the Intel CPU nodes of Pleiades,
Electra, and Aitken.
The AMD EPYC Rome nodes use a different operating system that
requires loading different modules, see
`NASA HECC Aitken Rome Node documentation`_, and potentially
recompiling the software stack for different library versions
and architecture.

Load GNU environment with MPI backend and install
Python based``scons`` build tool required by Cello/Enzo-E

.. code-block:: bash

   $ module load gcc/8.2 boost/1.62 mpi-hpe/mpt python3/3.5.2
   $ pip install --user scons
   # Add to path so that scons is found 
   $ export PATH=${HOME}/.local/bin:$PATH

Note, building with Intel is not straightforward (and thus not shown/recommended
here) as ``scons`` does not properly detect NASA's Intel compiler modules.
Also note, the newer ``python3/3.7.0`` module causes some problems when
configuring ``Charm++`` so using ``3.5.2`` as above is recommended/required.

Now build HDF5 (required because the provided HDF5 modules are only supported with
Intel modules).
Working in the ``src`` directory is of course optional but done in this guide
for consistency.

.. code-block:: bash

   $ mkdir -p ~/src
   $ wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.gz
   $ tar xzf hdf5-1.10.7.tar.gz
   $ cd hdf5-1.10.7
   $ ./configure --prefix=${PWD}/install --enable-build-mode=production
   $ make -j8
   $ make install


.. _NASA HECC Aitken Rome Node documentation: https://www.nas.nasa.gov/hecc/support/kb/preparing-to-run-on-aitken-rome-nodes_657.html 


Install ``Charm++``
-------------------

``Charm++`` can be `downloaded here <http://charm.cs.illinois.edu/software>`_.

Now get and build Charm++ (version 6.10.2 was the latest one when writing this guide
and verified to work).

.. code-block:: bash

   $ cd ~/src
   $ wget http://charm.cs.illinois.edu/distrib/charm-6.10.2.tar.gz
   $ tar xzf charm-6.10.2.tar.gz
   $ cd charm-v6.10.2/
   $ ./build charm++ mpi-linux-x86_64 gcc gfortran --with-production --enable-tracing -j8

   
Install Cello/Enzo-E
--------------------


---------------------
Environment variables
---------------------

Several environment variables should be initialized before compiling
and running Enzo-E / Cello.

First, as in the general "getting started" page, the ``CELLO_ARCH``,
``CELLO_PREC``, and ``CHARM_HOME`` environment variables must be
initialized.  ``CELLO_PREC`` can be set to ``single`` or ``double`` as
usual, and ``CELLO_ARCH`` should be set to ``pleiades_intel``.  This tells
the build system to use the ``config/pleiades_intel.py`` configuration
file, which is used for Comet.

.. code-block:: bash

   $ export CELLO_PREC=double
   $ export CELLO_ARCH=pleiades_gcc
   $ export CHARM_HOME=$HOME/src/charm-v6.10.2

--------
Building
--------

One option to obtain the source is to directly clone from GitHub.
By default the most recent version of the ``master`` branch will be checked out.

.. code-block:: bash

   $ cd ~/src
   $ git clone https://github.com/enzo-project/enzo-e.git
   $ cd enzo-e
   # (OPTIONAL) if you do not have grackle installed edit config/pleiades_gcc.py
   #   and set use_grackle = 1 and the corresponding path at the bottom
   $ make 

If everything went well the build process finishes with

.. code-block:: bash
  
  Install file: "build/Enzo/enzo-p" as "bin/enzo-p"
  Success
  done
  END   Enzo-P/Cello ./build.sh: arch = pleiades_gcc  prec = double  target = bin/enzo-p time = 1.53 min
  Done.

--------------
Running Enzo-E
--------------

To test whether Enzo-E works across nodes an interactive session could be used.
For example

.. code-block:: bash

  # get an interactive job.
  $ qsub -I  -l select=1:ncpus=40:model=sky_ele -l walltime=0:30:0
  
  # now within the job first prepare the environment
  $ module load gcc/8.2 boost/1.62 mpi-hpe/mpt python3/3.5.2
  $ export LD_LIBRARY_PATH=${HOME}/src/hdf5-1.10.7/install/lib:$LD_LIBRARY_PATH

  # Now run HelloWorld test (using 4 processes here):
  $ mpiexec -np 4 ./bin/enzo-p input/HelloWorld/Hi.in


