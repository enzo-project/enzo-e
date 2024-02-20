.. _getting_started:


----------------------------
Getting started using Enzo-E
----------------------------

This page will help you get Enzo-E and Cello up and running.  It
covers downloading the source code, porting the code to new platforms,
configuring and compiling the code, and running a sample test problem.

Other pages are available for helping to get started on specific
architectures, including the "Frontera" supercomputer at TACC and the
"Pleiades" supercomputer at NASA.

.. toctree::
   :maxdepth: 1
	   
   getting_started_pleiades
   getting_started_frontera

Dependency Installation
=======================

Before compiling ``Enzo-E / Cello``, you may need to download
and install 1. ``CMake``, 2. ``Charm++``, 3. ``HDF5``, 4. ``libpng``, 5. ``libboost``, and (optionally) 6. ``Grackle``:

1. Install ``CMake``
--------------------

Most systems nowaways have ``CMake`` already installed.
If not, you can get the binary distribution from the
`CMake <https://cmake.org/download/>`_ website.


2. Install ``Charm++``
----------------------

We generally recommend that you download `Charm++`, from the GitHub repository and then retrieve the latest version (7.0.0) known to build `Enzo-E/Cello`.
This is demonstrated by the following command:

..  code-block:: bash

  git clone https://github.com/UIUC-PPL/charm.git
  cd charm
  git checkout v7.0.0

Alternatively, if you don't have `git` installed, you can also download the appropriate release version of the code from the main `website <http://www.charmplusplus.org/>`_ or from the `Release page on GitHub <https://github.com/UIUC-PPL/charm/releases>`_.

`Charm++` can be configured to use different "backends" for handling communication, which include:
  * ``mpi``: our recommended default choice for distributed machines.
    Under this configuration, communication occurs via calls to MPI commands (this obviously requires an MPI installation).
  * ``netlrts``: our recommended choice for development/testing on a local machine (this can make debugging far simpler).
    Under this configuration, communication occurs via standard networking protocols (that all modern computers are equipped with).
    This is generally slower than other options.

`Charm++` uses `cmake` as the default build system.
To build `Charm++` start from directory holding all downloaded files.

Invoke the following command to build with the `mpi` backend:

..  code-block:: bash

  mkdir build-mpi # this directory name is arbitrary
  cd build-mpi
  cmake -DNETWORK=mpi -DSMP=OFF ..
  make # you can specify make -j4 to compile with 4 threads

The following commands to build with the `netlrts` backend:

..  code-block:: bash

  mkdir build-netlrts # this directory name is arbitrary
  cd build-netlrts
  cmake -DNETWORK=netlrts -DSMP=OFF ..
  make # you can specify make -j4 to compile with 4 threads


As an aside `Charm++` can also be configured to directly use a machine's low-level communication primitives (that ``MPI`` implementations wrap).
It can be complicated to compile with these backends, so additional details are provided on an architechture-specific basis.


3. Install ``HDF5``
-------------------

"`HDF5 <https://www.hdfgroup.org/HDF5/>`_ is a "data model, library, and
file format for storing and managing data", and is the primary library
used by Enzo-E / Cello for data output.

If HDF5 is not already installed on your machine, it may be available
through your operating system distribution, otherwise it can be
downloaded from the `HDF5 <https://www.hdfgroup.org/HDF5/>`_ website.
Enzo-E / Cello currently uses the "serial" (non-MPI) version of HDF5.

4. Install ``libpng``
---------------------

"`libpng <https://www.libpng.org/pub/png/libpng.html>`_ is the official
PNG reference library", and is the image format used by Enzo-E / Cello.

If ``libpng`` is not already installed on your machine, it may be
available through your operating system distribution, otherwise it can
be downloaded from the `libpng
<https://www.libpng.org/pub/png/libpng.html>`_ website.

5. Install ``libboost-dev``
---------------------------

"`Boost <https://www.boost.org/>`_ provides free peer-reviewed portable C++ source libraries."

If ``libboost-dev`` is not already installed on your machine, it may be
available through your operating system distribution, otherwise it can
be downloaded from the `libboost <https://www.boost.org/>`_ website.

6. Install Grackle  (Optional)
------------------------------

By default, Enzo-E requires the Grackle chemistry and cooling library.
If you do not need to use Grackle, you can simple disabling it by setting
``-DUSE_GRACKLE=OFF`` when you configure Enzo-E.
See the `Grackle documentation <https://grackle.readthedocs.io>`__ for installation
instructions.

7. Install yt (Optional)
------------------------

If you want to use the yt python package to analyse Enzo-E output data, you should install the latest version from source.
This can be done with the following commands:

.. code-block:: bash

    git clone https://github.com/yt-project/yt.git
    cd yt
    pip install -e .

This is **NOT** a requirement for building and running Enzo-E, but it is used in some tests.

.. _how-to_configure_build:

Configuring/Building
====================

``Enzo-E / Cello`` is currently hosted on github.com (previously bitbucket.com).
To obtain the latest version of the source code, you may clone it from the repository `Enzo-E / Cello github repository <https://github.com/enzo-project/enzo-e.git>`_:

.. code-block:: bash

    git clone https://github.com/enzo-project/enzo-e.git

For a basic configuration on a linux system with the GNU compiler stack, the following command should work out of the box.
If not, please report any problems.
See the following subsection for more configuration options.

..  code-block:: bash

  mkdir my-first-build # Again, the directory name can be anything
  cd my-first-build
  # replace <PATH/TO/charm/build-dir> with the path to the directory created
  # when you built charm++ (either build-mpi or build-netlrts if you've been
  # following along)
  cmake -DCHARM_ROOT=<PATH/TO/charm/build-dir> \
        -DEnzo-E_CONFIG=linux_gcc -DUSE_GRACKLE=OFF ..
  make -j4 # -j4 tells make to execute up to 4 commands in parallel

To build on a Mac, you should only need to replace ``linux_gcc`` with ``darwin_clang``.

Note, if ``ninja`` is installed, the ``ninja`` build system can be used for faster build times.
This is done by adding ``-GNinja`` to the ``cmake`` command (before the ``..``) and calling ``ninja`` afterwards instead of ``make``.

The Enzo-E executable is built within ``bin/``.

Configuration options
---------------------

Current ``cmake`` options are listed in the following subsubsections.
Skip ahead to :ref:`how_to_specify_the_configuration` for details about how to specify the configuration.

General Configuration
^^^^^^^^^^^^^^^^^^^^^

These options are related to Enzo-E's general configuration.

.. list-table:: General Configuration
   :widths: 10 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``USE_DOUBLE_PREC``
     - Use double precision. Turn off for single precision.
     - ON
   * - ``new_output``
     - Temporary setting for using new Output implementation
     - OFF
   * - ``node_size``
     - Maximum number of procesess per shared-memory node (can be larger than needed)
     - 64
   * - ``smp``
     - Use Charm++ in SMP mode (Charm++ must have been compiled to support SMP mode).
     - OFF
   * - ``balance``
     - Enable charm++ dynamic load balancing
     - ON
   * - ``balancer_included``
     - Charm++ load balancer included
     - "CommonLBs"
   * - ``balancer_default``
     - Charm++ load balancer to use by default
     - "TreeLB"

Configuring Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

The parameters in the following table control whether ``Enzo-E / Cello`` should make use of certain external dependencies.
These options all assume that ``cmake`` is succesfully able to find the location of these dependencies.
In some cases, you may need to provide additional hints about the location of the dependencies (see :ref:`how_to_specify_the_configuration` for more details)

.. list-table:: External Dependencies
   :widths: 10 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``USE_GRACKLE``
     - Use Grackle Chemistry
     - ON
   * - ``use_jemalloc``
     - Use the jemalloc library for memory allocation
     - OFF
   * - ``use_papi``
     - Use the PAPI performance API
     - OFF

Floating-Point Optimization Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generally, ``Enzo-E / Cello`` will not be compiled with value-unsafe floating point optimizations unless they are explicitly enabled.
This default behavior was chosen to support automated tests that check symmetry preservation in certain components of the code (mostly related to hydro/MHD).
While this feature is not essential for all applications, it is very useful when debugging new features.

Note that the Intel compilers deviate from this behavior. By default, they enable more aggressive floating point optimizations by default (consequently the tests checking symmetry may fail). This is a large part of the reason that the Intel compilers have a reputation for producing faster code (in practice, this is not always the case).

.. list-table:: Floating point optimizations
   :widths: 10 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``OPTIMIZE_FP``
     - Enable value-unsafe floating point optimizations (for some compilers, this enables ``-ffast-math``).
     - OFF
   * - ``USE_SIMD``
     - Enables compiler support for OpenMP SIMD directives (for ``gcc``, ``icc``, and ``clang`` compilers this will NOT enable other OpenMP directives and should not link the openmp runtime library). Enabling this requires that ``OPTIMIZE_FP=ON``.
     - OFF


Profiling Options
^^^^^^^^^^^^^^^^^
These options control compilation choices that can be used to facillitate profiling in ``Enzo-E / Cello``

.. list-table:: Profile-Related Configuration
   :widths: 10 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``memory``
     - Track dynamic memory statistics.  Can be useful, but can cause problems on some systems that also override ``new [] ()`` / ``delete [] ()``
     - OFF
   * - ``use_gprof``
     - Compile with -pg to use gprof for performance profiling
     - OFF
   * - ``use_performance``
     - Use Cello Performance class for collecting performance data (currently requires global reductions, and may not be fully functional) (basic time data on root processor is still output)
     - ON
   * - ``use_projections``
     - Compile the CHARM++ version for use with the Projections performance tool.
     - OFF

Testing Options
^^^^^^^^^^^^^^^

These options configure properties of automated tests. These options currently just affect tests in the :ref:`ctest framework <ctest>` and don't affect tests in the :ref:`pytest framework <pytest>`.

.. list-table:: Testing-Related Configuration
   :widths: 10 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``PARALLEL_LAUNCHER``
     - Launcher to use for parallel tests
     - charmrun
   * - ``PARALLEL_LAUNCHER_NPROC_ARG``
     - Argument to set number of processing elements for parallel launcher (for use with ``charmrun``)
     - +p
   * - ``PARALLEL_LAUNCHER_NPROC``
     - Number of processors to run parallel unit tests
     - 4
   * - ``BUILD_TESTING``
     - Whether to setup the CTest infrastructure and build unit test binaries (which are primarily built to be executed by the CTest infrastructure). This has no effect on the pytest infrastructure.
     - "ON"


Debugging Options
^^^^^^^^^^^^^^^^^

The following options are useful for debugging.
       
.. list-table:: Debug Options
   :widths: 10 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``check``
     - Do extra run-time checking.  Useful for debugging, but can potentially slow calculations down
     - OFF
   * - ``debug``
     - Whether to enable displaying messages with the DEBUG() series of statements. Also writes messages to out.debug.<P> where P is the (physical) process rank. Still requires the "DEBUG" group to be enabled in ``Monitor`` (that is ``Monitor::is_active("DEBUG")`` must be true for any output)
     - OFF
   * - ``debug_field``
     -
     - OFF
   * - ``debug_field_face``
     -
     - OFF

   * - ``debug_verbose``
     - Print periodically all field values.  See ``src/Field/field_FieldBlock.cpp``
     - OFF
   * - ``trace``
     - Print out detailed messages with the TRACE() series of statements
     - OFF
   * - ``trace_charm``
     - Print out messages with the TRACE_CHARM() and TRACEPUP() series of statements
     - OFF
   * - ``verbose``
     - Trace main phases
     - OFF

Misc Options
^^^^^^^^^^^^

The following options don't really belong in any other category

.. list-table:: Misc Options
   :widths: 10 30 5
   :header-rows: 1

   * - Name
     - Description
     - Default
   * - ``USE_PRECOMPILED_HEADERS``
     - Precompile headers to try to reduce compile time
     - ON

.. _how_to_specify_the_configuration:

Specifying Configuration Options
--------------------------------

All variables can be set either on the command line by ``-D<variable>=<value>`` or
in a machine config, see below.
For example, a configure line may look like

..  code-block:: bash

  cmake -DCHARM_ROOT=$(pwd)/../../charm/build-gcc-mpi-proj -DEnzo-E_CONFIG=msu_hpcc_gcc -DGrackle_ROOT=${HOME}/src/grackle/build-gcc -Duse_projections=ON -Duse_jemalloc=ON -Dbalance=ON  ..

To see all available (and selected) options you can also run ``ccmake .`` in the
build directory (after running ``cmake`` in first place), or use the ``ccmake`` GUI 
directly to interactively configure Enzo-E by calling ``ccmake ..`` in an empty build
directory.

If packages (external libraries) are not found automatically or if the wrong one was
picked up, you can specify the search path by
``-D<package_name>_ROOT=/PATH/TO/PACKAGE/INSTALL``,
cf., the ``cmake`` example command just above.
Note:

* these package locations are also picked up from the environment, i.e., an alternative option  is ``export <package_name>_ROOT=/PATH/TO/PACKAGE/INSTALL`` .

* to specify the path to a ``libpng`` install, use ``-DPNG_ROOT=/PATH/TO/LIBPNG`` instead of ``-DLIBPNG_ROOT=...``.

The last option is a machine specific configuration file (see below).

In addition, the general `cmake` option is available to set basic optimization flags via ``CMAKE_BUILD_TYPE`` with values of:

* ``Release`` (typically ``-O3``)
* ``RelWithDebInfo`` (typically ``-O2 -g``)
* ``Debug`` (typically ``-O0 -g``)

The first choice is generally fastest, while the second is a sensible choice during development (the compiler performs most optimizations and includes debugging information in the executable)

Note that setting ``-DOPTIMIZE_FP=ON`` will modify each of these build types to include value-unsafe floating point optimizations.


Machine files
-------------

Finally, for convenience we provide the option to set default value for your own
machine/setup, see the ``*.cmake`` files in the ``config`` directory.

You can specify compilers and option in there that will be used a default when
``cmake`` is called with ``-DEnzo-E_CONFG=my_config_name`` where ``my_config_name`` requires
a corresponding ``config/my_config_name.cmake`` to exist.
The alternative directory for the machine configuration files is a ``.enzo-e`` directory
in your home directory.
If a file with the same name exists in your ``${HOME}/.enzo-e`` directory and in the ``config``
directory only the first one will be used.
*Note*, all command line parameter take precedence over the default options.
In other words, if ``USE_DOUBLE_PREC`` is ``ON`` in the machine file (or even automatically
through the global default), the running
``cmake -DEnzo-E_CONFIG=my_config_name -DUSE_DOUBLE_PREC=OFF ..`` will result in a single
precision version of Enzo-E.

Options in the machine file can also include the paths to external libraries and
can be set via a "cached string", i.e., via

..  code-block:: cmake

  set(CHARM_ROOT "/home/user/Charm/charm/build-mpi" CACHE STRING "my charm build")

Running
=======

In this section we run Enzo-E on a simple "Hello World" test program
and take a look at Enzo-E's output.

1. Run Enzo-E
-------------

An included "Hello World" problem can be run using the following
from the ``$CELLO_HOME`` directory:

     ``charmrun +p4 bin/enzo-e input/HelloWorld/Hi.in``

This assumes that the ``charmrun`` command is in your path.  If it
is not, then you will need to include the path name as well, e.g.:

     ``~/Charm/bin/charmrun +p4 bin/enzo-e input/HelloWorld/Hi.in``

This also assumes that local connections can be established passwordless.
If errors like

..  code-block:: bash

    Permission denied (publickey,gssapi-keyex,gssapi-with-mic,password,hostbased).
    Charmrun> Error 255 returned from remote shell (localhost:0)

are displayed, a node local run (i.e., no "remote" connections even to the local host) could be used instead by adding ``++local`` to ``charmrun``, e.g.:

     ``~/Charm/bin/charmrun ++local +p4 bin/enzo-e input/HelloWorld/Hi.in``

If you receive an error like

..  code-block:: bash

    Charmrun> Timeout waiting for node-program to connect

try running ``./bin/enzo-e`` without ``charmrun`` as crashes due to, e.g.,
libraries not being found may not be displaying.

If all goes well, Enzo-E will run the Hello World problem.  Below are
some of the generated images from the longer-running "HelloWorld.in"
problem (note "HelloWorld.in" runs for about an hour, compared to a
couple minutes for the shorter "Hi.in" input parameter file).  These
images show density and mesh hierarchy structure with blocks colored
by level and by age.

----

Time = 0.00

.. image:: hello-de-0000.png
   :scale: 40 %

.. image:: hello-mesh-level-0000.png
   :scale: 40 %

.. image:: hello-mesh-age-0000.png
   :scale: 40 %

----------------------

Time = 0.05

.. image:: hello-de-0086.png
   :scale: 40 %

.. image:: hello-mesh-level-0086.png
   :scale: 40 %

.. image:: hello-mesh-age-0086.png
   :scale: 40 %

----------------------

Time = 0.10

.. image:: hello-de-0165.png
   :scale: 40 %                   

.. image:: hello-mesh-level-0165.png
   :scale: 40 %

.. image:: hello-mesh-age-0165.png
   :scale: 40 %


If you look at the ``Hi.in`` parameter file contents, you will notice that there are some ``"include"`` directives that include other files.  When Enzo-E / Cello runs, it will generate a ``"parameters.out"`` file, which is the input file but with the included files inlined.  This ``"parameters.out"`` file is itself a valid Enzo-E / Cello parameter file (though you may wish to rename it before using it as a parameter file to avoid it being overwritten.)

If you encounter any problems in getting Enzo-E to compile or run,
please contact the Enzo-E / Cello community at the `users' mailing list <https://groups.google.com/g/enzo-e-users>`_ and someone will be happy to help resolve the problems.

----

2020-04-10: Updated with corrections from Joshua Smith.
