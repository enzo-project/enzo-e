.. toctree::
   :maxdepth: 1

**************************
Frequently Asked Questions
**************************


Building
********

Q: Charm++ build fails with ``error: C compiler cannot create executables``
---------------------------------------------------------------------------

If the Charm++ build command fails during the configure step with the error above
*and* the mentioned `config.log` (e.g., in the `mpi-linux-x86_64-gfortran-smp-gcc/tmp`
directory) contains an error like

.. code-block:: bash

  conftest.c:133:1: error: expected identifier or '(' before numeric constant
   4.8
   ^~~

then the error could be caused by the Python version of the build environment.
The error has been observed for Python 3.7.0 on NASA's Pleiades supercomputer
and could be circumvented by using the Python 3.5.2 module.

Q: Enzo-E build fails with ``Function ... has no IMPLICIT type``
----------------------------------------------------------------

If the Enzo-E build fails with error messages like

.. code-block:: bash

    build/Enzo/flux_hll.F:351:50:

        write (6,*) 'e', char2(i-1), char2(i), char2(i+1)
    Error: Function 'char2' at (1) has no IMPLICIT type
   
and you are using a newer ``gfortran`` (problem confirmed to show for
at least major releases 8 and 9) try adding

.. code-block:: bash

  args_arch_fortran = '-ffixed-line-length-132'

to your configuration file (cf., ``config/linux_gcc_9.py``).


Running
*******

Q: Starting Enzo-E fails with ``... error while loading shared libraries ...``
------------------------------------------------------------------------------

If a custom library (e.g., one not provided by the system or the module environment)
was used to compile Enzo-E, the library must also be available at the time of execution.

A common error seen is

.. code-block:: bash

  ./bin/enzo-p: error while loading shared libraries: libhdf5.so.103: cannot open shared object file: No such file or directory

To fix this add the path to the ``lib`` directory of the HDF5 library used when building Enzo, i.e.,
the on specified in the config file of the chosen architecture, cf., ``config/linux_gcc_9.py``.

For example

.. code-block:: bash

  export LD_LIBRARY_PATH=${HOME}/src/hdf5-1.10.7/install/lib:$LD_LIBRARY_PATH
