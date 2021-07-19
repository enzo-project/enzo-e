# Enzo-E/Cello
[![CircleCI](https://circleci.com/gh/enzo-project/enzo-e.svg?style=svg)](https://circleci.com/gh/enzo-project/enzo-e)
[![Documentation Status](https://readthedocs.org/projects/enzo-e/badge/?version=latest)](https://enzo-e.readthedocs.io/en/latest/?badge=latest)
 
## What is this repository for?

   The purpose of the Enzo-E/Cello project is to develop an
   astrophysics and cosmology software application "Enzo-E", which is
   built on a highly scalable parallel adaptive mesh refinement (AMR)
   software framework "Cello," which is being developed
   concurrently. The Enzo-E application will be capable of running
   extreme scale numerical simulations to address current scientific
   questions in astrophysics and cosmology, including star formation,
   molecular cloud turbulence, interstellar medium, galaxy formation,
   interglactic medium, formation of the first stars and galaxies,
   galaxy clusters, and cosmic reionization.  The Cello AMR framework
   can be used independently of Enzo-E, thus enabling researchers in
   other scientific fields to develop AMR applications capable of
   running on "Petascale-and-beyond" HPC platforms.

   Enzo-E/Cello is currently under active development.  It still has
   several known issues that need to be addressed before it can be
   used for scientific research, though it is still usable as a
   hydrodynamics solver, especially in an educational setting.

## Resources

   * [Documentation](https://enzo-e.readthedocs.io/)
   * [Users Mailing List](https://groups.google.com/forum/#!forum/enzo-e-users)
   * [Developers Mailing List](https://groups.google.com/forum/#!forum/enzo-e-dev)   

The Enzo-E project is fully open-source and contributions are very
welcome. If you're interested in getting involved, come say hello on
the developers mailing list!

## Branch Notes:

This branch uses  `cmake` to build. In source builds (e.g., calling `make` from within
the folder containing the sources it not supported and a separate folder (either
within the repository or anywhere else on the filesystem) must be used.

First, we install the latest `Charm++`, which now also uses `cmake` as default
build system (as example we use the pure MPI backend).

```bash
git clone https://github.com/UIUC-PPL/charm.git
cd charm
git checkout v7.0.0-rc1
# Note, the directory name can be anything
mkdir build-mpi
cd build-mpi
cmake -DNETWORK=mpi -DSMP=OFF ..
make
```

Second, we install `Enzo-E` (in a different directory)

```bash
git clone https://github.com/forrestglines/enzo-e.git
cd enzo-e
git checkout cmake
# Again, the directory name can be anything
mkdir build-mpi
cd build-mpi
```

For a basic configuration on a linux system with GNU compiler stack the following
command should work out of the box. If not please report.
Also see following subsection for more configuration options.

```
cmake -DCHARM_ROOT=/PATH/TO/charm/build-mpi -DEnzo-E_CONFIG=linux_gcc -DUSE_GRACKLE=OFF ..
make
```

Note, if `ninja` is installed, the `ninja` build system can be used for faster build times by adding `-GNinja` (before the `..`) to the `cmake` command and calling `ninja` afterwards instead of `make.

The Enzo-E executable is built within `bin/`. 

#### Configuration options

Current `cmake` options follow the notation of the SCons options (potentially 
subject to change), i.e., the following (default value at the end of the line)
is available:
- `USE_GRACKLE` "Use Grackle Chemistry" ON
- `USE_DOUBLE_PREC` "Use double precision. Turn off for single precision." ON
- `new_output` "Temporary setting for using new Output implementation" OFF
- `node_size` "Maximum number of procesess per shared-memory node (can be larger than needed)" 64
- `trace` "Print out detailed messages with the TRACE() series of statements" OFF
- `verbose` "Trace main phases" OFF
- `trace_charm` "Print out messages with the TRACE_CHARM() and TRACEPUP() series of statements" OFF
- `debug` "Whether to enable displaying messages with the DEBUG() series of
statements. Also writes messages to out.debug.<P> where P is the
(physical) process rank. Still requires the \"DEBUG\" group to be
enabled in Monitor (that is Monitor::is_active(\"DEBUG\") must be true for any output) OFF
- `debug_field` "" OFF
- `debug_field_face` "" OFF
- `check` "Do extra run-time checking.  Useful for debugging, but can potentially slow     calculations down" OFF
- `debug_verbose` "Print periodically all field values.  See src/Field/field_FieldBlock.cpp" OFF
- `memory` "Track dynamic memory statistics.  Can be useful, but can cause problems on some systems that also override new [] () / delete [] ()" OFF
- `balance` "Enable charm++ dynamic load balancing" ON`
- `balancer` "Charm++ load balancer to use" "CommonLBs"
- `use_gprof` "Compile with -pg to use gprof for performance profiling" OFF
- `use_performance` "Use Cello Performance class for collecting performance
data (currently requires global reductions, and may not be fully
functional) (basic time data on root processor is still output)" ON
- `use_projections` "Compile the CHARM++ version for use with the Projections performance tool." OFF
-` use_jemalloc "Use the jemalloc library for memory allocation" OFF
- `smp` "Use Charm++ in SMP mode." OFF
- `use_papi` "Use the PAPI performance API" OFF

All variables can be set either on the commad line by `-D<variable>=<value>` or
in a machine config, see below.
For example, a configure line may look like
```bash
cmake -DCHARM_ROOT=$(pwd)/../../charm/build-gcc-mpi-proj -DEnzo-E_CONFIG=msu_hpcc_gcc -DGrackle_ROOT=${HOME}/src/grackle/build-gcc -Duse_projections=ON -Duse_jemalloc=ON -Dbalance=ON  ..
```
To see all available (and selected) options you can also run `ccmake .` in the
build directory (after running `cmake` in first place), or use the `ccmake` GUI 
directly to interactively configure Enzo-E by calling `ccmake ..` in an empty build
directory.

If packages (external libraries) are not found automatically or if the wrong one was
picked up, you can specify the search path by
`-D<package_name>_ROOT=/PATH/TO/PACKAGE/INSTALL`,
cf., the `cmake` example command just above.
Note, these package location are also picked up from the environment, i.e., an alternative
option  is `export <package_name>_ROOT=/PATH/TO/PACKAGE/INSTALL` .

The last option is a machine specific configuration file (see below).

In addition, the general `cmake` option to set basic optimization flags via
`CMAKE_BUILD_TYPE` with values of 
- `Release` (typically `-O3`),
- `RelWithDebInfo` (typically `-O2 -g`), and
- `Debug` (typically `-O0 -g`)

are available.


##### Machine files
Finally, for convenience we provide the option to set default value for your own
machine/setup, see the `*.cmake` files in the `config` directory.

You can specify compilers and option in there that will be used a default when
`cmake` is called with `-DEnzo-E_CONFG=my_config_name` where `my_config_name` requires
a corresponding `config/my_config_name.cmake` to exist.
The alternative directory for the machine configuration files is a `.enzo-e` directory
in your home directory.
If a file with the same name exists in your `${HOME}/.enzo-e` directory and in the `config`
directory only the first one will be used.
*Note*, all command line parameter take precedence over the default options.
In other words, if `USE_DOUBLE_PREC` is `ON` in the machine file (or even automatically
through the global default), the running
`cmake -DEnzo-E_CONFIG=my_config_name -DUSE_DOUBLE_PREC=OFF ..` will result in a single
precision version of Enzo-E.

Options in the machine file can also include the paths to external libraries and
can be set via a "cached string", i.e., via

```cmake
set(CHARM_ROOT "/home/user/Charm/charm/build-mpi" CACHE STRING "my charm build")
```

### Pleiades with Intel compilers

```bash
# Load more up-to-date software stack
module use -a /nasa/modulefiles/testing
module purge
module load pkgsrc/2021Q1 gcc/10.2 mpi-hpe/mpt.2.23 boost/1.75 comp-intel/2020.4.304 hdf5/1.12.0_serial pkgsrc/2021Q1

# Build Grackle (optional)
# Following https://grackle.readthedocs.io/en/latest/Installation.html
mkdir -p ~/src
cd ~/src
git clone https://github.com/grackle-project/grackle
cd grackle
git submodule update --init
./configure
# create build directory as install target later
mkdir build-icc
cd src/clib
# Adjust config to loaded modules, set install path, and optimization options
sed -i 's/1.8.18/1.12.0/;s/2020.2.254/2020.4.304/;s#$(HOME)/local#$(PWD)/../../build-icc#;s/-axAVX -xSSE4.1 -ip -ipo/-axCORE-AVX512,CORE-AVX2/' Make.mach.nasa-pleiades
make machine-nasa-pleiades
make
make install

# Build Charm++
cd ~/src
git clone https://github.com/UIUC-PPL/charm.git
cd charm
git checkout v7.0.0-rc1
mkdir build-icc-mpi
cd build-icc-mpi
cmake -DNETWORK=mpi -DSMP=OFF -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort ..
make

# Build Enzo-E
cd ~/src
git clone https://github.com/forrestglines/enzo-e.git
cd enzo-e
git checkout cmake

# Custom environment override for the cmake call specific to Pleiades system as
# during the linking step (done with the `charmc` wrapper) the `mpicxx` wrapper is called,
# which, in turn, by default calls `g++` (but we need to link with `icpc` at the lowest level).
export MPICXX_CXX=icpc

mkdir build-icc-mpi
cd build-icc-mpi
cmake -DCHARM_ROOT=${HOME}/src/charm/build-icc-mpi -DGrackle_ROOT=${HOME}/src/grackle/build-icc -DEnzo-E_CONFIG=pleiades_icc ..
make

# To run Enzo-E simply call `mpiexec ./build-icc-mpi/bin/enzo-e input/HelloWorld/Hi.in` as usual

```

### TACC Frontera with Intel compilers

```bash
# Assuming the default Intel stack is loaded only boost and hdf5 modules are required
module load boost hdf5

# Build Grackle (optional)
# Following https://grackle.readthedocs.io/en/latest/Installation.html
mkdir -p ~/src
cd ~/src
git clone https://github.com/grackle-project/grackle
cd grackle
git submodule update --init
./configure
# create build directory as install target later
mkdir build-icc
cd src/clib
# Adjust config to  set install path, and optimization options
# TACC Stampede config works just fine for TACC Frontera
sed -i 's#$(HOME)/local#$(PWD)/../../build-icc#;s/-xCORE-AVX2/-xCORE-AVX512/' Make.mach.tacc-stampede-intel
make machine-tacc-stampede-intel
make
make install

# Build Charm++
cd ~/src
git clone https://github.com/UIUC-PPL/charm.git
cd charm
git checkout v7.0.0-rc1
mkdir build-icc-mpi
cd build-icc-mpi
cmake -DNETWORK=mpi -DSMP=OFF -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort ..
make

# Build Enzo-E
cd ~/src
git clone https://github.com/forrestglines/enzo-e.git
cd enzo-e
git checkout cmake
mkdir build-icc-mpi
cd build-icc-mpi
cmake -DCHARM_ROOT=${HOME}/src/charm/build-icc-mpi -DGrackle_ROOT=${HOME}/src/grackle/build-icc -DEnzo-E_CONFIG=frontera_icc ..
make

# To run Enzo-E simply call `ibrun ./build-icc-mpi/bin/enzo-e input/HelloWorld/Hi.in +balancer GreedyLB` as usual

```

#### FAQ

##### Q: `/usr/bin/ld: libenzo.a(enzo_EnzoInitialInclinedWave.cpp.o): undefined reference to symbol '__libm_sse2_sincos'`

Double check that you're linking with Intel compiler, e.g., whether `export MPICXX_CXX=icpc` is set.

##### Q: `ld: ./src/grackle/build-icc/lib/libgrackle.so: undefined reference to '__must_be_linked_with_icc_or_xild'`

Make sure `-ip and -ipo` flags between grackle and Cello/Enzo-E build are consistent.

##### Q: `enzo_EnzoMethodGrackle.cpp:(.text+0x1f5): undefined reference to `vtable for EnzoMethodGrackle'` or similar

Potentially, there are preprocessor defines missing when processing `*.ci` files so that there is a mismatch between header declaration and definitions.
Double check that all `#ifdefs` in the `src/*/*.ci` are also set in the `CHARM_PREPROC_DEFS` `CMake` variable.
