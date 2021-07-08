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
# Note, the Fortran flags have just been tested for gfortran so far
cmake -DCHARM_ROOT=/PATH/TO/charm/build-mpi -DCMAKE_CXX_COMPILER=/PATH/TO/charm/build-mpi/bin/charmc -DCMAKE_Fortran_FLAGS="-fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132" ..
make
```

Note, if `ninja` is installed, the `ninja` build system can be used for faster build times by adding `-GNinja` (before the `..`) to the `cmake` command and calling `ninja` afterwards instead of `make.

The Enzo-E executable is built within `bin/`. Example inputs are copied into build directory as part of the `cmake` configuration (TODO). The executable can be run using
```
/PATH/TO/charm/build-mpi/bin/charmrun +p4 src/enzo_e_exe input/test_cosmo-bcg.in
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

#### FAQ

##### Q: `/usr/bin/ld: libenzo.a(enzo_EnzoInitialInclinedWave.cpp.o): undefined reference to symbol '__libm_sse2_sincos'`

Double check that you're linking with Intel compiler, e.g., whether `export MPICXX_CXX=icpc` is set.

##### Q: `ld: ./src/grackle/build-icc/lib/libgrackle.so: undefined reference to '__must_be_linked_with_icc_or_xild'`

Make sure `-ip and -ipo` flags between grackle and Cello/Enzo-E build are consistent.

##### Q: `enzo_EnzoMethodGrackle.cpp:(.text+0x1f5): undefined reference to `vtable for EnzoMethodGrackle'` or similar

Potentially, there are preprocessor defines missing when processing `*.ci` files so that there is a mismatch between header declaration and definitions.
Double check that all `#ifdefs` in the `src/*/*.ci` are also set in the `CHARM_PREPROC_DEFS` `CMake` variable.
