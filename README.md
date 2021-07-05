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

