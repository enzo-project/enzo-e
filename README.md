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

This branch uses  `cmake` to build. Note that the project should not be
compiled from within the repository by in another build repository. For cloning
the repository, we suggest using
```
mkdir enzo-e-project
cd enzo-e-project
git clone https://github.com/forrestglines/enzo-e.git
cd enzo-e
git checkout cmake
cd ../
```
For building, we create a new directory inside of `enzo-e-project`, export a `CHARM_HOME` environment variable, then use cmake to build.
```
mkdir build
cd build
cmake  ../enzo-e
make
```
If `ninja` is installed, the `ninja` build system can be used for faster build times
```
cmake  -GNinja ../enzo-e
ninja
```

The Enzo-E executable is built within `src/`. Example inputs are copied into build directory as part of the `cmake` configuration. The executable can be run using
```
$CHARM_HOME/bin/charmrun +p4 src/enzo_e_exe input/test_cosmo-bcg.in
```

