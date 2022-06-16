.. _Frontera:

Getting Started: Frontera
=========================

..  code-block:: bash

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
  git checkout v7.0.0
  mkdir build-icc-mpi
  cd build-icc-mpi
  cmake -DNETWORK=mpi -DSMP=OFF -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort ..
  make

  # Build Enzo-E
  cd ~/src
  git clone https://github.com/enzo-project/enzo-e.git
  cd enzo-e
  mkdir build-icc-mpi
  cd build-icc-mpi
  cmake -DCHARM_ROOT=${HOME}/src/charm/build-icc-mpi -DGrackle_ROOT=${HOME}/src/grackle/build-icc -DEnzo-E_CONFIG=frontera_icc ..
  make
  # To run Enzo-E simply call `ibrun ./build-icc-mpi/bin/enzo-e input/HelloWorld/Hi.in +balancer GreedyLB` as usual
