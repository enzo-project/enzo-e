.. _Pleiades:

Getting Started: Pleiades
=========================

..  code-block:: bash

  # Load more up-to-date software stack
  module use -a /nasa/modulefiles/testing
  module purge
  module load pkgsrc/2021Q1 gcc/10.2 mpi-hpe/mpt.2.23 boost/1.76 comp-intel/2020.4.304 hdf5/1.12.0_serial pkgsrc/2021Q1
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
  git checkout v7.0.0
  mkdir build-icc-mpi
  cd build-icc-mpi
  cmake -DNETWORK=mpi -DSMP=OFF -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort ..
  make

  # Build Enzo-E
  cd ~/src
  git clone https://github.com/enzo-project/enzo-e.git
  cd enzo-e
  # Custom environment override for the cmake call specific to Pleiades system as
  # during the linking step (done with the `charmc` wrapper) the `mpicxx` wrapper is called,
  # which, in turn, by default calls `g++` (but we need to link with `icpc` at the lowest level).
  export MPICXX_CXX=icpc
  mkdir build-icc-mpi
  cd build-icc-mpi
  cmake -DCHARM_ROOT=${HOME}/src/charm/build-icc-mpi -DGrackle_ROOT=${HOME}/src/grackle/build-icc -DEnzo-E_CONFIG=pleiades_icc ..
  make
  # To run Enzo-E simply call `mpiexec ./build-icc-mpi/bin/enzo-e input/HelloWorld/Hi.in` as usual
