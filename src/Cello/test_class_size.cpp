// See LICENSE_CELLO file for license and copyright information

/// @file     test_class_size.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Class_Size class

#include "cello.hpp"
#include "main.hpp"
#include "test.hpp"

#include "enzo.hpp"

//#include "simulation.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  printf ("%4ld sizeof(Block)\n",sizeof(Block));
  printf ("%4ld sizeof(Boundary)\n",sizeof(Boundary));
  printf ("%4ld sizeof(EnzoBlock)\n",sizeof(EnzoBlock));
  printf ("%4ld sizeof(Factory)\n",sizeof(Factory));
  printf ("%4ld sizeof(FieldData)\n",sizeof(FieldData));
  printf ("%4ld sizeof(FieldDescr)\n",sizeof(FieldDescr));
  printf ("%4ld sizeof(FieldFace)\n",sizeof(FieldFace));
  // the precise size of FileHdf5 is hidden (otherwise, lots of unrelated
  // components would have a transitive dependence of hdf5)
  //printf ("%4ld sizeof(FileHdf5)\n",sizeof(FileHdf5));
  printf ("%4ld sizeof(Hierarchy)\n",sizeof(Hierarchy));
  printf ("%4ld sizeof(Method)\n",sizeof(Method));
  printf ("%4ld sizeof(Parameters)\n",sizeof(Parameters));
  printf ("%4ld sizeof(Problem)\n",sizeof(Problem));
  printf ("%4ld sizeof(Simulation)\n",sizeof(Simulation));
  printf ("%4ld sizeof(Stopping)\n",sizeof(Stopping));

  //--------------------------------------------------

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

