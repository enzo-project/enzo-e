// See LICENSE_CELLO file for license and copyright information

/// @file     test_Class_Size.cpp
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

  printf ("%4ld sizeof(CommBlock)\n",sizeof(CommBlock));
  printf ("%4ld sizeof(Boundary)\n",sizeof(Boundary));
  printf ("%4ld sizeof(EnzoCommBlock)\n",sizeof(EnzoCommBlock));
  printf ("%4ld sizeof(Factory)\n",sizeof(Factory));
  printf ("%4ld sizeof(FieldBlock)\n",sizeof(FieldBlock));
  printf ("%4ld sizeof(FieldDescr)\n",sizeof(FieldDescr));
  printf ("%4ld sizeof(FieldFace)\n",sizeof(FieldFace));
  printf ("%4ld sizeof(FileHdf5)\n",sizeof(FileHdf5));
  printf ("%4ld sizeof(GroupProcess)\n",sizeof(GroupProcess));
  printf ("%4ld sizeof(Hierarchy)\n",sizeof(Hierarchy));
  printf ("%4ld sizeof(Layout)\n",sizeof(Layout));
  printf ("%4ld sizeof(Method)\n",sizeof(Method));
  printf ("%4ld sizeof(Node)\n",sizeof(Node));
  printf ("%4ld sizeof(Parameters)\n",sizeof(Parameters));
  printf ("%4ld sizeof(Patch)\n",sizeof(Patch));
  printf ("%4ld sizeof(Problem)\n",sizeof(Problem));
  printf ("%4ld sizeof(Simulation)\n",sizeof(Simulation));
  printf ("%4ld sizeof(Stopping)\n",sizeof(Stopping));
  printf ("%4ld sizeof(Timestep)\n",sizeof(Timestep));
  printf ("%4ld sizeof(Tree)\n",sizeof(Tree));

  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

