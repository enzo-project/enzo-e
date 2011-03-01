// $Id: test_DataBlock.cpp 1696 2010-08-04 05:56:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_DataBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar  1 20:42:57 UTC 2011
/// @brief    Test program for the DataBlock class

#include "test.hpp"

#include "mesh.hpp"

#include PARALLEL_CHARM_INCLUDE(test_DataBlock.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  unit_class ("DataBlock");

  DataDescr * data_descr = new DataDescr;
  DataBlock * data_block = new DataBlock (data_descr, 3,4,5);

  unit_func("DataBlock");
  unit_assert (data_descr != NULL);

  //----------------------------------------------------------------------
  unit_func("extent");

  data_block->set_extent(-1, 1, -2, 2, -3, 3);
  double lower[3];
  double upper[3];
  data_block->extent(&lower[0],&upper[0],
		     &lower[1],&upper[1],
		     &lower[2],&upper[2]);


  unit_assert(lower[0] == -1.0);
  unit_assert(upper[0] ==  1.0);
  unit_assert(lower[1] == -2.0);
  unit_assert(upper[1] ==  2.0);
  unit_assert(lower[2] == -3.0);
  unit_assert(upper[2] ==  3.0);

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_DataBlock.def.h)
