// See LICENSE_CELLO file for license and copyright information

/// @file     test_Block.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar  1 20:42:57 UTC 2011
/// @brief    Test program for the Block class

#include "main.hpp" 
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Block");

#ifdef CONFIG_USE_CHARM
  CProxy_Patch proxy_patch;
#endif
  Factory factory;
  int patch_id = 0;
  int patch_rank = 0;
  
  Block * block = factory.create_block
    (0,0,0, 
     1,1,1,
     3,4,5,
     -1.0,-2.0,-3.0,
     2.0,  4.0, 6.0,
     patch_id,
     patch_rank,
#ifdef CONFIG_USE_CHARM
     proxy_patch,
#endif
     1);

  unit_func("Block");
  unit_assert (block != NULL);

  //----------------------------------------------------------------------

  double lower[3];
  double upper[3];

  block->lower(&lower[0],&lower[1],&lower[2]);
  block->upper(&upper[0],&upper[1],&upper[2]);

  unit_func("lower");
  unit_assert(lower[0] == -1.0);
  unit_assert(lower[1] == -2.0);
  unit_assert(lower[2] == -3.0);

  unit_func("upper");
  unit_assert(upper[0] ==  1.0);
  unit_assert(upper[1] ==  2.0);
  unit_assert(upper[2] ==  3.0);

  //----------------------------------------------------------------------

  unit_func("write");


  // Block
  //   index
  //   size
  //   lower
  //   upper
  //   cycle
  //   time
  //   dt

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END
