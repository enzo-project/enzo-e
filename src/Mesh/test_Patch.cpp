// $Id: test_Patch.cpp 1791 2010-10-29 22:54:06Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Patch class
 
#include "test.hpp"

#include "mesh.hpp"

#include PARALLEL_CHARM_INCLUDE(test_Patch.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();
  unit_class ("Patch");

  unit_func("Patch");

  DataDescr data_descr;

  //--------------------------------------------------
  // np = 1
  // ip = 0 [default]
  //--------------------------------------------------

  Patch * patch = new Patch;
  unit_assert(patch != NULL);

  unit_func("set_data_descr");

  patch->set_data_descr(&data_descr);
  unit_assert(patch->data_descr()==&data_descr);

  unit_func("set_size");

  // Patch = (7,3,2) cells

  patch->set_size(7,3,2);
  int npx,npy,npz;
  patch->size(&npx,&npy,&npz);
  unit_assert(npx==7 && npy==3 && npz==2);

  Layout * layout = patch->layout();

  unit_assert(layout != NULL);

  layout->set_process_range(0,1);
  layout->set_block_count(1,1,1);

  unit_func("set_extents");

  patch->set_extents(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

  double xm,xp,ym,yp,zm,zp;

  patch->extents(&xm,&xp,&ym,&yp,&zm,&zp);

  unit_assert(xm==0.0 && xp==1.0 &&
	      ym==0.0 && yp==1.0 &&
	      zm==0.0 && zp==1.0);
  
  unit_func("is_allocated");

  unit_assert(patch->blocks_allocated() == false);
  
  unit_func("allocate");

  patch->allocate_blocks();

  unit_assert(patch->blocks_allocated() == true);

  unit_func("block_count");

  unit_assert(patch->num_blocks()==1);
  
  unit_func("block");

  unit_assert(patch->block(0) != 0);

  // TEST BLOCK PROPERTIES
  unit_assert(false);

  // Deallocate local blocks

  unit_func("deallocate");

  patch->deallocate_blocks();

  unit_assert(patch->blocks_allocated() == false);

  // TEST MEMORY SIZE RESTORED--NO MEMORY LEAKS
  unit_assert(false)


  //--------------------------------------------------


  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Patch.def.h)
