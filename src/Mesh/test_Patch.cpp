// $Id: test_Patch.cpp 1791 2010-10-29 22:54:06Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @brief    Program implementing unit tests for the Patch class
 
#include <stdio.h>
#include <string>

#include "cello.hpp"
#include "test.hpp"
#include "mesh.hpp"

#include "parallel.def"

#include PARALLEL_CHARM_INCLUDE(test_Patch.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();
  unit_class ("Patch");

  unit_func("Patch");

  DataDescr data_descr;
  Patch * patch = new Patch;
  unit_assert(patch != NULL);

  unit_func("set_data_descr");

  patch->set_data_descr(&data_descr);
  unit_assert(patch->data_descr()==&data_descr);

  unit_func("set_size");

  patch->set_size(7,3,2);
  int npx,npy,npz;
  patch->size(&npx,&npy,&npz);
  unit_assert(npx==7 && npy==3 && npz==2);

  unit_func("set_layout");

  Layout * layout = new Layout;
  layout->set_process_range(0,1);
  layout->set_block_count(1,1,1);

  patch->set_layout (layout);
  unit_assert(patch->layout() == layout);

  unit_func("set_extents");

  patch->set_extents(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

  double xm,xp,ym,yp,zm,zp;

  patch->extents(&xm,&xp,&ym,&yp,&zm,&zp);

  unit_assert(xm==0.0 && xp==1.0 &&
	      ym==0.0 && yp==1.0 &&
	      zm==0.0 && zp==1.0);
  
//   void p_evolve();


  
// void allocate_blocks() ();

// /// Deallocate local blocks
// void deallocate_blocks() ();


  int block_count()  throw();
  

  DataBlock * block(int i) throw();

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Patch.def.h)
