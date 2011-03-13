// $Id: test_Patch.cpp 1791 2010-10-29 22:54:06Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Patch.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-06
/// @todo     Add 1D and 2D Patch tests
/// @todo     Add parallel Patch tests
/// @brief    Program implementing unit tests for the Patch class
 
#include "test.hpp"

#include "mesh.hpp"

#include PARALLEL_CHARM_INCLUDE(test_Patch.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  GroupProcess * group_process = GroupProcess::create();

  unit_init();
  unit_class ("Patch");

  //======================================================================
  // np = 1
  // ip = 0 [default]
  //======================================================================

  unit_func("Patch");

  FieldDescr * field_descr = new FieldDescr;

  // Set Patch size (12,12,12)

  int patch_size[] = {12,12,12};

  int patch_blocking[] = {3,3,3};

  Patch * patch = new Patch
    (NULL, group_process,
     patch_size[0],     patch_size[1],     patch_size[2],
     patch_blocking[0], patch_blocking[1], patch_blocking[2]);

  unit_assert(patch != NULL);

  //--------------------------------------------------
  unit_func("size");

  // Test that patch size is correct

  int npx,npy,npz;

  patch->size(&npx,&npy,&npz);

  unit_assert(patch_size[0]==npx && 
	      patch_size[1]==npy && 
	      patch_size[2]==npz);

  //--------------------------------------------------

  unit_func("blocking");

  // Test that patch blocking is correct

  int nbx,nby,nbz;

  patch->blocking(&nbx,&nby,&nbz);

  unit_assert(patch_blocking[0]==nbx && 
	      patch_blocking[1]==nby && 
	      patch_blocking[2]==nbz);

  //--------------------------------------------------

  unit_func("set_extents");

  // Set domain extents

  double domain_lower[] = {0.0, 0.0, 0.0};
  double domain_upper[] = {1.0, 1.0, 1.0};

  patch->set_extents(domain_lower[0], domain_upper[0], 
		     domain_lower[1], domain_upper[1], 
		     domain_lower[2], domain_upper[2]);

  // Test that the domain extents are correct

  double xm,xp,ym,yp,zm,zp;

  patch->extents(&xm,&xp,&ym,&yp,&zm,&zp);

  unit_assert(xm==domain_lower[0] && xp==domain_upper[0] &&
	      ym==domain_lower[1] && yp==domain_upper[1] &&
	      zm==domain_lower[2] && zp==domain_upper[2]);
  
  // Initialize how the Layout distributes the Patch data

  Layout * layout = patch->layout();

  unit_assert(layout != NULL);

  layout->set_process_range(0,1);

  // Test allocation of Patch into Blocks

  unit_func("is_allocated");

  unit_assert(patch->blocks_allocated() == false);
  
  unit_func("allocate");

  patch->allocate_blocks(field_descr);

  unit_assert(patch->blocks_allocated() == true);

  // Test that the allocated Blocks were initialized correctly

  unit_func("num_blocks");

  unit_assert(patch->num_blocks()==(size_t)nbx*nby*nbz);

  // loop over local data blocks and test their existence and properties

  ItBlock itBlocks (patch);

  Block *  block = 0;
  FieldBlock * field_block = 0;

  size_t count = 0;

  while ((block = ++itBlocks)) {

    ++count;

    unit_assert_quiet(block != NULL);

    if (block) {
      field_block = block->field_block();
      unit_assert_quiet(field_block != NULL);
      
    }

    // Test Block
    if (block) {
      // NO TESTS
    }

    // Test FieldBlock
    if (block && field_block) {

      // Test block size
      int nfx, nfy, nfz;

      unit_func("FieldBlock::size");
      field_block->size(&nfx,&nfy,&nfz);
      unit_assert_quiet (nfx == patch_size[0] / patch_blocking[0]);
      unit_assert_quiet (nfy == patch_size[1] / patch_blocking[1]);
      unit_assert_quiet (nfz == patch_size[2] / patch_blocking[2]);

      // Get block position in the Patch

      int ibx,iby,ibz;

      itBlocks.index(&ibx,&iby,&ibz);
      
      size_t ib = ibx + nbx*(iby + nby*ibz);
      unit_assert_quiet (count == ib);

      // Test block extents

      double xm,ym,zm;
      block->extent (&xm,&xp,&ym,&yp,&zm,&zp);

      
    }

    
    if (! block) {
      WARNING("test_Patch","Block tests skipped since Block not allocated");
    }

    if (block && ! field_block) {
      WARNING("test_Patch","Block tests skipped since FieldBlock not allocated");
    }

    // TEST BLOCK PROPERTIES
    //    unit_assert(unit_incomplete);

  }

  unit_func("num_blocks");
  unit_assert(count == patch->num_blocks());

  unit_assert(patch->block(0) != NULL);

  // Deallocate local blocks

  //--------------------------------------------------
  unit_func("deallocate");

  patch->deallocate_blocks();

  unit_assert(patch->blocks_allocated() == false);

  //--------------------------------------------------

  delete patch;
  delete field_descr;

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_Patch.def.h)
