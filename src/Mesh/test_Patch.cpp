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

  //======================================================================
  // np = 1
  // ip = 0 [default]
  //======================================================================

  unit_func("Patch","Patch");

  FieldDescr * field_descr = new FieldDescr;

  // Set Patch size (12,12,12)

  int patch_size[] = {12,12,12};

  int patch_blocking[] = {3,3,3};

  Patch * patch = new Patch
    (new Factory, group_process,
     patch_size[0],     patch_size[1],     patch_size[2],
     patch_blocking[0], patch_blocking[1], patch_blocking[2]);

  unit_assert(patch != NULL);

  //--------------------------------------------------
  unit_func("Patch","size");

  // Test that patch size is correct

  int npx,npy,npz;

  patch->size(&npx,&npy,&npz);

  unit_assert(patch_size[0]==npx && 
	      patch_size[1]==npy && 
	      patch_size[2]==npz);

  //--------------------------------------------------

  unit_func("Patch","blocking");

  // Test that patch blocking is correct

  int nbx,nby,nbz;

  patch->blocking(&nbx,&nby,&nbz);

  unit_assert(patch_blocking[0]==nbx && 
	      patch_blocking[1]==nby && 
	      patch_blocking[2]==nbz);

  //--------------------------------------------------


  // Set domain extents

  double domain_lower[] = {0.0, 0.0, 0.0};
  double domain_upper[] = {1.0, 1.0, 1.0};

  patch->set_lower(domain_lower[0], domain_lower[1], domain_lower[2]);
  patch->set_upper(domain_upper[0], domain_upper[1], domain_upper[2]);

  // Test that the domain extents are correct

  unit_func("Patch","set_lower");

  double xm,ym,zm;
  patch->lower(&xm,&ym,&zm);

  unit_assert(xm==domain_lower[0]);
  unit_assert(ym==domain_lower[1]);
  unit_assert(zm==domain_lower[2]);

  unit_func("Patch","set_upper");

  double xp,yp,zp;
  patch->upper(&xp,&yp,&zp);

  unit_assert(xp==domain_upper[0]);
  unit_assert(yp==domain_upper[1]);
  unit_assert(zp==domain_upper[2]);

  // Initialize how the Layout distributes the Patch data

  Layout * layout = patch->layout();

  unit_func("Patch","layout");
  unit_assert(layout != NULL);

  layout->set_process_range(0,1);

  // Test allocation of Patch into Blocks

  unit_func("Patch","blocks_allocated");

  unit_assert(patch->blocks_allocated() == false);
  
  unit_func("Patch","allocate_blocks");

  patch->allocate_blocks(field_descr);

  unit_assert(patch->blocks_allocated() == true);

  // Test that the allocated Blocks were initialized correctly

  unit_func("Patch","num_blocks");

  unit_assert(patch->num_blocks()==(size_t)nbx*nby*nbz);

  // loop over local data blocks and test their existence and properties

  ItBlock itBlocks (patch);

  Block *  block = 0;
  FieldBlock * field_block = 0;

  size_t block_counter = 0;

  while ((block = ++itBlocks)) {

    unit_func("Patch","allocate_blocks");
    unit_assert_quiet(block != NULL);

    unit_func("Block","field_block");
    field_block = block ? block->field_block() : NULL;
    unit_assert_quiet(field_block != NULL);

    // Test Block
    if (block) {
      // NO TESTS
    }

    // Test FieldBlock
    if (block && field_block) {

      // Test block size
      int nfx, nfy, nfz;

      unit_func("FieldBlock","size");
      field_block->size(&nfx,&nfy,&nfz);
      unit_assert_quiet (nfx == patch_size[0] / patch_blocking[0]);
      unit_assert_quiet (nfy == patch_size[1] / patch_blocking[1]);
      unit_assert_quiet (nfz == patch_size[2] / patch_blocking[2]);

      // Get block position in the Patch

      int ibx,iby,ibz;

      GroupProcess * group = patch->group();
      Layout      * layout = patch->layout();

      int ip = group->rank();

      int index_local = block_counter;
      int index_global = layout->global_index(ip,index_local);
      layout->block_indices(index_global,&ibx,&iby,&ibz);
      
      //      size_t ib = ibx + nbx*(iby + nby*ibz);

      unit_func("Layout","block_indices");
      unit_assert_quiet (unit_incomplete);

      // Test block extents

      double xmb,ymb,zmb;
      double xpb,ypb,zpb;

      block->lower (&xmb,&ymb,&zmb);
      block->upper (&xpb,&ypb,&zpb);

      // Not very rigorous
      unit_assert (xmb < xpb);
      unit_assert (ymb < ypb);
      unit_assert (zmb < zpb);

      // More rigorous tests below, but require global indices
      unit_func("Block","lower");
      unit_assert (unit_incomplete);

      unit_func("Block","upper");
      unit_assert (unit_incomplete);

      // Need ib? which is not implemented yet; note comparing floating point

      // unit_assert(xm + ibx*(xpb-xmb) == xmb);
      // unit_assert(ym + iby*(ypb-ymb) == ymb);
      // unit_assert(zm + ibz*(zpb-zmb) == zmb);

      // unit_assert(xm + (ibx+1)*(xpb-xmb) == xpb);
      // unit_assert(ym + (iby+1)*(ypb-ymb) == ypb);
      // unit_assert(zm + (ibz+1)*(zpb-zmb) == zpb);

    }

    unit_func("Block","index_patch");
    for (int ibz=0; ibz<nbz; ibz++) {
      for (int iby=0; iby<nby; iby++) {
	for (int ibx=0; ibx<nbx; ibx++) {

	  int jbx,jby,jbz;
	  block->index_patch(&jbx,&jby,&jbz);
	  printf ("%d %d %d\n",jbx,jby,jbz);
	  unit_assert_quiet(ibx==jbx);
	  unit_assert_quiet(iby==jby);
	  unit_assert_quiet(ibz==jbz);
	}
      }
    }
      
    unit_func("Block","neighbor");
    // Block::neighbor()
    unit_assert(false);

    // TEST BLOCK PROPERTIES
    //    unit_assert(unit_incomplete);

    ++block_counter;

  }

  unit_func("Patch","num_blocks");
  unit_assert(block_counter == patch->num_blocks());

  unit_func("Patch","allocate_blocks");

  // Deallocate local blocks

  //--------------------------------------------------
  unit_func("Patch","deallocate_blocks");

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
