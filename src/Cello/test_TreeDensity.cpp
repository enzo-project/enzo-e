// See LICENSE_CELLO file for license and copyright information

/// @file     test_TreeDensity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-23
/// @brief    Test program for the Tree class

#include "main.hpp"
#include "test.hpp"
#include "test_mesh.hpp"

#include "mesh.hpp"
#include "disk.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  //--------------------------------------------------
  // Read input HDF5 field into level array
  //--------------------------------------------------

  
  if (PARALLEL_ARGC != 4) {
    PARALLEL_PRINTF("Usage: %s <file_name> <field_name> <max_level>\n\n",PARALLEL_ARGV[0]);
    unit_assert(false);
    PARALLEL_EXIT;
  }

  const char * file_name  = PARALLEL_ARGV[1];
  const char * field_name = PARALLEL_ARGV[2];
  int max_level = atoi(PARALLEL_ARGV[3]);

  int nx,ny,nz;
  int * levels = create_levels_from_hdf5 (file_name, field_name,
					  &nx, &ny, &nz,  max_level);

			  
  // --------------------------------------------------
  // Write count of tagged zones in each level
  // --------------------------------------------------

  int * sum_field = new int [max_level+1];

  for (int i=0; i<=max_level; i++) sum_field[i]=0;

  for (int i=0; i<nx*ny*nz; i++) {
    sum_field[levels[i]]++;
  }

  int * zones_per_block = new int [max_level+1];
  // compute number of zones in a block at each level

  int n=nx*ny*nz;
  for (int i=0; i<=max_level; i++) {
    zones_per_block[i]=n;
    n/=8;
  }

  for (int i=0; i<=max_level; i++) TRACE2("Level %d  zones %d\n",i,sum_field[i]);

  // --------------------------------------------------
  // Create tree from level array
  // --------------------------------------------------

  Timer timer;
  timer.start();
  int d=3;
  int r=2;
  Tree tree (d,r);
  create_tree_from_levels (&tree, levels,nx,ny,nz);

  TRACE1 ("Initial time = %f",timer.value());
  TRACE1 ("Initial nodes = %d",tree.num_nodes());
  TRACE1 ("Initial depth = %d",tree.max_level());

  int * sum_mesh = new int [max_level+1];

  {
    for (int i=0; i<=tree.max_level(); i++) sum_mesh[i]=0;
    ItNode it_node(&tree);
    while (it_node.next_leaf()) {
      ++sum_mesh[it_node.node_trace()->level()];
    }
    for (int i=0; i<=max_level; i++) TRACE2("Initial level %d  blocks %d\n",i,sum_mesh[i]);
    for (int i=0; i<=max_level; i++) TRACE2("Initial level %d  zones %d\n",i,zones_per_block[i]*sum_mesh[i]);
  }

  // --------------------------------------------------
  // Write tree to file
  // --------------------------------------------------

  int mx=2048,my=2048;
  double th= 0.3*M_PI; // spin
  double ph= 0.1*M_PI;
  double ps= -0.06*M_PI;
  int falloff = 3;

  create_image_from_hdf5 (file_name,field_name,
			  "density_field.png",
			  mx,my, 0,max_level, ph,th,ps, 0.5, false, falloff);

  create_image_from_tree (&tree,"density_3d_1-initial.png",
			  mx,my, 0,max_level, ph,th,ps, 0.5, false, falloff);

  for (int i=0; i<=tree.max_level(); i++) {
    char filename[40];
    sprintf (filename,"density_3d_1-initial-L%d.png",i);
    create_image_from_tree (&tree,filename,
			  mx,my, i,i, ph,th,ps, 0.5, false, falloff);
  }

  double a90 = 0.5*M_PI;
  create_image_from_tree (&tree,"density_xy_1-initial.png",
			  mx,my, 0,max_level, 0.0,0.0,0.0, 1.0, true,0);
  create_image_from_tree (&tree,"density_yz_1-initial.png",
			  mx,my, 0,max_level, 0.0,-a90,a90, 1.0, true,0);
  create_image_from_tree (&tree,"density_zx_1-initial.png",
			  mx,my, 0,max_level, a90,a90,0.0, 1.0, true,0);

  // --------------------------------------------------
  // Balance tree
  // --------------------------------------------------

  timer.clear();
  timer.start();

  tree.balance();
  
  TRACE1 ("Balanced time = %f",timer.value());
  TRACE1 ("Balanced nodes = %d",tree.num_nodes());
  TRACE1 ("Balanced depth = %d",tree.max_level());

  {
    for (int i=0; i<=tree.max_level(); i++) sum_mesh[i]=0;
    ItNode it_node(&tree);
    while (it_node.next_leaf()) {
      ++sum_mesh[it_node.node_trace()->level()];
    }
    for (int i=0; i<=max_level; i++) TRACE2("Balanced level %d  blocks %d\n",i,sum_mesh[i]);
    for (int i=0; i<=max_level; i++) TRACE2("Balanced level %d  zones %d\n",i,zones_per_block[i]*sum_mesh[i]);
  }
  // --------------------------------------------------
  // Write tree to file
  // --------------------------------------------------

  create_image_from_tree (&tree,"density_3d_2-balanced.png",
			  mx,my,  0,max_level, ph,th,ps, 0.5, false, falloff);

  create_image_from_tree (&tree,"density_xy_2-balanced.png",
			  mx,my, 0,max_level, 0.0,0.0,0.0, 1.0, true,0);
  create_image_from_tree (&tree,"density_yz_2-balanced.png",
			  mx,my, 0,max_level, 0.0,-a90,a90, 1.0, true,0);
  create_image_from_tree (&tree,"density_zx_2-balanced.png",
			  mx,my, 0,max_level, a90,a90,0.0, 1.0, true,0);

  // --------------------------------------------------
  // Patch coalescing
  // --------------------------------------------------

  timer.clear();
  timer.start();


  tree.coalesce();

  TRACE1 ("Coalesced time = %f",timer.value());
  TRACE1 ("Coalesced nodes = %d",tree.num_nodes());
  TRACE1 ("Coalesced depth = %d",tree.max_level());

  {
    for (int i=0; i<=tree.max_level(); i++) sum_mesh[i]=0;
    ItNode it_node(&tree);
    while (it_node.next_leaf()) {
      ++sum_mesh[it_node.node_trace()->level()];
    }
    for (int i=0; i<=max_level; i++) TRACE2("Coalesced level %d  blocks %d\n",i,sum_mesh[i]);
    for (int i=0; i<=max_level; i++) TRACE2("Coalesced level %d  zones %d\n",i,zones_per_block[i]*sum_mesh[i]);
  }

  // --------------------------------------------------
  // Write tree to file
  // --------------------------------------------------

  create_image_from_tree (&tree,"density_3d_3-coalesced.png",
			  mx,my,  0,max_level, ph,th,ps, 0.5, false, falloff);

  create_image_from_tree (&tree,"density_xy_3-coalesced.png",
			  mx,my, 0,max_level, 0.0,0.0,0.0, 1.0, true,0);
  create_image_from_tree (&tree,"density_yz_3-coalesced.png",
			  mx,my, 0,max_level, 0.0,-a90,a90, 1.0, true,0);
  create_image_from_tree (&tree,"density_zx_3-coalesced.png",
			  mx,my, 0,max_level, a90,a90,0.0, 1.0, true,0);

  // --------------------------------------------------
  // --------------------------------------------------
  // Create 4-tree from level array
  // --------------------------------------------------

  {
    Timer timer;
    timer.start();
    int d=3;
    int r=4;
    Tree tree (d,r);
    create_tree_from_levels (&tree, levels,nx,ny,nz);

    TRACE1 ("Targeted time = %f",timer.value());
    TRACE1 ("Targeted nodes = %d",tree.num_nodes());
    TRACE1 ("Targeted depth = %d",tree.max_level());

    int * sum_mesh = new int [max_level+1];

    {
      for (int i=0; i<=tree.max_level(); i++) sum_mesh[i]=0;
      ItNode it_node(&tree);
      while (it_node.next_leaf()) {
	++sum_mesh[it_node.node_trace()->level()];
      }
      for (int i=0; i<=max_level; i++) TRACE2("Targeted level %d  blocks %d\n",i,sum_mesh[i]);
      for (int i=0; i<=max_level; i++) TRACE2("Targeted level %d  zones %d\n",i,zones_per_block[i]*sum_mesh[i]);
    }

    // --------------------------------------------------
    // Write tree to file
    // --------------------------------------------------

    create_image_from_tree (&tree,"density_3d_3-targeted.png",
			    mx,my,  0,max_level, ph,th,ps, 0.5, false, falloff);

    create_image_from_tree (&tree,"density_xy_3-targeted.png",
			    mx,my, 0,max_level, 0.0,0.0,0.0, 1.0, true,0);
    create_image_from_tree (&tree,"density_yz_3-targeted.png",
			    mx,my, 0,max_level, 0.0,-a90,a90, 1.0, true,0);
    create_image_from_tree (&tree,"density_zx_3-targeted.png",
			    mx,my, 0,max_level, a90,a90,0.0, 1.0, true,0);

    // --------------------------------------------------
  }

  unit_finalize();

  delete [] levels;
  delete [] sum_field;
  delete [] sum_mesh;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

