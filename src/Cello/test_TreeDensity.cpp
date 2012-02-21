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

/* #define OUTPUT_DENSITY 0 */

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  //--------------------------------------------------
  // Read input HDF5 field into level array
  //--------------------------------------------------

  
  if (!(PARALLEL_ARGC == 5 || PARALLEL_ARGC == 6)) {
    PARALLEL_PRINTF("Usage: %s <file_name> [ <group_name> ] <field_name> <min_level> <max_level>\n\n",PARALLEL_ARGV[0]);
    unit_assert(false);
    PARALLEL_EXIT;
  }

  char * file_name = 0;
  char * group_name = 0;
  char * field_name = 0;
  int min_level = 0;
  int max_level = 0;

  int arg=1;
  file_name = PARALLEL_ARGV[arg++];
  if (PARALLEL_ARGC == 6) {
    group_name = PARALLEL_ARGV[arg++];
  } else {
    group_name = strdup("");
  }
  field_name = PARALLEL_ARGV[arg++];
  min_level = atoi(PARALLEL_ARGV[arg++]);
  max_level = atoi(PARALLEL_ARGV[arg++]);

  int nx,ny,nz;
  double tol = 1e-4;
  refine_type refine = refine_log;
  // refine_type refine = refine_slope;
  int * levels = hdf5_to_levels
    (file_name, group_name, field_name, &nx, &ny, &nz, min_level,max_level,
     refine,tol);

  int dimension = nz==1?2:3;
  int refinement=2;

  // --------------------------------------------------
  // // Write count of tagged zones in each level
  // // --------------------------------------------------

  // int * sum_field = new int [max_level+1];

  // for (int i=0; i<=max_level; i++) sum_field[i]=0;

  // int n=nx*ny*nz;
  // for (int i=0; i<n; i++) {
  //   sum_field[levels[i]]++;
  // }

  // int * zones_per_block = new int [max_level+1];
  // // compute number of zones in a block at each level

  // int r2d = (dimension == 2) ? 4 : 8;
  // for (int i=0; i<=max_level; i++) {
  //   zones_per_block[i]=n;
  //   n/=r2d;
  // }

  // --------------------------------------------------
  // Create tree from level array
  // --------------------------------------------------

  Timer timer;
  timer.start();
  Tree tree (dimension,refinement);

  levels_to_tree (&tree, levels,nx,ny,nz);

  timer.stop();

  printf ("tree initial time      = %f\n",timer.value());
  printf ("tree initial num_nodes = %d\n",tree.num_nodes());
  printf ("tree initial max_level = %d\n",tree.max_level());


  int * sum_mesh = new int [max_level+1];

  {
    for (int i=0; i<=tree.max_level(); i++) sum_mesh[i]=0;
    ItNode it_node(&tree);
    timer.clear();
    timer.start();
    while (it_node.next_leaf()) {
      ++sum_mesh[it_node.node_trace()->level()];
    }
    timer.stop();
    printf ("traversed %d nodes in %f Mnodes/s\n",tree.num_nodes(),
	    1e-6*tree.num_nodes()/timer.value());
  }

  // --------------------------------------------------
  // Write tree to file
  // --------------------------------------------------

  //  int mx=1024,my=1024;
  int mx=512,my=512;
  double th= 0.1*M_PI;
  double ph= 0.3*M_PI; // spin
  double ps= -0.06*M_PI;
  int falloff = 3;
  double scale = 1.0;
  double a90 = 0.5*M_PI;

  if (dimension == 2) {
    timer.clear();  timer.start();
    tree_to_png (&tree,"density_2d_1-initial.png",
		 mx,my, min_level,max_level, 0.0,0.0,0.0, scale, true, falloff);
    TRACE1 ("Time = %f",timer.value());
  } else  {
    timer.clear();  timer.start();
    tree_to_png (&tree,"density_3d_1-initial.png",
		 mx,my, min_level,max_level, th,ph,ps, scale, false, falloff);
    TRACE1 ("Time = %f",timer.value());
  }

#ifdef OUTPUT_DENSITY
    timer.clear();  timer.start();
  hdf5_to_png (file_name,group_name,field_name,"density_xy_field.png",
	       mx,my, min_level,max_level, 0.0,0.0,0.0, 1.0, true,falloff);
    TRACE1 ("Time = %f",timer.value());

  if (dimension == 3) {
    timer.clear();  timer.start();
    hdf5_to_png (file_name,group_name,field_name,"density_field.png",
		 mx,my, min_level,max_level, th,ph,ps, scale, false, falloff);
    TRACE1 ("Time = %f",timer.value());

    timer.clear();  timer.start();
    hdf5_to_png (file_name,group_name,field_name,"density_yz_field.png",
		 mx,my, min_level,max_level, 0.0,-a90,a90, 1.0, true,falloff);
    TRACE1 ("Time = %f",timer.value());
    timer.clear();  timer.start();
    hdf5_to_png (file_name,group_name,field_name,"density_zx_field.png",
		 mx,my, min_level,max_level, a90,a90,0.0, 1.0, true,falloff);
    TRACE1 ("Time = %f",timer.value());
  }
#endif


  for (int i=0; i<=tree.max_level(); i++) {
    char filename[40];
    sprintf (filename,"density_3d_1-initial-L%d.png",i);
    timer.clear();  timer.start();
    tree_to_png (&tree,filename,
		 mx,my, i,i, th,ph,ps, scale, false, falloff);
    TRACE1 ("Time = %f",timer.value());
  }

    timer.clear();  timer.start();
  tree_to_png (&tree,"density_xy_1-initial.png",
	       mx,my, min_level,max_level, 0.0,0.0,0.0, 1.0, true,0);
    TRACE1 ("Time = %f",timer.value());
  if (dimension == 3) {
    timer.clear();  timer.start();
    tree_to_png (&tree,"density_yz_1-initial.png",
		 mx,my, min_level,max_level, 0.0,-a90,a90, 1.0, true,0);
    TRACE1 ("Time = %f",timer.value());
    timer.clear();  timer.start();
    tree_to_png (&tree,"density_zx_1-initial.png",
		 mx,my, min_level,max_level, a90,a90,0.0, 1.0, true,0);
    TRACE1 ("Time = %f",timer.value());
  }


  // --------------------------------------------------
  // Balance tree
  // --------------------------------------------------

  timer.clear();
  timer.start();

  tree.balance();
  
  timer.stop();

  printf ("tree balanced time      = %f\n",timer.value());
  printf ("tree balanced num_nodes = %d\n",tree.num_nodes());
  printf ("tree balanced max_level = %d\n",tree.max_level());

  {
    for (int i=0; i<=tree.max_level(); i++) sum_mesh[i]=0;
    ItNode it_node(&tree);
    timer.clear();
    timer.start();
    while (it_node.next_leaf()) {
      ++sum_mesh[it_node.node_trace()->level()];
    }
    timer.stop();
    printf ("traversed %d nodes in %f Mnodes/s\n",tree.num_nodes(),
	    1e-6*tree.num_nodes()/timer.value());
  }
  // --------------------------------------------------
  // Write tree to file
  // --------------------------------------------------

    timer.clear();  timer.start();
  tree_to_png (&tree,"density_3d_2-balanced.png",
	       mx,my,  min_level,max_level, th,ph,ps, scale, false, falloff);
    TRACE1 ("Time = %f",timer.value());

    timer.clear();  timer.start();
  tree_to_png (&tree,"density_xy_2-balanced.png",
	       mx,my, min_level,max_level, 0.0,0.0,0.0, 1.0, true,0);
    TRACE1 ("Time = %f",timer.value());
    timer.clear();  timer.start();
  tree_to_png (&tree,"density_yz_2-balanced.png",
	       mx,my, min_level,max_level, 0.0,-a90,a90, 1.0, true,0);
    TRACE1 ("Time = %f",timer.value());
    timer.clear();  timer.start();
  tree_to_png (&tree,"density_zx_2-balanced.png",
	       mx,my, min_level,max_level, a90,a90,0.0, 1.0, true,0);
    TRACE1 ("Time = %f",timer.value());

  // --------------------------------------------------
  // Patch coalescing
  // --------------------------------------------------

  timer.clear();
  timer.start();


  tree.coalesce();
  timer.stop();

  printf ("tree merged time      = %f\n",timer.value());
  printf ("tree merged num_nodes = %d\n",tree.num_nodes());
  printf ("tree merged max_level = %d\n",tree.max_level());

  // --------------------------------------------------
  // Write tree to file
  // --------------------------------------------------

    timer.clear();  timer.start();
  tree_to_png (&tree,"density_xy_3-coalesced.png",
	       mx,my, min_level,max_level, 0.0,0.0,0.0, 1.0, true,0);
    TRACE1 ("Time = %f",timer.value());

  if (dimension == 3) {
    timer.clear();  timer.start();
    tree_to_png (&tree,"density_yz_3-coalesced.png",
		 mx,my, min_level,max_level, 0.0,-a90,a90, 1.0, true,0);
    TRACE1 ("Time = %f",timer.value());
    timer.clear();  timer.start();
    tree_to_png (&tree,"density_zx_3-coalesced.png",
		 mx,my, min_level,max_level, a90,a90,0.0, 1.0, true,0);
    TRACE1 ("Time = %f",timer.value());

    timer.clear();  timer.start();
    tree_to_png (&tree,"density_3d_3-coalesced.png",
		 mx,my,  min_level,max_level, th,ph,ps, scale, false, falloff);
    TRACE1 ("Time = %f",timer.value());
  }


  // --------------------------------------------------
  // Targeted Refineent
  // --------------------------------------------------

  {
    Timer timer;
    timer.start();
    int d=3;
    int r=2;
    Tree tree (d,r);
    timer.clear();  timer.start();
    bool target = true;
    levels_to_tree (&tree, levels,nx,ny,nz,target);
    TRACE1 ("Time = %f",timer.value());

    // --------------------------------------------------
    // Write tree to file
    // --------------------------------------------------

    timer.clear();  timer.start();
    tree_to_png (&tree,"density_xy_4-targeted.png",
		 mx,my, min_level,max_level, 0.0,0.0,0.0, 1.0, true,0,target);


    if (dimension == 3) {
      timer.clear();  timer.start();
      tree_to_png (&tree,"density_yz_4-targeted.png",
		   mx,my, min_level,max_level, 0.0,-a90,a90, 1.0, true,0);
      TRACE1 ("Time = %f",timer.value());
      timer.clear();  timer.start();
      tree_to_png (&tree,"density_zx_4-targeted.png",
		   mx,my, min_level,max_level, a90,a90,0.0, 1.0, true,0);
      TRACE1 ("Time = %f",timer.value());

      timer.clear();  timer.start();
      tree_to_png (&tree,"density_3d_4-targeted.png",
		   mx,my,  min_level,max_level, th,ph,ps, scale, false, falloff);
      TRACE1 ("Time = %f",timer.value());

    }
    ItNode it_node(&tree);
    int num_nodes = 0;
    int max_level = 0;
    int * level_count = new int [tree.max_level()+1];
    Node * node = 0;
    while ((node = it_node.next_node())) {
      if (node->data() != 0) {
	int level = *((int *)node->data());
	++ num_nodes;
	++ level_count[level];
	if (level > max_level) max_level = level;
      }
    }
    printf ("tree targeted time      = %f\n",timer.value());
    printf ("tree targeted num_nodes = %d\n",num_nodes);
    printf ("tree targeted max_level = %d\n",max_level);
    delete [] level_count;
    // --------------------------------------------------
  }

  unit_finalize();

  delete [] levels;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

