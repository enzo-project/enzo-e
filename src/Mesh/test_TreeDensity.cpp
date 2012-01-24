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
  // Read density_512.h5 file into density[] array
  //--------------------------------------------------

  FileHdf5 file_density ("input","density_512.h5");

  file_density.file_open();

  // H5T_IEEE_F32BE
  int nx,ny,nz;
  scalar_type type = scalar_type_unknown;
  file_density.data_open ("Density",&type,&nx,&ny,&nz);

  float * density = new float [nx*ny*nz];
  file_density.data_read(density);

  unit_assert (false);

  //--------------------------------------------------
  // Refine on density
  //--------------------------------------------------

  int d=3;
  int r=2;
  int min_level = 0;
  int max_level = 10;
  Tree tree (d,r);

  // find the min and max
  float dmin   =1.0e37;
  float dmax = -1.0e37;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + nx*(iy + ny*iz);
	if (density[i] < dmin) dmin = density[i];
	if (density[i] > dmax) dmax = density[i];
      }
    }
  }
  TRACE2 ("min = %f  max = %f",dmin,dmax);

  // create level array from density

  int * levels = new int [nx*ny*nz];

  // linear interpolate log density between minimum level and maximum level

  float lg_dmin = log(dmin);
  float lg_dmax = log(dmax);
  float mult = 1.0*(max_level - min_level) / (lg_dmax - lg_dmin);

  int c=0;
  int imx = -1000, imn=1000;
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + nx*(iy + ny*iz);
	float lg_d = log (density[i]);
	levels[i] = min_level + (lg_d - lg_dmin) * mult;
	if (levels[i] < imn) imn = levels[i];
	if (levels[i] > imx) imx = levels[i];
      }
    }
  }
  TRACE3 ("%d %d %d ",imn,imx,c);

  Timer timer;
  timer.start();
  create_tree_from_levels (&tree, levels,nx,ny,nz);

  TRACE1 ("Tree create time = %f",timer.value());
  TRACE1 ("Tree (unbalanced) nodes = %d",tree.num_nodes());

  //  int mx=2048,my=2048;
  int mx=2048,my=2048;
  double th= 0.3*M_PI; // spin
  double ph= 0.2*M_PI;
  double ps= 0.0*M_PI;

  create_image_from_tree (&tree,"density.png",
			  mx,my, 0,max_level, th,ph,ps, 0.5, true);

  timer.clear();
  timer.start();

  tree.balance();
  
  TRACE1 ("Tree balance time = %f",timer.value());
  TRACE1 ("Tree (balanced) nodes = %d",tree.num_nodes());


  create_image_from_tree (&tree,"density_balanced.png",
			  mx,my,  0,max_level, th,ph,ps, 0.5, true);

  unit_finalize();

  delete [] levels;
  delete [] density;

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

