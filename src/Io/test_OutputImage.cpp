// See LICENSE_CELLO file for license and copyright information

/// @file     test_OutputImage.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the OutputImage class

#include "main.hpp" 
#include "test.hpp"

#include "io.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  unit_class("OutputImage");

  OutputImage * output = new OutputImage;

  unit_func ("OutputImaged");

  unit_assert (output != NULL)

    // CODE FORMERLY IN test_Monitor.cpp

  // int n = 128;

  // // Allocate array

  // if (parallel->rank()==0) {
  //   PARALLEL_PRINTF ("pngwriter version = %g\n",pngwriter::version());
  // }

  // float * array = new float [n*n*n];
  // for (int i=0; i<n*n*n; i++) array[i]=0;

  // // Set values to 1 for radius between n/8 and n/4 (note n/2 is boundary)
  // for (int iz=0; iz<n; iz++) {
  //   double z = (iz - 0.5*n) / n;
  //   for (int iy=0; iy<n; iy++) {
  //     double y = (iy - 0.5*n) / n;
  //     for (int ix=0; ix<n; ix++) {
  // 	int i = ix + n*(iy + n*iz);
  // 	double x = (ix - 0.5*n) / n;
  // 	//	array[i] = y;
  // 	double r = sqrt(x*x + y*y + z*z);
  // 	if (0.25 < r && r < 0.5) array[i] = 1;
  //     }
  //   }
  // }


  // unit_func("image");

  // {
  //   output->image("output_image_1.png",
  // 		   n,n,
  // 		   array,
  // 		   n,n,n,
  // 		   n,n,n,
  // 		   0,0,0, axis_x,reduce_sum,0,1);
  //   unit_assert(true);
  // }

  // {
  //   double map_r[] = {0.0, 1.0};
  //   double map_g[] = {0.0, 0.5};
  //   double map_b[] = {0.5, 1.0};
  //   output->image_set_map(2,map_r,map_g,map_b);
  //   output->image("output_image_2.png",
  // 		   n,n,
  // 		   array,
  // 		   n,n,n,
  // 		   n,n,n,
  // 		   0,0,0, axis_y,reduce_sum,0,1);
  //   unit_assert(true);
  // }

  // {
  //   double map_r[] = {0.0, 1.0, 0.0, 0.0};
  //   double map_g[] = {0.0, 0.0, 1.0, 0.0};
  //   double map_b[] = {0.0, 0.0, 0.0, 1.0};
  //   output->image_set_map(4,map_r,map_g,map_b);
  //   output->image("output_image_3.png",
  // 		   n,n,
  // 		   array,
  // 		   n,n,n,
  // 		   n,n,n,
  // 		   0,0,0, axis_z,reduce_sum,0,1);
  //   unit_assert(true);
  // }

  // {
  //   double map_r[] = {0.0, 1.0};
  //   double map_g[] = {0.0, 1.0};
  //   double map_b[] = {0.0, 1.0};
  //   output->image_set_map(2,map_r,map_g,map_b);

  //   output->image_open("output_image_4.png",n,n);

  //   for (int iz=0; iz<2; iz++) {
  //     int iz0 = iz*n/2;
  //     for (int iy=0; iy<2; iy++) {
  // 	int iy0 = iy*n/2;
  // 	for (int ix=0; ix<2; ix++) {
  // 	  int ix0 = ix*n/2;
  // 	  int i = ix0 + n*(iy0 + n*iz0);
  // 	  output->image_reduce(array+i,
  // 			       n,n,n,
  // 			       n/2,n/2,n/2,
  // 			       ix0,iy0,iz0,
  // 			       axis_z,reduce_avg);
  // 	}
  //     }
  //   }

  //   output->image_close(0,1);

  //   unit_assert(true);
  // }


  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END


