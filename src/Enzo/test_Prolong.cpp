// See LICENSE_CELLO file for license and copyright information

/// @file     test_Prolong.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-11
/// @brief    Test program for the Prolong class

#define CHARM_ENZO

#include "test.hpp"

#include "main.hpp"
#include "enzo.hpp"

#include "charm_enzo.hpp"

#define CLEAR -1111.0

void print_grid(double * grid, int nd3[3])
{
  PARALLEL_PRINTF("FINE\n");
  for (int iy=0; iy<nd3[1]; iy++) {
    for (int ix=0; ix<nd3[0]; ix++) {
      int i = ix + nd3[0]*iy;
      if (grid[i] == CLEAR)
	PARALLEL_PRINTF ("  -  ");
      else
	PARALLEL_PRINTF ("%4.1f ",grid[i]);
    }
    PARALLEL_PRINTF("\n");
  }
}

void set_grid(double * grid, int nd3[3], 
	      double c, double x, double y)
{
  for (int iy=0; iy<nd3[1]; iy++) {
    for (int ix=0; ix<nd3[0]; ix++) {
      int i = ix + nd3[0]*iy;
      grid[i] = c + ix*x + iy*y;
    }
  }
}

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  //--------------------------------------------------

  int n3[3] = {4,4,4};
  int g3[3] = {3,3,3};
  int nd3[3] = {n3[0]+2*g3[0],
		n3[1]+2*g3[1],
		n3[2]+2*g3[2]};

  double * grid_fine = new double [nd3[0]*nd3[1]];
  double * grid_coarse = new double [nd3[0]*nd3[1]];

  set_grid  (grid_coarse,nd3,1,1,2);     // 1 + 1*ix + 2* iy
  print_grid(grid_coarse,nd3);

  set_grid  (grid_fine,  nd3,CLEAR,0,0); // = CLEAR
  print_grid(grid_fine,nd3);

  // FORTRAN_NAME(interpolate)(&ndim, grid_coarse, nd3, pstart, pend,
  // 			    refine,
  // 			    grid_fine, nd3, gstart, work, &imethod,
  // 			    &iposflag, &error);

  ProlongLinear prolong_linear;
  prolong_linear.apply(precision_double,
		       grid_fine,
		       nd3[0],nd3[1],nd3[2], n3[0],n3[1],n3[2],
		       grid_coarse,
		       nd3[0],nd3[1],nd3[2], n3[0]/2,n3[1]/2,n3[2]/2);

  print_grid(grid_fine,nd3);
		       
		       
  //======================================================================
  PARALLEL_EXIT;
  //======================================================================

  {
    FieldDescr field_descr;

    int i_d = field_descr.insert_field("density");
    int i_v = field_descr.insert_field("velocity_x");

    field_descr.set_precision(i_d,precision_double);
    field_descr.set_precision(i_v,precision_single);

    // ghost zone depths
    const int gd = 3;
    const int gv = 2;
    field_descr.set_ghosts(i_d,gd,gd);
    field_descr.set_ghosts(i_v,gv,gv);

    // block size
    const int nx = 4;
    const int ny = 4;

    // array dimensions
    const int mxd = nx+2*gd;
    const int myd = ny+2*gd;
    const int mxv = nx+2*gv;
    const int myv = ny+2*gv;
    const int md = mxd*myd;
    const int mv = mxv*myv;

    FieldBlock field_block_f (nx,ny);
    FieldBlock field_block_c (nx,ny);

    field_block_f.allocate_array(&field_descr, true);

    double * df = (double *) field_block_f.field_values(i_d);
    float  * vf = (float  *) field_block_f.field_values(i_v);

    field_block_c.allocate_array(&field_descr, true);

    double * dc = (double *) field_block_c.field_values(i_d);
    float * vc  = (float *) field_block_c.field_values(i_v);

    //--------------------------------------------------

    // initialize coarse values

    for (int ix=0; ix<mxd; ix++) {
      for (int iy=0; iy<myd; iy++) {
	int id = ix + mxd * iy;
	dc[id] = -1.0*(ix + 2*iy);
      }
    }
    for (int ix=0; ix<mxv; ix++) {
      for (int iy=0; iy<myv; iy++) {
	int iv = ix + mxv * iy;
	vc[iv] = -1.0*(ix + 2*iy);
      }
    }


    //--------------------------------------------------

    unit_class("ProlongLinear");

    Prolong * prolong_linear = new ProlongLinear;
    unit_assert (prolong_linear != NULL);

    unit_class("EnzoProlong");

    Prolong * enzo_prolong = new EnzoProlong("SecondOrderA");
    unit_assert (enzo_prolong != NULL);

    unit_func ("apply()");

    double lower[3] = { 0.0, 0.0, 0.0 };
    double upper[3] = { 1.0, 1.0, 1.0 };
    field_block_f.print(&field_descr,"AF-1", lower,upper);
    field_block_c.print(&field_descr,"AC", lower,upper);

    char buffer[80];
    for (int icx=0; icx<2; icx++) {
      for (int icy=0; icy<2; icy++) {
	for (int id=0; id<md; id++) df[id]=11111.0;
	for (int iv=0; iv<mv; iv++) vf[iv]=22222.0;
	// prolong_linear->apply
	// 	(&field_block_f,&field_block_c, &field_descr, icx,icy,0);
	sprintf (buffer,"ProlongLinear-%d%d0.out",icx,icy);
	field_block_f.print(&field_descr,buffer, lower,upper);

	for (int id=0; id<md; id++) df[id]=11111.0;
	for (int iv=0; iv<mv; iv++) vf[iv]=22222.0;
	// enzo_prolong->apply
	// 	(&field_block_f,&field_block_c, &field_descr, icx,icy,0);
	sprintf (buffer,"EnzoProlong-%d%d0.out",icx,icy);
	field_block_f.print(&field_descr,buffer, lower,upper);
      }
    }

    unit_assert (false);

    //--------------------------------------------------

    unit_finalize();

    PARALLEL_EXIT;
  }
}

PARALLEL_MAIN_END


//======================================================================
#ifdef CONFIG_USE_CHARM
#  include "enzo.def.h"
#endif
//======================================================================
