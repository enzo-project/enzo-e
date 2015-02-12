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

void print_grid(const char * name, double * grid, int nd3[3])
{
  for (int iy=0; iy<nd3[1]; iy++) {
    for (int iz=0; iz<nd3[2]; iz++) {
      for (int ix=0; ix<nd3[0]; ix++) {
	int i = ix + nd3[0]*(iy + nd3[1]*iz);
	if (grid[i] == CLEAR)
	  PARALLEL_PRINTF ("  -  ");
	else
	  PARALLEL_PRINTF ("%s %d %d %d %10.5f ",name,ix,iy,iz,grid[i]);
	PARALLEL_PRINTF("\n");
      }
    }
  }
}

void set_grid(double * grid, int nd3[3], 
	      double c, double x, double y, double z)
{
  for (int iz=0; iz<nd3[2]; iz++) {
    for (int iy=0; iy<nd3[1]; iy++) {
      for (int ix=0; ix<nd3[0]; ix++) {
	int i = ix + nd3[0]*(iy + nd3[1]*iz);
	grid[i] = c + ix*x + iy*y + iz*z;
      }
    }
  }
}

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  //--------------------------------------------------

  const int g=3;
  int nc3[3] = {8,8,8};
  int nf3[3] = {12,12,12};
  int mf3[3] = {12,12,12};

  double * ac = new double [nc3[0]*nc3[1]*nc3[2]];
  double * af = new double [mf3[0]*mf3[1]*mf3[2]];

  int dc = g + nc3[0]*(g + nc3[1]*g);
  int df = g + mf3[0]*(g + mf3[1]*g);

  set_grid  (ac,nc3,1,10,100,1);    // 1 + 1*ix + 2* iy + 3*iz
  print_grid("coarse before interp3d",ac,nc3);
  set_grid  (af,mf3,CLEAR,0,0,0); // CLEAR + 0*ix + 0*iy + 0*iz


  int dim[3]   = {nc3[0],nc3[1],nc3[2]};
  int start[3] = {3,3,3};
  int end[3]   = {2+2*(nc3[0]-2),
		  2+2*(nc3[1]-2),
		  2+2*(nc3[2]-2)};
  int refine[3] = {2,2,2};
  int gdim[3] = {end[0]-2,end[1]-2,end[2]-2};
  int gstart[3] = {1,1,1};
  int wdim[3] = {dim[0]-1,dim[1]-1,dim[2]-1};
  double * w = new double[wdim[0]*wdim[1]*wdim[2]];
  int ierror = 0;

  FORTRAN_NAME(interp3d)(ac,w,
			 &dim[0],&dim[1],&dim[2],
			 &start[0],&start[1],&start[2],
			 &end[0],&end[1],&end[2],
			 &refine[0],&refine[1],&refine[2],
			 af,
			 &gdim[0],&gdim[1],&gdim[2],
			 &gstart[0],&gstart[1],&gstart[2],
			 &wdim[0],&wdim[1],&wdim[2],
			 &ierror);

  print_grid("fine after interp3d",af,mf3);

  //======================================================================
  unit_finalize();
  PARALLEL_EXIT;
  //======================================================================


  set_grid  (af,nf3,CLEAR,0,0,0); // CLEAR + 0*ix + 0*iy + 0*iz
  
  ProlongLinear prolong_linear;
  int im3[3] = {0,0,0};
  prolong_linear.apply(precision_double,
		       af, nf3, im3,nf3,
		       ac, nc3, im3,nc3);

  //  print_grid("fine after ProlongLinear",af,mf3);
		       
		       
  //======================================================================
  unit_finalize();
  PARALLEL_EXIT;
  //======================================================================

  {
    FieldDescr field_descr;

    int i_d = field_descr.insert_permanent("density");
    int i_v = field_descr.insert_permanent("velocity_x");

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

    FieldData field_data_f (&field_descr);
    FieldData field_data_c (&field_descr);

    field_data_f.allocate_permanent(true);

    double * df = (double *) field_data_f.values(i_d);
    float  * vf = (float  *) field_data_f.values(i_v);

    field_data_c.allocate_permanent(true);

    double * dc = (double *) field_data_c.values(i_d);
    float * vc  = (float *) field_data_c.values(i_v);

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
    field_data_f.print("AF-1");
    field_data_c.print("AC");

    char buffer[80];
    for (int icx=0; icx<2; icx++) {
      for (int icy=0; icy<2; icy++) {
	for (int id=0; id<md; id++) df[id]=11111.0;
	for (int iv=0; iv<mv; iv++) vf[iv]=22222.0;
	// prolong_linear->apply
	// 	(&field_data_f,&field_data_c, &field_descr, icx,icy,0);
	sprintf (buffer,"ProlongLinear-%d%d0.out",icx,icy);
	field_data_f.print(buffer);

	for (int id=0; id<md; id++) df[id]=11111.0;
	for (int iv=0; iv<mv; iv++) vf[iv]=22222.0;
	// enzo_prolong->apply
	// 	(&field_data_f,&field_data_c,  icx,icy,0);
	sprintf (buffer,"EnzoProlong-%d%d0.out",icx,icy);
	field_data_f.print(buffer);
      }
    }

    unit_assert (false);

    //--------------------------------------------------

    unit_finalize();

    //    exit_();
  }
}

PARALLEL_MAIN_END


//======================================================================
#include "enzo.def.h"
//======================================================================
