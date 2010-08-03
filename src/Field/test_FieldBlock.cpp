// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_FieldBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the FieldBlock class

#include <math.h>

#include "cello.hpp"

#include "error.hpp"
#include "test.hpp"
#include "field.hpp"

int main()
{

  //----------------------------------------------------------------------
  unit_init();
  //----------------------------------------------------------------------

  FieldDescr field_descr;

  int index_f1 = field_descr.insert_field("f1");
  int index_f2 = field_descr.insert_field("f2");
  int index_f3 = field_descr.insert_field("f3");
  int index_f4 = field_descr.insert_field("f4");
  int index_f5 = field_descr.insert_field("f5");

  field_descr.set_precision(index_f1, precision_single);
  field_descr.set_precision(index_f2, precision_double);
  field_descr.set_precision(index_f3, precision_double);
  field_descr.set_precision(index_f4, precision_double);
  field_descr.set_precision(index_f5, precision_quadruple);

  int g1[3] = {1,1,1};
  int g2[3] = {2,2,2};
  int g3[3] = {1,0,2};
  int g4[3] = {0,0,0};
  int g5[3] = {3,3,3};
  field_descr.set_ghosts(index_f1, g1[0],g1[1],g1[2]);
  field_descr.set_ghosts(index_f2, g2[0],g2[1],g2[2]);
  field_descr.set_ghosts(index_f3, g3[0],g3[1],g3[2]);
  field_descr.set_ghosts(index_f4, g4[0],g4[1],g4[2]);
  field_descr.set_ghosts(index_f5, g5[0],g5[1],g5[2]);

  field_descr.set_centering(index_f2, false, true,  true);
  field_descr.set_centering(index_f3, true,  false, true);
  field_descr.set_centering(index_f4, true,  true,  false);

  unit_class ("[precision]");
  unit_func ("sizeof(float)");
  unit_assert (sizeof(float)       == 4);
  unit_func ("sizeof(double)");
  unit_assert (sizeof(double)      == 8);
  unit_func ("sizeof(long-double)");
  unit_assert (sizeof(long double) == 16);

  //----------------------------------------------------------------------
  unit_class ("FieldBlock");
  //----------------------------------------------------------------------

  FieldBlock field_block;

  //----------------------------------------------------------------------

  unit_func("field_descr");

  field_block.set_field_descr(&field_descr);
  unit_assert (field_block.field_descr() == & field_descr);

  //----------------------------------------------------------------------

  unit_func("dimensions");

  int nx,ny,nz;
  nx=4; ny=5; nz=6;
  field_block.set_dimensions(nx,ny,nz);
  int dimensions[3];
  field_block.dimensions(&dimensions[0],&dimensions[1],&dimensions[2]);
  unit_assert(dimensions[0]==nx && dimensions[1]==ny && dimensions[2]==nz);

  nx=5; ny=3; nz=4;
  field_block.set_dimensions(nx,ny,nz);
  field_block.dimensions(&dimensions[0],&dimensions[1],&dimensions[2]);
  unit_assert(dimensions[0]==nx && dimensions[1]==ny && dimensions[2]==nz);

  //----------------------------------------------------------------------
  // allocate / deallocate
  //----------------------------------------------------------------------

  unit_func("array_allocated");
  unit_assert( ! field_block.array_allocated());

  // Allocate

  unit_func("allocate_array");
  field_block.allocate_array();
  unit_assert(field_block.array() != 0);

  unit_func("array_allocated");
  unit_assert( field_block.array_allocated());

  // Deallocate

  unit_func("deallocate_array");
  field_block.deallocate_array();
  unit_assert(field_block.array() == 0);

  unit_func("array_allocate");
  unit_assert( ! field_block.array_allocated());

  // Allocate

  unit_func("allocate_array");
  field_block.allocate_array();
  unit_assert(field_block.array() != 0);
  
  unit_func("array_allocated");
  unit_assert( field_block.array_allocated());

  
  //----------------------------------------------------------------------

  float  *      values_f1;
  double *      values_f2;
  double *      values_f3;
  double *      values_f4;
  long double * values_f5;

  float *      unknowns_f1;
  double *     unknowns_f2;
  double *     unknowns_f3;
  double *     unknowns_f4;
  long double * unknowns_f5;

  int bytes_f1;
  int bytes_f2;
  int bytes_f3;
  int bytes_f4;

  unit_func("field_values");  // without ghosts
  
  unit_assert((values_f1 = (float *)      field_block.field_values(index_f1)));
  unit_assert((values_f2 = (double *)     field_block.field_values(index_f2)));
  unit_assert((values_f3 = (double *)     field_block.field_values(index_f3)));
  unit_assert((values_f4 = (double *)     field_block.field_values(index_f4)));
  unit_assert((values_f5 = (long double *) field_block.field_values(index_f5)));
  
  bytes_f1 = (char *)values_f2 - (char *)values_f1;
  bytes_f2 = (char *)values_f3 - (char *)values_f2;
  bytes_f3 = (char *)values_f4 - (char *)values_f3;
  bytes_f4 = (char *)values_f5 - (char *)values_f4;

  // field sizes without ghosts
  int num_unknowns_f1  = (nx)*(ny)*(nz);
  int num_unknowns_f2  = (nx+1)*(ny)*(nz);
  int num_unknowns_f3  = (nx)*(ny+1)*(nz);
  int num_unknowns_f4  = (nx)*(ny)*(nz+1);

  // field sizes with ghosts

  int num_values_f1  = (nx+2*g1[0])  *(ny+2*g1[1])  *(nz+2*g1[2]);
  int num_values_f2  = (nx+2*g2[0]+1)*(ny+2*g2[1])  *(nz+2*g2[2]);
  int num_values_f3   = (nx+2*g3[0])  *(ny+2*g3[1]+1)*(nz+2*g3[2]);
  int num_values_f4   = (nx+2*g4[0])  *(ny+2*g4[1])  *(nz+2*g4[2]+1);
  
  unit_assert (bytes_f1 == (int)sizeof (float) * num_unknowns_f1);
  unit_assert (bytes_f2 == (int)sizeof (double)* num_unknowns_f2);
  unit_assert (bytes_f3 == (int)sizeof (double)* num_unknowns_f3);
  unit_assert (bytes_f4 == (int)sizeof (double)* num_unknowns_f4);

  //----------------------------------------------------------------------

  unit_func("field_unknowns");  // without ghosts

  typedef float type_f1;
  unknowns_f1 = (type_f1 *)      field_block.field_unknowns(index_f1);
  unknowns_f2 = (double *)     field_block.field_unknowns(index_f2);
  unknowns_f3 = (double *)     field_block.field_unknowns(index_f3);
  unknowns_f4 = (double *)     field_block.field_unknowns(index_f4);
  unknowns_f5 = (long double *) field_block.field_unknowns(index_f5);

  unit_assert(unknowns_f1 != 0);
  unit_assert(unknowns_f2 != 0);
  unit_assert(unknowns_f3 != 0);
  unit_assert(unknowns_f4 != 0);
  unit_assert(unknowns_f5 != 0);

  bytes_f1 = (char *)unknowns_f2 - (char *)unknowns_f1;
  bytes_f2 = (char *)unknowns_f3 - (char *)unknowns_f2;
  bytes_f3 = (char *)unknowns_f4 - (char *)unknowns_f3;
  bytes_f4 = (char *)unknowns_f5-(char *)unknowns_f4;

  unit_assert (bytes_f1    == 
	       (int)sizeof (float) * num_unknowns_f1);
  unit_assert (bytes_f2 == 
	       (int)sizeof (double)* num_unknowns_f2);
  unit_assert (bytes_f3 == 
	       (int)sizeof (double)* num_unknowns_f3);
  unit_assert (bytes_f4 == 
	       (int)sizeof (double)* num_unknowns_f4);


  //----------------------------------------------------------------------

  unit_func("allocate_ghosts");  // with ghosts

  unit_assert ( ! field_block.ghosts_allocated());
  field_block.allocate_ghosts();
  unit_assert ( field_block.ghosts_allocated());
  
  values_f1    = 
    (float *) field_block.field_values(index_f1);
  values_f2 = 
    (double *) field_block.field_values(index_f2);
  values_f3 = 
    (double *) field_block.field_values(index_f3);
  values_f4 = 
    (double *) field_block.field_values(index_f4);
  values_f5 =
    (long double *) field_block.field_values(index_f5);
  
  unit_assert(values_f1 != 0);
  unit_assert(values_f2 != 0);
  unit_assert(values_f3 != 0);
  unit_assert(values_f4 != 0);
  unit_assert(values_f5 != 0);

  bytes_f1 = (char *)values_f2 - (char *)values_f1;
  bytes_f2 = (char *)values_f3 - (char *)values_f2;
  bytes_f3 = (char *)values_f4 - (char *)values_f3;
  bytes_f4 = (char *)values_f5-(char *)values_f4;

  unit_assert (bytes_f1    == 
	       (int)sizeof (float) * num_values_f1);
  unit_assert (bytes_f2 == 
	       (int)sizeof (double)* num_values_f2);
  unit_assert (bytes_f3 == 
	       (int)sizeof (double)* num_values_f3);
  unit_assert (bytes_f4 == 
	       (int)sizeof (double)* num_values_f4);

  unit_func("field_unknowns");  // with ghosts

  unknowns_f1    = 
    (float *) field_block.field_unknowns(index_f1);

  unknowns_f2 = 
    (double *) field_block.field_unknowns(index_f2);
  unknowns_f3 = 
    (double *) field_block.field_unknowns(index_f3);
  unknowns_f4 = 
    (double *) field_block.field_unknowns(index_f4);
  unknowns_f5 =
    (long double *) field_block.field_unknowns(index_f5);

  unit_assert(unknowns_f1 != 0);
  unit_assert(unknowns_f2 != 0);
  unit_assert(unknowns_f3 != 0);
  unit_assert(unknowns_f4 != 0);
  unit_assert(unknowns_f5 != 0);

  bytes_f1 =    (char *)unknowns_f2 - (char *)unknowns_f1;
  bytes_f2 = (char *)unknowns_f3 - (char *)unknowns_f2;
  bytes_f3 = (char *)unknowns_f4 - (char *)unknowns_f3;
  bytes_f4 = (char *)unknowns_f5-(char *)unknowns_f4;

  // a,b fields  u unknowns  v values  g ghosts
  // bu - au = (bv + bg) - (av + ag) 
  //         = (bv + (bu-bv)) - (av + (au-av))
  //         = (bv - av) + (bu-bv) - (au-av)


  int n1[3] = { nx+2*g1[0],   ny+2*g1[1],   nz+2*g1[2]};
  int n2[3] = { nx+2*g2[0]+1, ny+2*g2[1],   nz+2*g2[2]};
  int n3[3] = { nx+2*g3[0],   ny+2*g3[1]+1, nz+2*g3[2]};
  int n4[3] = { nx+2*g4[0],   ny+2*g4[1],   nz+2*g4[2]+1};
  int n5[3] = { nx+2*g5[0],   ny+2*g5[1],   nz+2*g5[2]};

  int go_f1 = (int)sizeof (float)      * (g1[0] + n1[0] *(g1[1] + n1[1]*g1[2]));
  int go_f2 = (int)sizeof (double)     * (g2[0] + n2[0] *(g2[1] + n2[1]*g2[2]));
  int go_f3 = (int)sizeof (double)     * (g3[0] + n3[0] *(g3[1] + n3[1]*g3[2]));
  int go_f4 = (int)sizeof (double)     * (g4[0] + n4[0] *(g4[1] + n4[1]*g4[2]));
  int go_f5 = (int)sizeof (long double)* (g5[0] + n5[0] *(g5[1] + n5[1]*g5[2]));

  unit_assert (bytes_f1 == (int)sizeof (float) * num_values_f1 + go_f2 - go_f1);
  unit_assert (bytes_f2 == (int)sizeof (double)* num_values_f2 + go_f3 - go_f2);
  unit_assert (bytes_f3 == (int)sizeof (double)* num_values_f3 + go_f4 - go_f3);
  unit_assert (bytes_f4 == (int)sizeof (double)* num_values_f4 + go_f5 - go_f4);

  int passed;

  //--------------------------------------------------
  // Fill values interior then test against unknowns

  for (int iz=0; iz<n1[2]; iz++) {
    for (int iy=0; iy<n1[1]; iy++) {
      for (int ix=0; ix<n1[0]; ix++) {
	int i = ix + n1[0]*(iy + n1[1]*iz);
	values_f1[i] = 
	  (g1[0] <= ix && ix < n1[0]-g1[0] &&
	   g1[1] <= iy && iy < n1[1]-g1[1] &&
	   g1[2] <= iz && iz < n1[2]-g1[2]) ? 10 : -10;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n1[0]*(iy + n1[1]*iz);
	if (unknowns_f1[i] != 10) passed = false;
      }
    }
  }

  unit_assert(passed);

  //--------------------------------------------------

  for (int iz=0; iz<n2[2]; iz++) {
    for (int iy=0; iy<n2[1]; iy++) {
      for (int ix=0; ix<n2[0]; ix++) {
	int i = ix + n2[0]*(iy + n2[1]*iz);
	values_f2[i] = 
	  (g2[0] <= ix && ix < n2[0]-g2[0] &&
	   g2[1] <= iy && iy < n2[1]-g2[1] &&
	   g2[2] <= iz && iz < n2[2]-g2[2]) ? 20 : -20;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n2[0]*(iy + n2[1]*iz);
	if (unknowns_f2[i] != 20) passed = false;
      }
    }
  }

  unit_assert(passed);

  //--------------------------------------------------

  for (int iz=0; iz<n3[2]; iz++) {
    for (int iy=0; iy<n3[1]; iy++) {
      for (int ix=0; ix<n3[0]; ix++) {
	int i = ix + n3[0]*(iy + n3[1]*iz);
	values_f3[i] = 
	  (g3[0] <= ix && ix < n3[0]-g3[0] &&
	   g3[1] <= iy && iy < n3[1]-g3[1] &&
	   g3[2] <= iz && iz < n3[2]-g3[2]) ? 30 : -30;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n3[0]*(iy + n3[1]*iz);
	if (unknowns_f3[i] != 30) passed = false;
      }
    }
  }

  unit_assert(passed);

  //--------------------------------------------------

  for (int iz=0; iz<n4[2]; iz++) {
    for (int iy=0; iy<n4[1]; iy++) {
      for (int ix=0; ix<n4[0]; ix++) {
	int i = ix + n4[0]*(iy + n4[1]*iz);
	values_f4[i] = 
	  (g4[0] <= ix && ix < n4[0]-g4[0] &&
	   g4[1] <= iy && iy < n4[1]-g4[1] &&
	   g4[2] <= iz && iz < n4[2]-g4[2]) ? 40 : -40;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n4[0]*(iy + n4[1]*iz);
	if (unknowns_f4[i] != 40) passed = false;
      }
    }
  }

  unit_assert(passed);

  //--------------------------------------------------


  for (int iz=0; iz<n5[2]; iz++) {
    for (int iy=0; iy<n5[1]; iy++) {
      for (int ix=0; ix<n5[0]; ix++) {
	int i = ix + n5[0]*(iy + n5[1]*iz);
	values_f5[i] = 
	  (g5[0] <= ix && ix < n5[0]-g5[0] &&
	   g5[1] <= iy && iy < n5[1]-g5[1] &&
	   g5[2] <= iz && iz < n5[2]-g5[2]) ? 50 : -50;

      }
    }
  }

  passed = true;

  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = ix + n5[0]*(iy + n5[1]*iz);
	if (unknowns_f5[i] != 50) passed = false;
      }
    }
  }

  unit_assert(passed);

  //--------------------------------------------------

  // test that first value after f1 is f2 ghost, etc.

  unit_assert(*((double *)(values_f1+n1[0]*n1[1]*n1[2]))     == -20);
  unit_assert(*((double *)(values_f2+n2[0]*n2[1]*n2[2]))     == -30);
  unit_assert(*((double *)(values_f3+n3[0]*n3[1]*n3[2]))     ==  40);  // g4 has no ghosts
  unit_assert(*((long double *)(values_f4+n4[0]*n4[1]*n4[2])) == -50);

  unit_assert(*((float *) (values_f2-1)) == -10);
  unit_assert(*((double *) (values_f3-1)) == -20);
  unit_assert(*((double *) (values_f4-1)) == -30);
  unit_assert(*((double *) (values_f5-1)) ==  40); // g4 has no ghosts

  //----------------------------------------------------------------------
  unit_func("box_extent");

  field_block.set_box_extent(-1, 1, -2, 2, -3, 3);
  double box_lower[3];
  double box_upper[3];
  field_block.box_extent(&box_lower[0],&box_upper[0],
			 &box_lower[1],&box_upper[1],
			 &box_lower[2],&box_upper[2]);

  unit_assert(box_lower[0] == -1.0);
  unit_assert(box_upper[0] == 1.0);
  unit_assert(box_lower[1] == -2.0);
  unit_assert(box_upper[1] == 2.0);
  unit_assert(box_lower[2] == -3.0);
  unit_assert(box_upper[2] == 3.0);

  //----------------------------------------------------------------------
  unit_func("cell_width");

  double hx=0,hy=0,hz=0;

  field_block.cell_width(&hx,&hy,&hz);

  unit_assert(fabs(hx-2.0/nx) < 1e-6);
  unit_assert(fabs(hy-4.0/ny) < 1e-6);
  unit_assert(fabs(hz-6.0/nz) < 1e-6);

	
  //----------------------------------------------------------------------
  unit_func("clear");

  field_block.clear(2.0);

  unit_assert(2.0 == values_f1[0]);
  unit_assert(2.0 == values_f5[n5[0]*n5[1]*n5[2]-1]);

  field_block.clear(3.0, 1 );
  
  unit_assert(2.0 == values_f1[n1[0]*n1[1]*n1[2]-1]);
  unit_assert(3.0 == values_f2[0] );
  unit_assert(3.0 == values_f2[n2[0]*n2[1]*n2[2]-1]);
  unit_assert(2.0 == values_f3[0] );
	
  field_block.clear(4.0, 2, 3 );

  unit_assert(3.0 == values_f2[n2[0]*n2[1]*n2[2]-1]);
  unit_assert(4.0 == values_f3[0] );
  unit_assert(4.0 == values_f3[n3[0]*n3[1]*n3[2]-1]);
  unit_assert(4.0 == values_f4[0] );
  unit_assert(4.0 == values_f4[n4[0]*n4[1]*n4[2]-1]);
  unit_assert(2.0 == values_f5[0] );

  //----------------------------------------------------------------------
  unit_func("deallocate_ghosts");
  unit_assert( field_block.ghosts_allocated());
  field_block.deallocate_ghosts();

  values_f1    = 
    (float *) field_block.field_values(index_f1);
  values_f2 = 
    (double *) field_block.field_values(index_f2);
  values_f3 = 
    (double *) field_block.field_values(index_f3);
  values_f4 = 
    (double *) field_block.field_values(index_f4);
  values_f5 =
    (long double *) field_block.field_values(index_f5);

  unit_assert( ! field_block.ghosts_allocated());

  unit_assert(3.0 == values_f2[(nx+1)*ny*nz-1]);
  unit_assert(4.0 == values_f3[0] );
  unit_assert(4.0 == values_f3[nx*(ny+1)*nz-1]);
  unit_assert(4.0 == values_f4[0] );
  unit_assert(4.0 == values_f4[nx*ny*(nz+1)-1]);
  unit_assert(2.0 == values_f5[0] );
  
  //----------------------------------------------------------------------
  unit_func("split");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("merge");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_func("read");
  unit_assert(false);
  //----------------------------------------------------------------------
  unit_func("write");
  unit_assert(false);
	
  //----------------------------------------------------------------------
  unit_finalize();
  //----------------------------------------------------------------------
}
