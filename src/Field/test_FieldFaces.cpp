// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_FieldFaces.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the FieldFaces class

#include "test.hpp"

#include "field.hpp"

#include PARALLEL_CHARM_INCLUDE(test_FieldFaces.decl.h)


//----------------------------------------------------------------------

bool is_ghost(int ix,int iy,int iz,
	      int gx,int gy,int gz,
	      int mdx,int mdy,int mdz)
{
  return (gx <= ix && ix < mdx - gx && 
	  gy <= iy && iy < mdy - gy && 
	  gz <= iz && iz < mdz - gz);
}

//----------------------------------------------------------------------

double value (int ix, int iy, int iz,
	      int gx,int gy,int gz,
	      int mdx, int mdy, int mdz,
	      int k, int field, int M3, int N3)
{
  int i = ix + mdx * (iy + mdy * iz);
  double v = i + M3 * (k + N3*field);
  return is_ghost(ix,iy,iz,gx,gy,gz,mdx,mdy,mdz) ? v : -v;
}
//----------------------------------------------------------------------

void init_field_1(FieldBlock ** field_block,
		  int k, int i, 
		  int mx, int my, int mz,
		  int M3, int N3)
{

  const FieldDescr * field_descr = field_block[k]->field_descr();

  int gx, gy, gz;

  field_descr->ghosts(i, &gx, &gy, &gz);

  int mdx, mdy, mdz;

  mdx = mx + 2*gx;
  mdy = my + 2*gy;
  mdz = mz + 2*gz;

  float * v = (float *) field_block[k]->field_values(i);

  for (int iz = 0; iz < mdz; iz++) {
    for (int iy = 0; iy < mdy; iy++) {
      for (int ix = 0; ix < mdx; ix++) {
	int i = ix + mdx * (iy + mdy * iz);
	v[i] = value(ix,iy,iz,gx,gy,gz,mdx,mdy,mdz,k,0,M3,N3);
      }
    }
  }
}

//----------------------------------------------------------------------

void init_field_2(FieldBlock ** field_block,
		  int k, int i, 
		  int mx, int my, int mz,
		  int M3, int N3)
{

  const FieldDescr * field_descr = field_block[k]->field_descr();

  int gx, gy, gz;

  field_descr->ghosts(i, &gx, &gy, &gz);

  int mdx, mdy, mdz;

  mdx = mx + 2*gx;
  mdy = my + 2*gy;
  mdz = mz + 2*gz;

  double * v = (double *) field_block[k]->field_values(i);

  for (int iz = 0; iz < mdz; iz++) {
    for (int iy = 0; iy < mdy; iy++) {
      for (int ix = 0; ix < mdx; ix++) {
	int i = ix + mdx * (iy + mdy * iz);
	v[i] = value(ix,iy,iz,gx,gy,gz,mdx,mdy,mdz,k,1,M3,N3);
      }
    }
  }

}

//----------------------------------------------------------------------

void init_field_3(FieldBlock ** field_block,
		  int k, int i, 
		  int mx, int my, int mz,
		  int M3, int N3)
{

  const FieldDescr * field_descr = field_block[k]->field_descr();

  int gx, gy, gz;

  field_descr->ghosts(i, &gx, &gy, &gz);

  int mdx, mdy, mdz;

  mdx = mx + 2*gx;
  mdy = my + 2*gy;
  mdz = mz + 2*gz;

  long double * v = (long double *) field_block[k]->field_values(i);

  for (int iz = 0; iz < mdz; iz++) {
    for (int iy = 0; iy < mdy; iy++) {
      for (int ix = 0; ix < mdx; ix++) {
	int i = ix + mdx * (iy + mdy * iz);
	v[i] = value(ix,iy,iz,gx,gy,gz,mdx,mdy,mdz,k,2,M3,N3);
      }
    }
  }

}

//======================================================================

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  //--------------------------------------------------
  // Initialize the global field descriptor object field_descr
  //--------------------------------------------------

  FieldDescr * field_descr = new FieldDescr;

  // insert fields

  int i1 = field_descr->insert_field("field_1");
  int i2 = field_descr->insert_field("field_2");
  int i3 = field_descr->insert_field("field_3");

  // initialize field precisions

  field_descr->set_precision(i1, precision_single);
  field_descr->set_precision(i2, precision_double);
  field_descr->set_precision(i3, precision_quadruple);

  // initialize field ghost zone depths

  field_descr->set_ghosts(i1, 1,1,1);
  field_descr->set_ghosts(i2, 0,1,2);
  field_descr->set_ghosts(i3, 1,0,3);

  const int mx = 4, my = 5, mz = 6;

  const int M3 = (mx+2*1)*(my+2*1)*(mz+2*3);  // maximum size with ghosts

  //--------------------------------------------------
  // Initialize the field blocks in 4 x 4 x 4 array
  //--------------------------------------------------

  const int N  = 4;
  const int N3 = N*N*N;

  FieldBlock * field_block[N3];

  for (int kz = 0; kz < N; kz++) {
    for (int ky = 0; ky < N; ky++) {
      for (int kx = 0; kx < N; kx++) {

	int k = kx + N * (ky + N * kz);

	// Create the FaceBlock object

	field_block[k] = new FieldBlock (field_descr, mx, my, mz);

	// Allocate field blocks including ghosts

	field_block[k]->allocate_array();
	field_block[k]->allocate_ghosts();

	// Initialize fields
	// (each field initialized separately since each has different
	// ghost zone depth)

	init_field_1(field_block,k,i1,mx,my,mz,M3,N3);
	init_field_2(field_block,k,i2,mx,my,mz,M3,N3);
	init_field_3(field_block,k,i3,mx,my,mz,M3,N3);

      }
    }
  }

  // Send/recv faces in xp direction

  INCOMPLETE("test_FieldFaces");
  // send xp faces

  for (int kz = 0; kz < N; kz++) {
    for (int ky = 0; ky < N; ky++) {
      for (int kx = 0; kx < N-1; kx++) {
	int k = kx + N * (ky + N * kz);
	if (field_block[k] != NULL) {
	  FieldFaces * ff = field_block[k]->field_faces();
	}

      }
    }
  }

  // receive xm ghosts
	
  for (int kz = 0; kz < N; kz++) {
    for (int ky = 0; ky < N; ky++) {
      for (int kx = 1; kx < N; kx++) {
      }
    }
  }

  // Test values

  for (int kz = 0; kz < N; kz++) {
    for (int ky = 0; ky < N; ky++) {
      for (int kx = 1; kx < N; kx++) {
      }
    }
  }

  // Send all faces

  //   send_init();
  //   send_begin();
  //   send_end();
  //   send_final();
	
  //   recv_init();
  //   recv_begin();
  //   recv_end();
  //   recv_final();
	
  //   sendrecv_init();
  //   sendrecv_begin();
  //   sendrecv_end();
  //   sendrecv_final();

  unit_finalize();
  //----------------------------------------------------------------------

  PARALLEL_EXIT;
}
PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_FieldFaces.def.h)
