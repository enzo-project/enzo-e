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

double init_value (int ix, int iy, int iz,
	      int gx,int gy,int gz,
	      int mdx, int mdy, int mdz,
	      int index_block, int index_field, 
	      int MD3, int ND3)
{
  int i = ix + mdx * (iy + mdy * iz);
  double value = i + MD3 * (index_block + ND3*index_field);
  return is_ghost(ix,iy,iz,gx,gy,gz,mdx,mdy,mdz) ? -value : value;
}

//----------------------------------------------------------------------

double test_value (int ix, int iy, int iz,
	      int gx,int gy,int gz,
	      int mdx, int mdy, int mdz,
	      int index_block, int index_field, 
	      int MD3, int ND3)
{
  int i = ix + mdx * (iy + mdy * iz);
  double value = i + MD3 * (index_block + ND3*index_field);
  return value;
}

//----------------------------------------------------------------------

template<class T>
void init_field(T * field_values,
		int index_block, int index_field, 
		int mx, int my, int mz,
		int gx, int gy, int gz,
		int ND3)
{

  const int mdx = mx + 2*gx;
  const int mdy = my + 2*gy;
  const int mdz = mz + 2*gz;
  const int MD3 = mdx*mdy*mdz;

  for (int iz = 0; iz < mdz; iz++) {
    for (int iy = 0; iy < mdy; iy++) {
      for (int ix = 0; ix < mdx; ix++) {
	int i = ix + mdx * (iy + mdy * iz);
	field_values[i] = init_value(ix,iy,iz,
				     gx,gy,gz,
				     mdx,mdy,mdz,
				     index_block,index_field,
				     MD3,ND3);
      }
    }
  }
}

//----------------------------------------------------------------------

template<class T>
bool test_field(T * field_values,
		int index_block, int index_field, 
		int mx, int my, int mz,
		int gx, int gy, int gz,
		int ND3)
{

  const int mdx = mx + 2*gx;
  const int mdy = my + 2*gy;
  const int mdz = mz + 2*gz;
  const int MD3 = mdx*mdy*mdz;

  bool result = true;
  for (int iz = 0; iz < mdz; iz++) {
    for (int iy = 0; iy < mdy; iy++) {
      for (int ix = 0; ix < mdx; ix++) {
	int i = ix + mdx * (iy + mdy * iz);
	if (field_values[i] != test_value(ix,iy,iz,
					  gx,gy,gz,
					  mdx,mdy,mdz,
					  index_block,index_field,
					  MD3,ND3))
	  result = false;
      }
    }
  }
  return result;
}

//----------------------------------------------------------------------

void init_fields
(
 FieldDescr * field_descr,
 FieldBlock * field_block[],
 int nx,int ny, int nz,
 int mx, int my, int mz)
{
  //--------------------------------------------------
  // Initialize the field blocks in 4 x 4 x 4 array
  //--------------------------------------------------

  const int ND3 = nx*ny*nz;

  for (int kz = 0; kz < nz; kz++) {
    for (int ky = 0; ky < ny; ky++) {
      for (int kx = 0; kx < nx; kx++) {

	int index_block = kx + nx * (ky + ny * kz);

	// Create the FaceBlock object

	field_block[index_block] = new FieldBlock (field_descr, mx, my, mz);

	// Allocate field blocks including ghosts

	field_block[index_block]->allocate_array();
	field_block[index_block]->allocate_ghosts();

	// Initialize fields

	int gx, gy, gz;

	// field 0
	field_descr->ghosts(0, &gx, &gy, &gz);
	float * v1 = (float *)
	  (field_block[index_block]->field_values(0));
	init_field(v1,index_block,0,mx,my,mz,gx,gy,gz,ND3);

	// field 1
	field_descr->ghosts(1, &gx, &gy, &gz);
	double * v2 = (double *)
	  (field_block[index_block]->field_values(1));
	init_field(v2,index_block,1,mx,my,mz,gx,gy,gz,ND3);

	// field 2
	field_descr->ghosts(2, &gx, &gy, &gz);
	long double * v3 = (long double *)
	  (field_block[index_block]->field_values(2));
	init_field(v3,index_block,2,mx,my,mz,gx,gy,gz,ND3);
 
      }
    }
  }
}

bool test_fields
(
 FieldDescr * field_descr,
 FieldBlock * field_block[],
 int nx,int ny, int nz,
 int mx, int my, int mz)
{
  const int ND3 = nx*ny*nz;

  bool result = true;

  for (int kz = 0; kz < nz; kz++) {
    for (int ky = 0; ky < ny; ky++) {
      for (int kx = 0; kx < nx; kx++) {

	int index_block = kx + nx * (ky + ny * kz);

	int gx, gy, gz;

	// field 0
	field_descr->ghosts(0, &gx, &gy, &gz);
	float * v1 = (float *)
	  (field_block[index_block]->field_values(0));
	result = result && 
	  unit_assert(test_field(v1,index_block,0,mx,my,mz,gx,gy,gz,ND3)); 

	// field 1
	field_descr->ghosts(1, &gx, &gy, &gz);
	double * v2 = (double *)
	  (field_block[index_block]->field_values(1));
	result = result &&
	  unit_assert(test_field(v2,index_block,1,mx,my,mz,gx,gy,gz,ND3));

	// field 2
	field_descr->ghosts(2, &gx, &gy, &gz);
	long double * v3 = (long double *)
	  (field_block[index_block]->field_values(2));
	result = result &&
	  unit_assert(test_field(v3,index_block,2,mx,my,mz,gx,gy,gz,ND3));
 
      }
    }
  }
  return result;
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

  field_descr->insert_field("field_1");
  field_descr->insert_field("field_2");
  field_descr->insert_field("field_3");

  // initialize field precisions

  field_descr->set_precision(0, precision_single);
  field_descr->set_precision(1, precision_double);
  field_descr->set_precision(2, precision_quadruple);

  // initialize field ghost zone depths

  field_descr->set_ghosts(0, 1,1,1);
  field_descr->set_ghosts(1, 0,1,2);
  field_descr->set_ghosts(2, 1,0,3);


  int nx=4, ny=4, nz=4;
  FieldBlock * field_block[nx*ny*nz];

  int mx=3, my=5, mz=6;

  init_fields(field_descr,field_block,nx,ny,nz,mx,my,mz);

  //----------------------------------------------------------------------
  // Refresh all lower x-axis ghosts
  //----------------------------------------------------------------------

  for (int kz = 0; kz < nz; kz++) {
    for (int ky = 0; ky < ny; ky++) {
      for (int kx = 0; kx < (nx-1); kx++) { // kx < (nx-1)

	int index_lower = kx + nx * (ky + ny * kz);
	FieldBlock * field_lower = field_block[index_lower];
	FieldFaces * faces_lower = field_lower->field_faces();
	int index_upper = (kx + 1) + nx * (ky + ny * kz);
	FieldBlock * field_upper = field_block[index_upper];
	FieldFaces * faces_upper = field_upper->field_faces();

	// -1 == all fields
	faces_lower -> load(field_lower, -1, axis_x,  face_upper);
	faces_lower -> copy(faces_upper, -1, axis_x,  face_upper);
	faces_upper -> store(field_upper, -1, axis_x, face_lower);

      }
    }
  }

  unit_func("FieldFaces","load/copy/store");
  unit_assert(test_fields(field_descr,field_block,nx,ny,nz,mx,my,mz));

  //----------------------------------------------------------------------	
  // clean up
  //----------------------------------------------------------------------	

  for (int kz = 0; kz < nz; kz++) {
    for (int ky = 0; ky < ny; ky++) {
      for (int kx = 1; kx < nx; kx++) {
	int index_block = kx + nx * (ky + ny * kz);
	delete field_block[index_block];
	field_block[index_block] = 0;
      }
    }
  }

  delete field_descr;

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
