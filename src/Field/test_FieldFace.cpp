// See LICENSE_CELLO file for license and copyright information

/// @file     test_FieldFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the FieldFace class

#include "main.hpp" 
#include "test.hpp"

#include "field.hpp"

//----------------------------------------------------------------------

bool is_ghost(int ix,int iy,int iz,
	      int gx,int gy,int gz,
	      int mdx,int mdy,int mdz)
{
  return ( ! ((gx <= ix) && (ix < (mdx - gx)) && 
	      (gy <= iy) && (iy < (mdy - gy)) && 
	      (gz <= iz) && (iz < (mdz - gz))) );
}

//----------------------------------------------------------------------

double test_value (int ix, int iy, int iz,
		   int gx,int gy,int gz,
		   int mdx, int mdy, int mdz,
		   int ibx, int iby, int ibz,
		   int nbx, int nby, int nbz,
		   int index_field, 
		   int MD3, int ND3)
{
  int mx = (mdx-2*gx);
  int my = (mdy-2*gy);
  int mz = (mdz-2*gz);
  int NX = nbx*mx;
  int NY = nby*my;
  int NZ = nbz*mz;
  double xg = (ibx*mx + ix + NX) % NX;
  double yg = (iby*my + iy + NY) % NY;
  double zg = (ibz*mz + iz + NZ) % NZ;
  double value = xg + 3*yg + 5*zg;
  return value;
}

//----------------------------------------------------------------------

double init_value (int ix, int iy, int iz,
		   int gx,int gy,int gz,
		   int mdx, int mdy, int mdz,
		   int ibx, int iby, int ibz,
		   int nbx, int nby, int nbz,
		   int index_field, 
		   int MD3, int ND3)
{
  double value = test_value (ix,iy,iz,gx,gy,gz,mdx,mdy,mdz,ibx,iby,ibz,nbx,nby,nbz,index_field,MD3,ND3);
  return is_ghost(ix,iy,iz,gx,gy,gz,mdx,mdy,mdz) ? -value : value;
}

//----------------------------------------------------------------------

template<class T>
void init_field(T * field_values,
		int ibx, int iby, int ibz,
		int nbx, int nby, int nbz,
		int index_field, 
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
				     ibx,iby,ibz,
				     nbx,nby,nbz,
				     index_field,
				     MD3,ND3);
      }
    }
  }
}

//----------------------------------------------------------------------

template<class T>
bool test_field(T * field_values,
		int ibx, int iby, int ibz,
		int nbx, int nby, int nbz,
		int index_field, 
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
	double value = test_value(ix,iy,iz,
				  gx,gy,gz,
				  mdx,mdy,mdz,
				  ibx,iby,ibz,
				  nbx,nby,nbz,
				  index_field,
				  MD3,ND3);
	//	printf ("%p test\n",field_values+i);
	if (field_values[i] != value) {
	  result = false;
	  PARALLEL_PRINTF ("mismatch block   (%d %d %d)\n",            ibx,iby,ibz);
	  PARALLEL_PRINTF ("mismatch field    %d\n",                   index_field);
	  PARALLEL_PRINTF ("mismatch index   (%d %d %d) %d\n",         ix,iy,iz,i);
	  PARALLEL_PRINTF ("mismatch size    (%d %d %d)\n",            mdx,mdy,mdz);
	  PARALLEL_PRINTF ("mismatch actual   %g\n",  (double) field_values[i]);
	  PARALLEL_PRINTF ("mismatch expected %g\n",  (double) value);
	}
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
 int nbx,int nby, int nbz,
 int mx, int my, int mz)
{
  //--------------------------------------------------
  // Initialize the field blocks in 4 x 4 x 4 array
  //--------------------------------------------------

  const int ND3 = nbx*nby*nbz;

  for (int ibz = 0; ibz < nbz; ibz++) {
    for (int iby = 0; iby < nby; iby++) {
      for (int ibx = 0; ibx < nbx; ibx++) {

	int index_block = ibx + nbx * (iby + nby * ibz);

	// Create the FaceBlock object

	field_block[index_block] = new FieldBlock (mx, my, mz);

	FieldBlock * block = field_block[index_block];

	// Allocate field blocks including ghosts

	block->allocate_array(field_descr);
	block->allocate_ghosts(field_descr);

	// Initialize fields

	int gx, gy, gz;

	// field 0
	field_descr->ghosts(0, &gx, &gy, &gz);
	float * v1 = (float *) (block->field_values(0));
	init_field(v1,ibx,iby,ibz,nbx,nby,nbz,0,mx,my,mz,gx,gy,gz,ND3);

	// field 1
	field_descr->ghosts(1, &gx, &gy, &gz);
	double * v2 = (double *) (block->field_values(1));
	init_field(v2,ibx,iby,ibz,nbx,nby,nbz,1,mx,my,mz,gx,gy,gz,ND3);

	// field 2
	field_descr->ghosts(2, &gx, &gy, &gz);
	long double * v3 = (long double *) (block->field_values(2));
	init_field(v3,ibx,iby,ibz,nbx,nby,nbz,2,mx,my,mz,gx,gy,gz,ND3);
 
      }
    }
  }
}

//----------------------------------------------------------------------

bool test_fields
(
 FieldDescr * field_descr,
 FieldBlock * field_block[],
 int nbx,int nby, int nbz,
 int mx, int my, int mz)
{
  const int ND3 = nbx*nby*nbz;

  bool result = true;

  for (int ibz = 0; ibz < nbz; ibz++) {
    for (int iby = 0; iby < nby; iby++) {
      for (int ibx = 0; ibx < nbx; ibx++) {

	int index_block = ibx + nbx * (iby + nby * ibz);

	FieldBlock * block = field_block[index_block];

	int gx, gy, gz;

	bool test_result;

	// field 0
	field_descr->ghosts(0, &gx, &gy, &gz);
	float * v1 = (float *) (block->field_values(0));
	test_result = test_field(v1,ibx,iby,ibz,nbx,nby,nbz,0,mx,my,mz,gx,gy,gz,ND3);
	unit_assert(test_result); 
	result = result && test_result;

	// field 1
	field_descr->ghosts(1, &gx, &gy, &gz);
	double * v2 = (double *) (block->field_values(1));
	test_result = test_field(v2,ibx,iby,ibz,nbx,nby,nbz,1,mx,my,mz,gx,gy,gz,ND3);
	unit_assert(test_result); 
	result = result && test_result;

	// field 2
	field_descr->ghosts(2, &gx, &gy, &gz);
	long double * v3 = (long double *) (block->field_values(2));
	test_result = test_field(v3,ibx,iby,ibz,nbx,nby,nbz,2,mx,my,mz,gx,gy,gz,ND3);
	unit_assert(test_result); 
	result = result && test_result;
 
      }
    }
  }
  return result;
}

//======================================================================
PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("FieldFace");

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
  field_descr->set_ghosts(1, 1,2,3);
  field_descr->set_ghosts(2, 3,2,1);


  int nbx=2, nby=3, nbz=4;
  FieldBlock * field_block[nbx*nby*nbz];

  int mx=5, my=6, mz=7;

  init_fields(field_descr,field_block,nbx,nby,nbz,mx,my,mz);

  //----------------------------------------------------------------------
  // Refresh ghosts
  //----------------------------------------------------------------------

  FieldFace face;

  for (int ia = 0; ia<3; ++ia) {

    for (int ibz = 0; ibz < nbz; ibz++) {
      for (int iby = 0; iby < nby; iby++) {
	for (int ibx = 0; ibx < nbx; ibx++) {
	  axis_enum axis = (axis_enum)(ia);
	  
	  int index_lower = ibx + nbx * (iby + nby * ibz);
	  FieldBlock * field_lower = field_block[index_lower];

	  int index_upper = 0;
	  if (axis==0) index_upper = ((ibx+1)%nbx) + nbx * (  iby        + nby * ibz);
	  if (axis==1) index_upper =   ibx        + nbx * (((iby+1)%nby) + nby * ibz);
	  if (axis==2) index_upper =   ibx        + nbx * (  iby        + nby * ((ibz+1)%nbz));

	  FieldBlock * field_upper = field_block[index_upper];

	  face.allocate   (field_descr, field_lower, axis);
	  face.load       (field_descr, field_lower, axis, face_upper);
	  face.store      (field_descr, field_upper, axis, face_lower);

	  face.load       (field_descr, field_upper, axis, face_lower);
	  face.store      (field_descr, field_lower, axis, face_upper);
	  face.deallocate ();

	}
      }
    }
  }

  unit_func("load/copy/store");
  unit_assert(test_fields(field_descr,field_block,nbx,nby,nbz,mx,my,mz));

  //----------------------------------------------------------------------	
  // clean up
  //----------------------------------------------------------------------	

  for (int ibz = 0; ibz < nbz; ibz++) {
    for (int iby = 0; iby < nby; iby++) {
      for (int ibx = 1; ibx < nbx; ibx++) {
	int index_block = ibx + nbx * (iby + nby * ibz);
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

