// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_FieldFaces.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Unit tests for the FieldFaces class

#include <math.h>

#include "cello.hpp"

#include "error.hpp"
#include "test.hpp"
#include "field.hpp"

#include "parallel.def"
#include PARALLEL_CHARM_INCLUDE(test_FieldFaces.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init();

  //--------------------------------------------------
  // Initialize the field descriptor object fd
  //--------------------------------------------------

  FieldDescr fd;

  // insert fields

  int i1 = fd.insert_field("f1");
  int i2 = fd.insert_field("f2");
  int i3 = fd.insert_field("f3");

  // initialize field precisions

  fd.set_precision(i1, precision_single);
  fd.set_precision(i2, precision_double);
  fd.set_precision(i3, precision_quadruple);

  // initialize field ghost zone depths

  fd.set_ghosts(i1, 1,1,1);
  fd.set_ghosts(i2, 0,1,2);
  fd.set_ghosts(i3, 1,0,3);

  //--------------------------------------------------
  // Initialize the field blocks fb[0:6]
  //--------------------------------------------------

  // Seven field blocks in 7-point stencil pattern

  FieldBlock fb[7];

  float       *v1[7];
  double      *v2[7];
  long double *v3[7];

  const int XM=1,XP=2,YM=3,YP=4,ZM=5,ZP=6;

  for (int ib=0; ib<7; ib++) {

    // Set field blocks' field descriptors

    fb[ib].set_field_descr(&fd);

    // Set field blocks' dimensions

    const int nx=4, ny=5, nz=6;

    fb[ib].set_dimensions(nx,ny,nz);

    // Allocate field blocks including ghosts

    fb[ib].allocate_array();
    fb[ib].allocate_ghosts();

    // Get aliases for each field

    v1[ib]   = (float *)       fb[ib].field_values(i1);
    v2[ib]   = (double *)      fb[ib].field_values(i2);
    v3[ib]   = (long double *) fb[ib].field_values(i3);

    // Initialize fields
    // (each field initialized separately since each has different
    // ghost zone depth)

    int ndx,ndy,ndz;
    int gx,gy,gz;

    // field 1

    fd.ghosts(i1,&gx,&gy,&gz);

    ndx = nx + 2*gx;
    ndy = ny + 2*gy;
    ndz = nz + 2*gz;

    for (int iz=0; iz<ndz; iz++) {
      for (int iy=0; iy<ndy; iy++) {
	for (int ix=0; ix<ndx; ix++) {
	  int    i = ix + ndx*(iy + ndy*iz);
	  double v = i + 1000*(ib + 10*(1));
	  v1[ib][i] = 
	    (gx <= ix && ix < ndx-gx &&
	     gy <= iy && iy < ndy-gy &&
	     gz <= iz && iz < ndz-gz) ? v : -v;
	}
      }
    }

    // field 2

    fd.ghosts(i2,&gx,&gy,&gz);

    ndx = nx + 2*gx;
    ndy = ny + 2*gy;
    ndz = nz + 2*gz;

    for (int iz=0; iz<ndz; iz++) {
      for (int iy=0; iy<ndy; iy++) {
	for (int ix=0; ix<ndx; ix++) {
	  int i = ix + ndx*(iy + ndy*iz);
	  double v = i + 1000*(ib + 10*(2));
	  v2[ib][i] = 
	    (gx <= ix && ix < ndx-gx &&
	     gy <= iy && iy < ndy-gy &&
	     gz <= iz && iz < ndz-gz) ? v : -v;
	}
      }
    }

    // field 3

    fd.ghosts(i3,&gx,&gy,&gz);

    ndx = nx + 2*gx;
    ndy = ny + 2*gy;
    ndz = nz + 2*gz;

    for (int iz=0; iz<ndz; iz++) {
      for (int iy=0; iy<ndy; iy++) {
	for (int ix=0; ix<ndx; ix++) {
	  int i = ix + ndx*(iy + ndy*iz);
	  double v = i + 1000*(ib + 10*(3));
	  v3[ib][i] = 
	    (gx <= ix && ix < ndx-gx &&
	     gy <= iy && iy < ndy-gy &&
	     gz <= iz && iz < ndz-gz) ? v : -v;

	}
      }
    }
  }


  //----------------------------------------------------------------------
  // PASS if ghosts allocated
  //----------------------------------------------------------------------

  // Send faces in xp direction

  
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
