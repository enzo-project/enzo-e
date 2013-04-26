// See LICENSE_CELLO file for license and copyright information

/// @file     charm_refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with initialization

#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

// void SimulationCharm::p_refresh()
// { refresh(); };

//----------------------------------------------------------------------

void SimulationCharm::refresh()
{
  TRACE("SimulationCharm::refresh");
}

//----------------------------------------------------------------------

void CommBlock::p_refresh() 
{
  TRACE("CommBlock::p_refresh");
  refresh(); 
}

//----------------------------------------------------------------------

void CommBlock::refresh ()
{
  TRACE ("CommBlock::refresh()");

  bool is_boundary[3][2];

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  Boundary * boundary = simulation->problem()->boundary();
  FieldDescr * field_descr = simulation->field_descr();
  
  bool periodic = boundary->is_periodic();

  CProxy_CommBlock block_array = thisProxy;

  //--------------------------------------------------
  // Refresh
  //--------------------------------------------------

  int ibx,iby,ibz;

  thisIndex.array(&ibx,&iby,&ibz);

  int nbx = size_[0];
  int nby = size_[1];
  int nbz = size_[2];
  
  bool fx3[3],fy3[3],fz3[3];
  determine_boundary_
    (is_boundary,&fx3[0],&fx3[2],&fy3[0],&fy3[2],&fz3[0],&fz3[2]);
  fx3[1]=true;
  fy3[1]=true;
  fz3[1]=true;
  fx3[0] = fx3[0] && (periodic || ! is_boundary[axis_x][face_lower]);
  fx3[2] = fx3[2] && (periodic || ! is_boundary[axis_x][face_upper]);
  fy3[0] = fy3[0] && (periodic || ! is_boundary[axis_y][face_lower]);
  fy3[2] = fy3[2] && (periodic || ! is_boundary[axis_y][face_upper]);
  fz3[0] = fz3[0] && (periodic || ! is_boundary[axis_z][face_lower]);
  fz3[2] = fz3[2] && (periodic || ! is_boundary[axis_z][face_upper]);
  int ix3[3],iy3[3],iz3[3];
  ix3[0] = (ibx - 1 + nbx) % nbx;
  iy3[0] = (iby - 1 + nby) % nby;
  iz3[0] = (ibz - 1 + nbz) % nbz;
  ix3[1] = ibx;
  iy3[1] = iby;
  iz3[1] = ibz;
  ix3[2] = (ibx + 1) % nbx;
  iy3[2] = (iby + 1) % nby;
  iz3[2] = (ibz + 1) % nbz;
  
  // Refresh face ghost zones

  bool gx,gy,gz;
  gx = false;
  gy = false;
  gz = false;

  int fxl = 1;
  int fyl = (nby==1 && ! periodic) ? 0 : 1;
  int fzl = (nbz==1 && ! periodic) ? 0 : 1;

  for (int fx=-fxl; fx<=fxl; fx++) {
    for (int fy=-fyl; fy<=fyl; fy++) {
      for (int fz=-fzl; fz<=fzl; fz++) {
	int sum = abs(fx)+abs(fy)+abs(fz);
	if ((fx3[fx+1] && fy3[fy+1] && fz3[fz+1]) &&
	    ((sum==1 && field_descr->refresh_face(2)) ||
	     (sum==2 && field_descr->refresh_face(1)) ||
	     (sum==3 && field_descr->refresh_face(0)))) {

	  FieldFace field_face (block_->field_block(),field_descr);

	  field_face.set_face(fx,fy,fz);
	  field_face.set_ghost(gx,gy,gz);
	  
	  DEBUG9("index %d %d %d  %d %d %d  %d %d %d",
		 index_[0],index_[1],index_[2],
		 ix3[fx+1],iy3[fy+1],iz3[fz+1],
		 fx,fy,fz);

	  int n; 
	  char * array;
	  field_face.load(&n, &array);

	  Index index;

	  index.set_array(ix3[fx+1],iy3[fy+1],iz3[fz+1]);
	  index.set_level(0);
	  index.clean();

	  thisProxy[index].x_refresh (n,array,-fx,-fy,-fz);
	}
      }
    }
  }

  // NOTE: x_refresh() calls compute, but if no incoming faces
  // it will never get called.  So every block also calls
  // x_refresh() itself with a null array

  x_refresh (0,0,0, 0, 0);

}

//----------------------------------------------------------------------

void CommBlock::x_refresh (int n, char * buffer, int fx, int fy, int fz)
{

  TRACE ("CommBlock::x_refresh()");
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();

  if ( n != 0) {

    // n == 0 is the call from self to ensure x_refresh()
    // always gets called at least once

    bool gx,gy,gz;
    gx = false;
    gy = false;
    gz = false;

    FieldFace field_face(block_->field_block(), field_descr);

    field_face.set_face(fx,fy,fz);
    field_face.set_ghost(gx,gy,gz);

    field_face.store (n, buffer);
  }

  //--------------------------------------------------
  // Compute
  //--------------------------------------------------

  if (sync_refresh_.done()) {
    TRACE ("CommBlock::x_refresh() calling prepare()");
    prepare();
  }
}

//----------------------------------------------------------------------

int CommBlock::count_refresh_()
{
  // Count incoming faces for face 

  int nx,ny,nz;
  block_->field_block()->size (&nx,&ny,&nz);

  // Determine axes that may be neighbors

  bool fxm = nx > 1;
  bool fxp = nx > 1;
  bool fym = ny > 1;
  bool fyp = ny > 1;
  bool fzm = nz > 1;
  bool fzp = nz > 1;

  // Adjust for boundary faces

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  bool periodic = simulation->problem()->boundary()->is_periodic();

  bool is_boundary[3][2];
  is_on_boundary (is_boundary);

  fxm = fxm && (periodic || ! is_boundary[axis_x][face_lower]);
  fxp = fxp && (periodic || ! is_boundary[axis_x][face_upper]);
  fym = fym && (periodic || ! is_boundary[axis_y][face_lower]);
  fyp = fyp && (periodic || ! is_boundary[axis_y][face_upper]);
  fzm = fzm && (periodic || ! is_boundary[axis_z][face_lower]);
  fzp = fzp && (periodic || ! is_boundary[axis_z][face_upper]);

  // count self

  int count = 1;

  // count faces

  FieldDescr * field_descr = simulation->field_descr();

  if (field_descr->refresh_face(2)) {
    if ( fxm ) ++count;
    if ( fxp ) ++count;
    if ( fym ) ++count;
    if ( fyp ) ++count;
    if ( fzm ) ++count;
    if ( fzp ) ++count;
  }

  // count edges

  if (field_descr->refresh_face(1)) {
    if ( fxm && fym ) ++count;
    if ( fxm && fyp ) ++count;
    if ( fxp && fym ) ++count;
    if ( fxp && fyp ) ++count;
    if ( fym && fzm ) ++count;
    if ( fym && fzp ) ++count;
    if ( fyp && fzm ) ++count;
    if ( fyp && fzp ) ++count;
    if ( fzm && fxm ) ++count;
    if ( fzm && fxp ) ++count;
    if ( fzp && fxm ) ++count;
    if ( fzp && fxp ) ++count;
  }

  // count corners

  if (field_descr->refresh_face(0)) {
    if ( fxm && fym && fzm ) ++count;
    if ( fxm && fym && fzp ) ++count;
    if ( fxm && fyp && fzm ) ++count;
    if ( fxm && fyp && fzp ) ++count;
    if ( fxp && fym && fzm ) ++count;
    if ( fxp && fym && fzp ) ++count;
    if ( fxp && fyp && fzm ) ++count;
    if ( fxp && fyp && fzp ) ++count;
  }

  return count;
}

#endif /* CONFIG_USE_CHARM */
