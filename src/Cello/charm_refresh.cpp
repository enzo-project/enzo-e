// See LICENSE_CELLO file for license and copyright information

/// @file     charm_refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with refreshing ghost zones

#ifdef CONFIG_USE_CHARM

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//----------------------------------------------------------------------

void CommBlock::p_refresh_begin() 
{
  TRACE("CommBlock::p_refresh_begin");
  simulation()->performance()->start_region(perf_refresh);
#ifdef TEMP_SKIP_REFRESH

  static int warning_displayed = false;

  if  ( ! warning_displayed ) {
    warning_displayed = true;
    WARNING("CommBlock::q_adapt_exit",
	    "CALLING p_output() INSTEAD OF REFRESH FOR IMAGE MESH CREATION");
  }

  //@@@@@@@@@@@@@@@@@@@@@@2
  double min_reduce[2];

  min_reduce[0] = 0.0;
  min_reduce[1] = 0.0;
  CkCallback callback (CkIndex_CommBlock::p_output(NULL), thisProxy);
  contribute( 2*sizeof(double), min_reduce, CkReduction::min_double, callback);

  //    thisProxy.p_output();
  //@@@@@@@@@@@@@@@@@@@@@@2
#else
  refresh(); 
#endif
}

//----------------------------------------------------------------------

void CommBlock::refresh ()
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  Boundary * boundary = simulation->problem()->boundary();

  TRACE ("CommBlock::refresh() setting QD q_refresh_end()");

  const Config * config = simulation->config();

  std::string refresh_type = config->field_refresh_type;

  if (refresh_type == "quiescence") {

    TRACE("Setting refresh quiescence detection callback");

    CkStartQD (CkCallback(CkIndex_CommBlock::q_refresh_end(),
			  thisProxy[thisIndex]));

    if (! is_leaf()) return;

  } 

  // which faces need to be refreshed?
  bool no_refresh[3][2];
  is_on_boundary (no_refresh);
  bool periodic = boundary->is_periodic();
  if (periodic) {
    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
	no_refresh[axis][face] = false;
      }
    }
  }

  // set face loop limits accordingly
  int ixm = no_refresh[0][0] ? 0 : -1;
  int iym = no_refresh[1][0] ? 0 : -1;
  int izm = no_refresh[2][0] ? 0 : -1;
  int ixp = no_refresh[0][1] ? 0 : 1;
  int iyp = no_refresh[1][1] ? 0 : 1;
  int izp = no_refresh[2][1] ? 0 : 1;

  // Include ghost zones in refresh?

  // Forest size needed for Index

  int n3[3];
  size_forest(&n3[0],&n3[1],&n3[2]);

  TRACE6 ("limits %d %d %d   %d %d %d",ixm,iym,izm,ixp,iyp,izp);

  int refresh_rank = config->field_refresh_rank;
  int rank         = simulation->dimension();
  bool refresh_type_counter = (refresh_type == "counter");
  if (refresh_type_counter) loop_refresh_.stop() = 0;

  // (icx0,icy0,icz0) are loop limits
  // e.g. 
  // if face (ix,iy,iz) = (0,0,1), then children are ([01], [01], 1)
  // if face (ix,iy,iz) = (0,-1,0), then children are ([01], 0, [01])
  // if face (ix,iy,iz) = (1,0,-1), then children are (1, [01], 0)
  for (int ix=ixm; ix<=ixp; ix++) {
    for (int iy=iym; iy<=iyp; iy++) {
      for (int iz=izm; iz<=izp; iz++) {
	int face_rank = rank - (abs(ix) + abs(iy) + abs(iz));
	if (face_rank >= refresh_rank) {
	  TRACE3("Calling refresh_same %d %d %d",ix,iy,iz);
	  Index index = index_.index_neighbor(ix,iy,iz,n3);

	  // update coarse neighbors
	  if (! is_neighbor(index)) {

	    if (refresh_type_counter) ++loop_refresh_.stop();
	    //	    refresh_coarse(index);

	  } else {

	    // update any fine neighbors...

	    int count_nibling = refresh_fine (index,ix,iy,iz,n3);

	    // ...update neighbor if none

	    if (count_nibling == 0 && index_ != index) {

	      if (refresh_type_counter) ++loop_refresh_.stop();

	      refresh_same(index,ix,iy,iz);

	    }
	  }
	}
      }
    }
  }
  TRACE1("CommBlock::refresh() face counter = %d",
	 loop_refresh_.index());
  if (refresh_type_counter) {
    ++loop_refresh_.stop();
    x_refresh_same (0,0,0,0,0);
  }
}

//----------------------------------------------------------------------

void CommBlock::refresh_coarse (Index index)
{
}

//----------------------------------------------------------------------

void CommBlock::refresh_same (Index index, 
			      int ix,  int iy,  int iz)
{

  
  TRACE3("ENTER refresh_same %d %d %d",ix,iy,iz);

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  FieldFace field_face (block_->field_block(),field_descr);

  field_face.set_face(ix,iy,iz);
	  
  int n; 
  char * array;
  field_face.load(&n, &array);

  thisProxy[index].x_refresh_same (n,array,-ix,-iy,-iz);
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_same (int n, char * buffer, int fx, int fy, int fz)
{
  TRACE3 ("CommBlock::x_refresh_same(%d %d %d)",fx,fy,fz);
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();

  if ( n != 0) {

    // n == 0 is the call from self to ensure x_refresh_same()
    // always gets called at least once

    FieldFace field_face(block_->field_block(), field_descr);

    field_face.set_face(fx,fy,fz);

    field_face.store (n, buffer);
  }

  const Config * config = simulation->config();
  std::string refresh_type = config->field_refresh_type;
  
  if (refresh_type == "counter") {
    if (loop_refresh_.done()) {
      q_refresh_end();
    }
  }
}

//----------------------------------------------------------------------

int CommBlock::refresh_fine 
(Index index, 
 int ix,  int iy,  int iz,
 int n3[3])
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();

  int rank = simulation->dimension();


  if (ix==0 && iy==0 && iz==0) return 0;

  int icxm = (ix == 0) ? 0 : (ix+1)/2;
  int icxp = (ix == 0) ? 1 : (ix+1)/2;
  int icym = (iy == 0) ? 0 : (iy+1)/2;
  int icyp = (iy == 0) ? 1 : (iy+1)/2;
  int iczm = (iz == 0) ? 0 : (iz+1)/2;
  int iczp = (iz == 0) ? 1 : (iz+1)/2;
  if (rank < 3) {iczm=0;  iczp=0; }
  if (rank < 2) {icym=0;  icyp=0; }

  TRACE9 ("NIBLING: Face (%d %d %d) child (%d:%d %d:%d %d:%d)",
	  ix,iy,iz, icxm,icxp,icym, icyp,iczm,iczp);


  int count = 0;
  // loop over niblings along face
  int ic3[3];
  for (ic3[0]=icxm; ic3[0]<=icxp; ic3[0]++) {
    for (ic3[1]=icym; ic3[1]<=icyp; ic3[1]++) {
      for (ic3[2]=iczm; ic3[2]<=iczp; ic3[2]++) {

	Index index_nibling = index.index_nibling(ix,iy,iz,ic3,n3);

	// if there exists a neighboring nibling

	if (is_nibling(index_nibling)) {

	  // increment counter

	  ++ count;

	  // pack face data

	  FieldFace field_face (block_->field_block(),field_descr);

	  field_face.set_face(ix,iy,iz);
	  
	  int n; 
	  char * array;
	  field_face.load(&n, &array);

	  // send face data

	  thisProxy[index].x_refresh_fine (n,array,-ix,-iy,-iz);

	}
      }
    }
  }
  return count;
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_fine (int n, char * buffer, int fx, int fy, int fz)
 {
   TRACE3 ("CommBlock::x_refresh_fine(%d %d %d)",fx,fy,fz);

   Simulation * simulation = proxy_simulation.ckLocalBranch();

   FieldDescr * field_descr = simulation->field_descr();
   const Config * config    = simulation->config();

   FieldFace field_face(block_->field_block(), field_descr);

   field_face.set_face(fx,fy,fz);
   
   field_face.store (n, buffer);

   std::string refresh_type = config->field_refresh_type;
  
   if (refresh_type == "counter") {
     if (loop_refresh_.done()) {
       q_refresh_end();
     }
   }
}

//----------------------------------------------------------------------

void CommBlock::q_refresh_end()
{
  simulation()->performance()->stop_region(perf_refresh);
  TRACE ("CommBlock::q_refresh_end() calling prepare()");
  prepare();
}

//----------------------------------------------------------------------

#endif /* CONFIG_USE_CHARM */
