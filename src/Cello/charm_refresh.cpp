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
#ifdef TEMP_SKIP_REFRESH
  WARNING("CommBlock::q_adapt_exit",
	  "CALLING p_output() INSTEAD OF REFRESH FOR IMAGE MESH CREATION");

  //@@@@@@@@@@@@@@@@@@@@@@2
  double min_reduce[2];

  min_reduce[0] = 0.0;
  min_reduce[1] = 0.0;
  CkCallback callback (CkIndex_CommBlock::p_output(NULL), thisProxy);
  contribute( 2*sizeof(double), min_reduce, CkReduction::min_double, callback);

  //    thisProxy.p_output();
  //@@@@@@@@@@@@@@@@@@@@@@2
#else
  //  refresh(); 
#endif

  TRACE("CommBlock::p_refresh_begin");

  Performance * performance = simulation()->performance();
  if (! performance->is_region_active(perf_refresh)) {
    performance->start_region(perf_refresh);
  }

  refresh(); 
}

//----------------------------------------------------------------------

void CommBlock::refresh ()
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  TRACE ("CommBlock::refresh() setting QD q_refresh_end()");

  const Config * config = simulation->config();

  std::string refresh_type = config->field_refresh_type;

  if (refresh_type == "quiescence") {

    TRACE("Setting refresh quiescence detection callback");

    CkStartQD (CkCallback(CkIndex_CommBlock::q_refresh_end(),
			  thisProxy[thisIndex]));

    if (! is_leaf()) return;

#ifdef CELLO_TRACE
    index_.print ("leaf",-1,2);
#endif

  } 

  int ifxm,ifym,ifzm;
  int ifxp,ifyp,ifzp;
  loop_limits_refresh_(&ifxm,&ifym,&ifzm,&ifxp,&ifyp,&ifzp);

  // Include ghost zones in refresh?

  // Forest size needed for Index

  int n3[3];
  size_forest(&n3[0],&n3[1],&n3[2]);

  TRACE6 ("limits %d %d %d   %d %d %d",ifxm,ifym,ifzm,ifxp,ifyp,ifzp);

  int refresh_rank = config->field_refresh_rank;
  bool refresh_type_counter = (refresh_type == "counter");
  loop_refresh_.stop() = 0;

  int rank = simulation->dimension();

  for (int ifx=ifxm; ifx<=ifxp; ifx++) {
    for (int ify=ifym; ify<=ifyp; ify++) {
      for (int ifz=ifzm; ifz<=ifzp; ifz++) {
	int face_rank = rank - (abs(ifx) + abs(ify) + abs(ifz));
	if (face_rank >= refresh_rank) {
	  TRACE3("Calling refresh_same %d %d %d",ifx,ify,ifz);
	  Index index_neighbor = index_.index_neighbor(ifx,ify,ifz,n3);

	  // if no neighbor in level, refresh coarse neighbor

	  int if3[3] = {ifx,ify,ifz};

	  if (! is_neighbor(index_neighbor,if3)) {

#ifdef CELLO_TRACE
	    index_.print("not neighbor A");
	    index_neighbor.print("not neighbor A");
#endif
	    refresh_coarse(index_neighbor.index_parent(),ifx,ify,ifz);

	  } else {

	    // else check for adjacent fine neighbors

	    int icxm,icym,iczm;
	    int icxp,icyp,iczp;
	    loop_limits_nibling_(&icxm,&icym,&iczm,
				 &icxp,&icyp,&iczp,
				 ifx,ify,ifz);

	    int num_niblings = 0;
	    int icx,icy,icz;
	    for (icx=icxm; icx<=icxp; icx++) {
	      for (icy=icym; icy<=icyp; icy++) {
		for (icz=iczm; icz<=iczp; icz++) {
		  Index index_nibling = 
		    index_neighbor.index_child(icx,icy,icz);
		  if (is_nibling(index_nibling,if3)) {
		    ++num_niblings;
		    refresh_fine(index_nibling, ifx,ify,ifz, icx,icy,icz);
		  }
		}
	      }
	    }

	    // otherwise update neighbor

	    if (num_niblings == 0 && index_ != index_neighbor) {

	      refresh_same(index_neighbor,ifx,ify,ifz);

	    }
	  }
	}
      }
    }
  }

  if (refresh_type_counter) {

    // Prevent hang if single-CommBlock simulation
    ++loop_refresh_.stop();

    x_refresh_same (0,0,0,0,0);

  }
}

//----------------------------------------------------------------------

void CommBlock::loop_limits_refresh_
(int * ifxm, int * ifym, int * ifzm, int * ifxp, int * ifyp, int * ifzp)
  const throw()
{

  Boundary * boundary = simulation()->problem()->boundary();

  // which faces need to be refreshed?
  bool on_boundary[3][2];
  is_on_boundary (on_boundary);
  TRACE6("is_on_boundary %d %d %d  %d %d %d",
	 on_boundary[0][0],on_boundary[1][0],on_boundary[2][0],
	 on_boundary[0][1],on_boundary[1][1],on_boundary[2][1]);

  bool periodic = boundary->is_periodic();
  if (periodic) {
    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
	on_boundary[axis][face] = false;
      }
    }
  }

  // set face loop limits accordingly
  (*ifxm) = on_boundary[0][0] ? 0 : -1;
  (*ifym) = on_boundary[1][0] ? 0 : -1;
  (*ifzm) = on_boundary[2][0] ? 0 : -1;
  (*ifxp) = on_boundary[0][1] ? 0 : 1;
  (*ifyp) = on_boundary[1][1] ? 0 : 1;
  (*ifzp) = on_boundary[2][1] ? 0 : 1;

  int rank = simulation()->dimension();
  if (rank < 2) (*ifym) = (*ifyp) = 0;
  if (rank < 3) (*ifzm) = (*ifzp) = 0;

}

//----------------------------------------------------------------------

void CommBlock::refresh_same (Index index, 
			      int ifx,  int ify,  int ifz)
{

#ifdef CELLO_TRACE
  index_.print("refresh_same A",-1,2);
  index.print("refresh_same B",-1,2);
#endif

  TRACE3("ENTER refresh_same %d %d %d",ifx,ify,ifz);

  // <duplicated code: refactor me!>
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();

  FieldFace field_face (field_block,field_descr);

  field_face.set_face(ifx,ify,ifz);
	  
  int n; 
  char * array;
  field_face.load(&n, &array);
  // </duplicated code>

  ++loop_refresh_.stop();

  thisProxy[index].x_refresh_same (n,array,-ifx,-ify,-ifz);
}

//----------------------------------------------------------------------

void CommBlock::refresh_fine 
(Index index, 
 int ifx, int ify, int ifz, 
 int icx, int icy, int icz)
{
  TRACE3("ENTER refresh_fine %d %d %d",ifx,ify,ifz);

  // <duplicated code: refactor me!>
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();
  Prolong * prolong = simulation->problem()->prolong();

  // pack face data

  FieldFace field_face (field_block,field_descr);

  field_face.set_face(ifx,ify,ifz);
  field_face.set_prolong(prolong,icx,icy,icz);
	  
  int n; 
  char * array;
  field_face.load(&n, &array);
  // </duplicated code>

  // send face data

  TRACE6("Calling x_refresh_fine %d %d %d %d %d %d",
	 ifx,ify,ifz, icx,icy,icz);

  ++loop_refresh_.stop();

  thisProxy[index].x_refresh_fine (n,array,
				   -ifx,-ify,-ifz,
				   icx,icy,icz);
	  
}

//----------------------------------------------------------------------

void CommBlock::loop_limits_nibling_ 
(int *icxm, int *icym, int *iczm,
 int *icxp, int *icyp, int *iczp,
 int ifx, int ify, int ifz) const throw()
{
  int rank = simulation()->dimension();

  (*icxm) = (ifx == 0) ? 0 : (ifx+1)/2;
  (*icxp) = (ifx == 0) ? 1 : (ifx+1)/2;
  (*icym) = (ify == 0) ? 0 : (ify+1)/2;
  (*icyp) = (ify == 0) ? 1 : (ify+1)/2;
  (*iczm) = (ifz == 0) ? 0 : (ifz+1)/2;
  (*iczp) = (ifz == 0) ? 1 : (ifz+1)/2;
  if (rank < 2) (*icym) = (*icyp) = 0;
  if (rank < 3) (*iczm) = (*iczp) = 0;
}

//----------------------------------------------------------------------

void CommBlock::refresh_coarse 
(
 Index index, 
 int ifx,  int ify,  int ifz)
{

  TRACE3("ENTER refresh_coarse %d %d %d",ifx,ify,ifz);
#ifdef CELLO_TRACE
  index_.print("refresh_coarse A",-1,2);
  index.print ("refresh_coarse B",-1,2);
#endif

  // <duplicated code: refactor me!>
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();
  Restrict * restrict = simulation->problem()->restrict();

  int level = index_.level();
  int icx,icy,icz;
  index_.child(level,&icx,&icy,&icz);

  FieldFace field_face (field_block,field_descr);

  field_face.set_restrict(restrict,icx,icy,icz);
  field_face.set_face(ifx,ify,ifz);
	  
  int n; 
  char * array;
  field_face.load(&n, &array);
  // </duplicated code>

  ++loop_refresh_.stop();

  thisProxy[index].x_refresh_coarse (n,array,-ifx,-ify,-ifz,
				     icx,icy,icz);

}

//----------------------------------------------------------------------

void CommBlock::x_refresh_same (int n, char * buffer, int ifx, int ify, int ifz)
{
  TRACE3 ("CommBlock::x_refresh_same(%d %d %d)",ifx,ify,ifz);
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();

  if ( n != 0) {

    // n == 0 is the call from self to ensure x_refresh_same()
    // always gets called at least once

    FieldFace field_face(field_block, field_descr);

    field_face.set_face(ifx,ify,ifz);

    field_face.store (n, buffer);
  }

  std::string refresh_type = simulation->config()->field_refresh_type;
  
  if (refresh_type == "counter") {
    if (loop_refresh_.done()) {
      q_refresh_end();
    }
  }
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_fine (int n, char * buffer, 
				int ifx, int ify, int ifz,
				int icx, int icy, int icz)
{
#ifdef CELLO_TRACE
  index_.print("x_refresh_fine",-1,2);
#endif
  TRACE6 ("CommBlock::x_refresh_fine(face %d %d %d  child  %d %d %d)",
	  ifx,ify,ifz,icx,icy,icz);

  // <duplicated code: refactor me!>
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();

  FieldFace field_face(field_block, field_descr);

  field_face.set_face(ifx,ify,ifz);

  Prolong *    prolong = simulation->problem()->prolong();
  field_face.set_prolong(prolong,icx,icy,icz);
   
  field_face.store (n, buffer);
  // </duplicated code>

  std::string refresh_type = simulation->config()->field_refresh_type;
  
  if (refresh_type == "counter") {
    if (loop_refresh_.done()) {
      q_refresh_end();
    }
  }
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_coarse (int n, char * buffer, 
				int ifx, int ify, int ifz,
				int icx, int icy, int icz)
{
#ifdef CELLO_TRACE
  index_.print("x_refresh_coarse",-1,2);
#endif
  TRACE6 ("CommBlock::x_refresh_coarse(face %d %d %d  child  %d %d %d)",
	  ifx,ify,ifz,icx,icy,icz);

  // <duplicated code: refactor me!>
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();

  FieldFace field_face(field_block, field_descr);

  field_face.set_face(ifx,ify,ifz);

  Restrict *   restrict = simulation->problem()->restrict();
  field_face.set_restrict(restrict,icx,icy,icz);
   
  field_face.store (n, buffer);
  // </duplicated code>

  std::string refresh_type = simulation->config()->field_refresh_type;
  
  if (refresh_type == "counter") {
    if (loop_refresh_.done()) {
      q_refresh_end();
    }
  }
}

//----------------------------------------------------------------------

void CommBlock::q_refresh_end()
{
  Performance * performance = simulation()->performance();
  if (performance->is_region_active(perf_refresh))
    performance->stop_region(perf_refresh);
  //  simulation()->performance()->stop_region(perf_refresh);

  TRACE ("CommBlock::q_refresh_end() calling prepare()");
  prepare();
}

//----------------------------------------------------------------------

#endif /* CONFIG_USE_CHARM */
