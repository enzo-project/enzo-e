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

  Performance * performance = simulation()->performance();
  if (! performance->is_region_active(perf_refresh)) {
    performance->start_region(perf_refresh);
  }
  //  simulation()->performance()->start_region(perf_refresh);

  refresh(); 
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

#ifdef CELLO_TRACE
    index_.print ("leaf",-1,2);
#endif

  } 

  // which faces need to be refreshed?
  bool no_refresh[3][2];
  is_on_boundary (no_refresh);
  TRACE6("is_on_boundary %d %d %d  %d %d %d",
	 no_refresh[0][0],no_refresh[1][0],no_refresh[2][0],
	 no_refresh[0][1],no_refresh[1][1],no_refresh[2][1]);

  bool periodic = boundary->is_periodic();
  if (periodic) {
    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
	no_refresh[axis][face] = false;
      }
    }
  }

  // set face loop limits accordingly
  int ifxm = no_refresh[0][0] ? 0 : -1;
  int ifym = no_refresh[1][0] ? 0 : -1;
  int ifzm = no_refresh[2][0] ? 0 : -1;
  int ifxp = no_refresh[0][1] ? 0 : 1;
  int ifyp = no_refresh[1][1] ? 0 : 1;
  int ifzp = no_refresh[2][1] ? 0 : 1;

  int rank         = simulation->dimension();
  if (rank < 2) ifym = ifyp = 0;
  if (rank < 3) ifzm = ifzp = 0;

  // Include ghost zones in refresh?

  // Forest size needed for Index

  int n3[3];
  size_forest(&n3[0],&n3[1],&n3[2]);

  TRACE6 ("limits %d %d %d   %d %d %d",ifxm,ifym,ifzm,ifxp,ifyp,ifzp);

  int refresh_rank = config->field_refresh_rank;
  bool refresh_type_counter = (refresh_type == "counter");
  if (refresh_type_counter) loop_refresh_.stop() = 0;

  // (icx0,icy0,icz0) are loop limits
  // e.g. 
  // if face (ifx,ify,ifz) = (0,0,1), then children are ([01], [01], 1)
  // if face (ifx,ify,ifz) = (0,-1,0), then children are ([01], 0, [01])
  // if face (ifx,ify,ifz) = (1,0,-1), then children are (1, [01], 0)
  for (int ifx=ifxm; ifx<=ifxp; ifx++) {
    for (int ify=ifym; ify<=ifyp; ify++) {
      for (int ifz=ifzm; ifz<=ifzp; ifz++) {
	int face_rank = rank - (abs(ifx) + abs(ify) + abs(ifz));
	if (face_rank >= refresh_rank) {
	  TRACE3("Calling refresh_same %d %d %d",ifx,ify,ifz);
	  Index index_neighbor = index_.index_neighbor(ifx,ify,ifz,n3);

	  // if no neighbor, update coarse

	  if (! is_neighbor(index_neighbor)) {

	    if (refresh_type_counter) ++loop_refresh_.stop();

	    refresh_coarse(index_neighbor.index_parent(),ifx,ify,ifz);

	  } else {

	    // if neighbor check for niblings...

	    int num_niblings = 0;

	    int icxm,icym,iczm;
	    int icxp,icyp,iczp;
	    loop_limits_nibling_(&icxm,&icym,&iczm,
				 &icxp,&icyp,&iczp,
				 ifx,ify,ifz);

	    int ic3[3];
	    for (ic3[0]=icxm; ic3[0]<=icxp; ic3[0]++) {
	      for (ic3[1]=icym; ic3[1]<=icyp; ic3[1]++) {
		for (ic3[2]=iczm; ic3[2]<=iczp; ic3[2]++) {
		  Index index_nibling = index_neighbor.index_child(ic3);
		  if (is_nibling(index_nibling)) {
		    ++num_niblings;
		    refresh_fine(index_nibling,ifx,ify,ifz,ic3);
		  }
		}
	      }
	    }

	    // otherwise update neighbor

	    if (num_niblings == 0 && index_ != index_neighbor) {

	      if (refresh_type_counter) ++loop_refresh_.stop();

	      refresh_same(index_neighbor,ifx,ify,ifz);

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

void CommBlock::refresh_same (Index index, 
			      int ifx,  int ify,  int ifz)
{

#ifdef CELLO_TRACE
  index_.print("refresh_same A",-1,2);
  index.print("refresh_same B",-1,2);
#endif

  TRACE3("ENTER refresh_same %d %d %d",ifx,ify,ifz);

  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();

  FieldFace field_face (field_block,field_descr);

  field_face.set_face(ifx,ify,ifz);
	  
  int n; 
  char * array;
  field_face.load(&n, &array);

  thisProxy[index].x_refresh_same (n,array,-ifx,-ify,-ifz);
}

//----------------------------------------------------------------------

void CommBlock::refresh_fine 
(Index index, int ifx,  int ify,  int ifz, int ic3[3])
{
  TRACE3("ENTER refresh_fine %d %d %d",ifx,ify,ifz);
  Simulation * simulation = proxy_simulation.ckLocalBranch();
  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();
  Prolong * prolong = simulation->problem()->prolong();

  // pack face data

  FieldFace field_face (field_block,field_descr);

  field_face.set_face(ifx,ify,ifz);
  field_face.set_prolong(prolong,ic3[0],ic3[1],ic3[2]);
	  
  int n; 
  char * array;
  field_face.load(&n, &array);

  // send face data

  TRACE6("Calling x_refresh_fine %d %d %d %d %d %d",
	 ifx,ify,ifz, ic3[0],ic3[1],ic3[2]);
  thisProxy[index].x_refresh_fine (n,array,
				   -ifx,-ify,-ifz,
				   ic3[0],ic3[1],ic3[2]);
	  
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

  const Config * config = simulation->config();
  std::string refresh_type = config->field_refresh_type;
  
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

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();
  Prolong * prolong = simulation->problem()->prolong();

  FieldFace field_face(field_block, field_descr);

  field_face.set_face(ifx,ify,ifz);
  field_face.set_prolong(prolong,icx,icy,icz);
   
  field_face.store (n, buffer);

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

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  FieldDescr * field_descr = simulation->field_descr();
  FieldBlock * field_block = block_->field_block();
  Restrict * restrict = simulation->problem()->restrict();

  FieldFace field_face(field_block, field_descr);

  field_face.set_face(ifx,ify,ifz);
  field_face.set_restrict(restrict,icx,icy,icz);
   
  field_face.store (n, buffer);

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
