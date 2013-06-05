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
  Performance * performance = simulation()->performance();
  if (! performance->is_region_active(perf_refresh)) {
    performance->start_region(perf_refresh);
  }

  refresh(); 
}

//----------------------------------------------------------------------

void CommBlock::refresh ()
{

  TRACE("BEGIN PHASE REFRESH");

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  const Config * config = simulation->config();

  std::string refresh_type = config->field_refresh_type;

  if (refresh_type == "quiescence") {

    TRACE("REFRESH Setting refresh quiescence detection callback");

    CkStartQD (CkCallback(CkIndex_CommBlock::q_refresh_end(),
			  thisProxy[thisIndex]));

    if (! is_leaf()) return;

  } 

  int ifxm,ifym,ifzm;
  int ifxp,ifyp,ifzp;
  loop_limits_refresh_(&ifxm,&ifym,&ifzm,&ifxp,&ifyp,&ifzp);

  // Forest size needed for Index

  int n3[3];
  size_forest(&n3[0],&n3[1],&n3[2]);

  int refresh_rank = config->field_refresh_rank;
  bool refresh_type_counter = (refresh_type == "counter");
  loop_refresh_.stop() = 0;

  int rank = simulation->dimension();

  for (int ifx=ifxm; ifx<=ifxp; ifx++) {
    for (int ify=ifym; ify<=ifyp; ify++) {
      for (int ifz=ifzm; ifz<=ifzp; ifz++) {
	int face_rank = rank - (abs(ifx) + abs(ify) + abs(ifz));
	if (refresh_rank <= face_rank && face_rank < rank) {

	  TRACE3("REFRESH Calling refresh_same %d %d %d",ifx,ify,ifz);
	  Index index_neighbor = index_.index_neighbor(ifx,ify,ifz,n3);

	  int if3[3] = {ifx,ify,ifz};

	  if (face_level(if3)==level_ - 1) {       // COARSE
	    refresh_coarse(index_neighbor.index_parent(),ifx,ify,ifz);

	  } else if (face_level(if3)==level_) {    // SAME
	    refresh_same(index_neighbor,ifx,ify,ifz);

	  } else if (face_level(if3)==level_+1) {  // FINE
	    
	    int icxm,icym,iczm;
	    int icxp,icyp,iczp;
	    loop_limits_nibling_(&icxm,&icym,&iczm,
				 &icxp,&icyp,&iczp,
				 ifx,ify,ifz);
	    int icx,icy,icz;
	    for (icx=icxm; icx<=icxp; icx++) {
	      for (icy=icym; icy<=icyp; icy++) {
		for (icz=iczm; icz<=iczp; icz++) {

		  int ic3[3] = {icx,icy,icz};
		  if (ifx == -1) ic3[0] = 1;
		  if (ify == -1) ic3[1] = 1;
		  if (ifz == -1) ic3[2] = 1;
		  if (ifx == +1) ic3[0] = 0;
		  if (ify == +1) ic3[1] = 0;
		  if (ifz == +1) ic3[2] = 0;

		  Index index_nibling = 
		    index_neighbor.index_child(ic3[0],ic3[1],ic3[2]);
		  
		  refresh_fine(index_nibling, ifx,ify,ifz, ic3[0],ic3[1],ic3[2]);
		}
	      }
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

void CommBlock::refresh_coarse 
(
 Index index, 
 int ifx,  int ify,  int ifz)
{
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

void CommBlock::x_refresh_coarse (int n, char * buffer, 
				int ifx, int ify, int ifz,
				int icx, int icy, int icz)
{
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

void CommBlock::refresh_same (Index index, 
			      int ifx,  int ify,  int ifz)
{
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

void CommBlock::x_refresh_same (int n, char * buffer, int ifx, int ify, int ifz)
{
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

void CommBlock::refresh_fine 
(Index index, 
 int ifx, int ify, int ifz, 
 int icx, int icy, int icz)
{
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

  ++loop_refresh_.stop();

  thisProxy[index].x_refresh_fine (n,array,
				   -ifx,-ify,-ifz,
				   icx,icy,icz);
	  
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_fine (int n, char * buffer, 
				int ifx, int ify, int ifz,
				int icx, int icy, int icz)
{
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

void CommBlock::q_refresh_end()
{
  Performance * performance = simulation()->performance();
  if (performance->is_region_active(perf_refresh))
    performance->stop_region(perf_refresh);

  char buffer[80];
  sprintf (buffer,"refresh-end-%d",cycle_);
  block()->field_block()->print(simulation()->field_descr(),buffer,true);
  TRACE("END   PHASE REFRESH");
  
  if (next_phase_ == phase_output) {
    TRACE("refresh calling output");
    prepare();
  } else if (next_phase_ == phase_adapt) {
    TRACE("refresh calling adapt");
    p_adapt_begin();
  } else {
    ERROR1 ("CommBlock::q_refresh_end()",
	    "Unknown next_phase %d",
	    next_phase_);
  }
}

//----------------------------------------------------------------------

#endif /* CONFIG_USE_CHARM */
