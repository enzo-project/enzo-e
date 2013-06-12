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

  int rank = simulation->dimension();

  // Forest size needed for Index

  int n3[3];
  size_forest(&n3[0],&n3[1],&n3[2]);

  int refresh_rank = config->field_refresh_rank;
  bool refresh_type_counter = (refresh_type == "counter");
  loop_refresh_.stop() = 0;

  ItFace it_face(rank,refresh_rank);
  int if3[3];
  while (it_face.next(if3)) {

    Index index_neighbor = index_.index_neighbor(if3[0],if3[1],if3[2],n3);

    if (face_level(if3)==level_ - 1) {       // COARSE

      int ic3[3];
      index_.child(level_,ic3+0,ic3+1,ic3+2);
      int ip3[3];
      parent_face_(ip3,if3,ic3);

      refresh_coarse(index_neighbor.index_parent(),ip3[0],ip3[1],ip3[2]);

    } else if (face_level(if3)==level_) {    // SAME

      refresh_same(index_neighbor,if3[0],if3[1],if3[2]);

    } else if (face_level(if3)==level_+1) {  // FINE
	    
      int ic3m[3];
      int ic3p[3];
      loop_limits_nibling_(ic3m,ic3m+1,ic3m+2,
			   ic3p,ic3p+1,ic3p+2,
			   if3[0],if3[1],if3[2]);
      int ic3[3];
      for (ic3[0]=ic3m[0]; ic3[0]<=ic3p[0]; ic3[0]++) {
	for (ic3[1]=ic3m[1]; ic3[1]<=ic3p[1]; ic3[1]++) {
	  for (ic3[2]=ic3m[2]; ic3[2]<=ic3p[2]; ic3[2]++) {

	    int jc3[3];
	    facing_child_ (jc3,ic3,if3);

	    Index index_nibling = 
	      index_neighbor.index_child(jc3[0],jc3[1],jc3[2]);
		  
	    refresh_fine(index_nibling, 
			 if3[0],if3[1],if3[2], 
			 ic3[0],ic3[1],ic3[2]);
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

  int level = index_.level();
  int ic3[3];
  index_.child(level,ic3,ic3+1,ic3+2);

   int n; 
   char * array;
   FieldFace * field_face = load_face_(&n,&array,
				       ifx,ify,ifz,
				       ic3[0],ic3[1],ic3[2],
				       false,false,false,
				       op_array_restrict);

  ++loop_refresh_.stop();

  thisProxy[index].x_refresh_coarse (n,array,-ifx,-ify,-ifz,
				     ic3[0],ic3[1],ic3[2]);
  delete field_face;
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_coarse (int n, char * buffer, 
				int ifx, int ify, int ifz,
				int icx, int icy, int icz)
{
  store_face_(n,buffer,
	      ifx,ify,ifz,
	      icx,icy,icz,
	      false,false,false,
	      op_array_restrict);

  std::string refresh_type = simulation()->config()->field_refresh_type;
  
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
  int n; 
  char * array;

  FieldFace * field_face = load_face_ (&n,&array,
				       ifx,ify,ifz,
				       0,0,0,
				       false,false,false,
				       op_array_copy);

  ++loop_refresh_.stop();

  thisProxy[index].x_refresh_same (n,array,-ifx,-ify,-ifz);

  delete field_face;
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_same (int n, char * buffer, int ifx, int ify, int ifz)
{
  Simulation * simulation = proxy_simulation.ckLocalBranch();

  if ( n != 0) {
    store_face_(n,buffer,
		ifx,ify,ifz,
		0,0,0,
		false,false,false,
		op_array_copy);
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
  int n; 
  char * array;

  FieldFace * field_face = load_face_ (&n, &array,
				       ifx,ify,ifz,
				       icx,icy,icz,
				       false,false,false,
				       op_array_prolong);

  ++loop_refresh_.stop();

  thisProxy[index].x_refresh_fine (n,array,
				   -ifx,-ify,-ifz,
				   icx,icy,icz);
  delete field_face;
	  
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_fine (int n, char * buffer, 
				int ifx, int ify, int ifz,
				int icx, int icy, int icz)
{
  store_face_(n,buffer,
	      ifx,ify,ifz,
	      icx,icy,icz,
	      false,false,false,
	      op_array_prolong);

  std::string refresh_type = simulation()->config()->field_refresh_type;
  
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
