// See LICENSE_CELLO file for license and copyright information

/// @file     charm_control.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-02-13
/// @brief    Functions controling control flow of charm entry functions

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

const char * phase_string [] = {
  "unknown",
  "adapt_called",
  "adapt_enter",
  "adapt_next",
  "adapt_exit",
  "refresh_enter",
  "refresh_exit"
};

//----------------------------------------------------------------------

void CommBlock::control_sync(int phase)
{
  
  std::string sync_type;

  if (phase == phase_sync_adapt_enter) {
    sync_type = simulation()->config()->mesh_sync_adapt_enter;
  } else if (phase == phase_sync_adapt_next) {
    sync_type = simulation()->config()->mesh_sync_adapt_next;
  } else if (phase == phase_sync_adapt_called) {
    sync_type = simulation()->config()->mesh_sync_adapt_called;
  } else if (phase == phase_sync_adapt_exit) {
    sync_type = simulation()->config()->mesh_sync_adapt_exit;
  } else if (phase == phase_sync_refresh_enter) {
    sync_type = simulation()->config()->mesh_sync_refresh_enter;
  } else if (phase == phase_sync_refresh_exit) {
    sync_type = simulation()->config()->mesh_sync_refresh_exit;
  } else {
    ERROR1 ("CommBlock::control_sync()",
	    "Unknown phase: phase %s",
	    phase_string[phase]);    
  }

  if (index().is_root()) {
    PARALLEL_PRINTF ("DEBUG: %s %s\n",phase_string[phase],sync_type.c_str());
  }

  if (sync_type == "contribute") {

    CkCallback cb;

    if (phase == phase_sync_adapt_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_enter(NULL), thisProxy);
    } else if (phase == phase_sync_adapt_next) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_next(NULL), thisProxy);
    } else if (phase == phase_sync_adapt_called) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_called(NULL), thisProxy);
    } else if (phase == phase_sync_adapt_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_exit(NULL), thisProxy);
    } else if (phase == phase_sync_refresh_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_refresh_enter(NULL), thisProxy);
    } else if (phase == phase_sync_refresh_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_refresh_exit(NULL), thisProxy);
    } else {
      ERROR2 ("CommBlock::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_string[phase],sync_type.c_str());    
    }

    contribute(cb);

  } else if (sync_type == "quiescence") {

    CkCallback cb;

    if (phase == phase_sync_adapt_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_enter(), thisProxy[thisIndex]);
    } else if (phase == phase_sync_adapt_next) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_next(), thisProxy[thisIndex]);
    } else if (phase == phase_sync_adapt_called) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_called(), thisProxy[thisIndex]);
    } else if (phase == phase_sync_adapt_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_exit(), thisProxy[thisIndex]);
    } else if (phase == phase_sync_refresh_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_refresh_enter(), thisProxy[thisIndex]);
    } else if (phase == phase_sync_refresh_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_refresh_exit(), thisProxy[thisIndex]);
    } else {
      ERROR2 ("CommBlock::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_string[phase],sync_type.c_str());    
    }

    CkStartQD (cb);

  } else if (sync_type == "neighbor") {

    control_sync_neighbor_(phase);

  } else if (sync_type == "array") {

    if (phase == phase_sync_adapt_enter) {
      if (index().is_root()) thisProxy.p_adapt_enter();
    } else if (phase == phase_sync_adapt_next) {
      if (index().is_root()) thisProxy.p_adapt_next();
    } else if (phase == phase_sync_adapt_called) {
      if (index().is_root()) thisProxy.p_adapt_called();
    } else if (phase == phase_sync_adapt_exit) {
      if (index().is_root()) thisProxy.p_adapt_exit();
    } else if (phase == phase_sync_refresh_enter) {
      if (index().is_root()) thisProxy.p_refresh_enter();
    } else if (phase == phase_sync_refresh_exit) {
      if (index().is_root()) thisProxy.p_refresh_exit();
    } else {
      ERROR2 ("CommBlock::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_string[phase],sync_type.c_str());    
    }

  } else if (sync_type == "none") {

    if (phase == phase_sync_adapt_enter) {
      adapt_enter_();
    } else if (phase == phase_sync_adapt_next) {
      adapt_next_();
    } else if (phase == phase_sync_adapt_called) {
      adapt_called_();
    } else if (phase == phase_sync_adapt_exit) {
      adapt_exit_();
    } else if (phase == phase_sync_refresh_enter) {
      refresh_enter_();
    } else if (phase == phase_sync_refresh_exit) {
      refresh_exit_();
    } else {
      ERROR2 ("CommBlock::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_string[phase],sync_type.c_str());    
    }

  } else {
    ERROR2 ("CommBlock::control_sync()",  
	    "Unknown sync type: phase %s sync type %s", 
	    phase_string[phase],sync_type.c_str());    
  }
}

//----------------------------------------------------------------------

void CommBlock::control_sync_count_ (int phase, int count)
{
  if (count != 0)  max_sync_[phase] = count;

  ++count_sync_[phase];

  if (0 < max_sync_[phase] && max_sync_[phase] <= count_sync_[phase]) {
    max_sync_[phase] = 0;
    count_sync_[phase] = 0;
    control_call_phase_(phase);
  }
}

//----------------------------------------------------------------------

void CommBlock::control_sync_neighbor_(int phase)
{
  if (!is_leaf()) {

    control_call_phase_ (phase);

  } else {

    const int level        = this->level();
    const int rank = simulation()->dimension();
    const int rank_refresh = simulation()->config()->field_refresh_rank;

    ItFace it_face(rank,rank_refresh);
    int of3[3];

    int num_neighbors = 0;

    while (it_face.next(of3)) {

      int ic3[3] = {0,0,0};

      Index index_neighbor = neighbor_(of3);

      const int level_face = face_level (of3);

      if (level_face == level) {

	// SEND-SAME: Face and level are sent to unique
	// neighboring block in the same level

	//--------------------------------------------------
	// ENTRY: #1 CommBlock::control_sync_neighbor_() -> p_control_sync_count()
	// ENTRY: same-level neighbor
	// ENTRY: adapt phase
	//--------------------------------------------------
	thisProxy[index_neighbor].p_control_sync_count(phase,0);
	//--------------------------------------------------

	++num_neighbors;

      } else if (level_face == level - 1) {

	// SEND-COARSE: Face, level, and child indices are sent to
	// unique neighboring block in the next-coarser level


	index_.child (level,&ic3[0],&ic3[1],&ic3[2]);

	int op3[3];
	parent_face_(op3,of3,ic3);

	// avoid redundant calls
	if (op3[0]==of3[0] && 
	    op3[1]==of3[1] && 
	    op3[2]==of3[2]) {

	  Index index_uncle = index_neighbor.index_parent();

	  //--------------------------------------------------
	  // ENTRY: #2 CommBlock::control_sync_neighbor_() -> p_control_sync_count()
	  // ENTRY: coarse-level neighbor
	  // ENTRY: adapt phase
	  //--------------------------------------------------
	  thisProxy[index_uncle].p_control_sync_count(phase,0);
	  //--------------------------------------------------

	  ++num_neighbors;

	}

      } else if (level_face == level + 1) {

	// SEND-FINE: Face and level are sent to all nibling
	// blocks in the next-finer level along the face.

	const int if3[3] = {-of3[0],-of3[1],-of3[2]};
	ItChild it_child(rank,if3);
	while (it_child.next(ic3)) {
	  Index index_nibling = index_neighbor.index_child(ic3);

	  // --------------------------------------------------
	  // ENTRY: #3 CommBlock::control_sync_neighbor_() -> p_control_sync_count()
	  // ENTRY: fine-level neighbor
	  // ENTRY: adapt phase
	  // --------------------------------------------------
	  thisProxy[index_nibling].p_control_sync_count(phase,0);
	  // --------------------------------------------------

	  ++num_neighbors;
	}

      } else {
	std::string bit_str = index_.bit_string(-1,2);
	WARNING4 ("CommBlock::control_sync_neighbor_()",
		  "phase %d %s level %d and face level %d differ by more than 1",
		  phase,bit_str.c_str(),level,level_face);
      }

    }
    control_sync_count_(phase,num_neighbors + 1);

  }
}

//----------------------------------------------------------------------

void CommBlock::control_call_phase_ (int phase)
{
  if (phase == phase_sync_adapt_called)  adapt_called_() ;
  if (phase == phase_sync_adapt_enter)   adapt_enter_() ;
  if (phase == phase_sync_adapt_next)    adapt_next_() ;
  if (phase == phase_sync_adapt_exit)    adapt_exit_() ;
  if (phase == phase_sync_refresh_enter) refresh_enter_() ;
  if (phase == phase_sync_refresh_exit)  refresh_exit_();
}

//======================================================================

