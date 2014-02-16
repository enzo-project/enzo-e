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

//----------------------------------------------------------------------

void CommBlock::control_sync_ (int phase, int count)
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
	// ENTRY: #1 CommBlock::control_sync_neighbor_() -> p_control_sync()
	// ENTRY: same-level neighbor
	// ENTRY: adapt phase
	//--------------------------------------------------
	thisProxy[index_neighbor].p_control_sync(phase,0);
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
	  // ENTRY: #2 CommBlock::control_sync_neighbor_() -> p_control_sync()
	  // ENTRY: coarse-level neighbor
	  // ENTRY: adapt phase
	  //--------------------------------------------------
	  thisProxy[index_uncle].p_control_sync(phase,0);
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
	  // ENTRY: #3 CommBlock::control_sync_neighbor_() -> p_control_sync()
	  // ENTRY: fine-level neighbor
	  // ENTRY: adapt phase
	  // --------------------------------------------------
	  thisProxy[index_nibling].p_control_sync(phase,0);
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
    control_sync_(phase,num_neighbors + 1);

  }
}

//----------------------------------------------------------------------

void CommBlock::control_call_phase_ (int phase)
{
    if (phase == phase_sync_adapt_called) adapt_called_() ;
    if (phase == phase_sync_adapt_next)   adapt_next_() ;
    if (phase == phase_sync_adapt_exit)   adapt_exit_() ;
    if (phase == phase_sync_refresh)      refresh_exit_();
}

//======================================================================

