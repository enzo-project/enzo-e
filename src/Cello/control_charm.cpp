// See LICENSE_CELLO file for license and copyright information

/// @file     control_charm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-02-13
/// @brief    Functions controling control flow of charm entry functions
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#ifdef CELLO_VERBOSE
#   define VERBOSE(A)						\
  fprintf (simulation()->fp_debug(),"[%s] %d TRACE %s\n",	\
	   name_.c_str(),__LINE__,A);				\
  if (index_.is_root()) {					\
    Monitor * monitor = simulation()->monitor();		\
    monitor->print("Control", A);				\
  } 
#else
#   define VERBOSE(A) ;
#endif

//----------------------------------------------------------------------

void Block::initial_exit_()
{  control_sync(phase_adapt_enter,"contribute"); }

//----------------------------------------------------------------------

void Block::adapt_enter_()
{
  VERBOSE("adapt_enter");

  performance_switch_ (perf_adapt,__FILE__,__LINE__);

  if ( do_adapt_()) {

    VERBOSE("adapt_enter");

    adapt_begin_();
    
  } else {

    adapt_exit_();

  }
}

//----------------------------------------------------------------------

void Block::adapt_exit_()
{
  VERBOSE("adapt_exit");

  CkCallback callback = CkCallback(CkIndex_Main::p_output_enter(), proxy_main);

  control_sync(callback,"quiescence");
}

//----------------------------------------------------------------------

void Block::output_enter_ ()
{
  VERBOSE("output_enter");

  performance_switch_ (perf_output,__FILE__,__LINE__);

  output_begin_();
}

//----------------------------------------------------------------------

void Block::output_exit_()
{

  VERBOSE("output_exit");

  
  if (index_.is_root()) {

    proxy_simulation[0].p_monitor();
  }

  control_sync(phase_stopping_enter,"contribute");
}

//----------------------------------------------------------------------

void Block::stopping_enter_()
{

  VERBOSE("stopping_enter");

  performance_switch_(perf_stopping,__FILE__,__LINE__);

  stopping_begin_();

}

//----------------------------------------------------------------------

void Block::stopping_exit_()
{

  VERBOSE("stopping_exit");

  if (stop_) {

    control_sync(phase_exit,"contribute");

  } else {

    compute_enter_();

  }

}

//----------------------------------------------------------------------

void Block::compute_enter_ ()
{
  VERBOSE("compute_enter");

  performance_switch_(perf_compute,__FILE__,__LINE__);

  compute_begin_();
}

//----------------------------------------------------------------------

void Block::compute_exit_ ()
{
  VERBOSE("compute_exit");

  adapt_enter_();
}

//----------------------------------------------------------------------

void Block::refresh_enter_() 
{
  VERBOSE("refresh_enter");

  performance_switch_(perf_refresh,__FILE__,__LINE__);

  refresh_begin_();
}

//----------------------------------------------------------------------

void Block::refresh_exit_()
{
  VERBOSE("refresh_exit");

  update_boundary_();

  control_sync(refresh_phase_, refresh_sync_);
}

//======================================================================

void Block::control_sync (CkCallback callback, std::string sync)
{

  if (sync == "contribute") {

    contribute(callback);

  } else if (sync == "quiescence") {

    if (index_.is_root()) CkStartQD(callback);

  } else if (sync == "neighbor") {

    ERROR("Block::control_sync (callback, sync)",
	  "Neighbor synchronization not implemented yet");

  } else {
    ERROR1 ("Block::control_sync()",  
	    "Unknown sync type: sync type %s", 
	    sync.c_str());    
  }
}

//----------------------------------------------------------------------

void Block::control_sync(int phase, std::string sync)
{
  if (sync == "contribute") {

    CkCallback callback;

    if (                     phase == phase_initial_exit) {
      callback = CkCallback (CkIndex_Block::r_initial_exit(NULL), thisProxy);
    } else if (              phase == phase_adapt_enter) {
      callback = CkCallback (CkIndex_Block::r_adapt_enter(NULL), thisProxy);
    } else if (              phase == phase_adapt_next) {
      callback = CkCallback (CkIndex_Block::r_adapt_next(NULL), thisProxy);
    } else if (              phase == phase_adapt_called) {
      callback = CkCallback (CkIndex_Block::r_adapt_called(NULL), thisProxy);
    } else if (              phase == phase_adapt_end) {
      callback = CkCallback (CkIndex_Block::r_adapt_end(NULL), thisProxy);
    } else if (              phase == phase_adapt_exit) {
      callback = CkCallback (CkIndex_Block::r_adapt_exit(NULL), thisProxy);
    } else if (              phase == phase_refresh_enter) {
      callback = CkCallback (CkIndex_Block::r_refresh_enter(NULL), thisProxy);
    } else if (              phase == phase_refresh_exit) {
      callback = CkCallback (CkIndex_Block::r_refresh_exit(NULL), thisProxy);
    } else if (              phase == phase_output_enter) {
      callback = CkCallback (CkIndex_Block::r_output_enter(NULL), thisProxy);
    } else if (              phase == phase_output_exit) {
      callback = CkCallback (CkIndex_Block::r_output_exit(NULL), thisProxy);
    } else if (              phase == phase_compute_enter) {
      callback = CkCallback (CkIndex_Block::r_compute_enter(NULL), thisProxy);
    } else if (              phase == phase_compute_continue) {
      callback = CkCallback (CkIndex_Block::r_compute_continue(NULL), thisProxy);
    } else if (              phase == phase_compute_exit) {
      callback = CkCallback (CkIndex_Block::r_compute_exit(NULL), thisProxy);
    } else if (              phase == phase_stopping_enter) {
      callback = CkCallback (CkIndex_Block::r_stopping_enter(NULL), thisProxy);
    } else if (              phase == phase_stopping_exit) {
      callback = CkCallback (CkIndex_Block::r_stopping_exit(NULL), thisProxy);
    } else if (              phase == phase_exit) {
      callback = CkCallback (CkIndex_Block::r_exit(NULL), thisProxy);
    } else {
      ERROR2 ("Block::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_name[phase],sync.c_str());    
    }

    contribute(callback);

  } else if (sync == "quiescence") {

    // Quiescence detection through the Main chare

    if (index_.is_root()) {
      if (                    phase == phase_initial_exit) {
	CkStartQD(CkCallback(CkIndex_Main::p_initial_exit(), proxy_main));
      } else if (             phase == phase_adapt_enter) {
	CkStartQD(CkCallback(CkIndex_Main::p_adapt_enter(), proxy_main));
      } else if (             phase == phase_adapt_next) {
	CkStartQD(CkCallback(CkIndex_Main::p_adapt_next(), proxy_main));
      } else if (             phase == phase_adapt_called) {
	CkStartQD(CkCallback(CkIndex_Main::p_adapt_called(), proxy_main));
      } else if (             phase == phase_adapt_end) {
	CkStartQD(CkCallback(CkIndex_Main::p_adapt_end(), proxy_main));
      } else if (             phase == phase_adapt_exit) {
	CkStartQD(CkCallback(CkIndex_Main::p_adapt_exit(), proxy_main));
      } else if (             phase == phase_refresh_enter) {
      	CkStartQD(CkCallback(CkIndex_Main::p_refresh_enter(), proxy_main));
      } else if (             phase == phase_refresh_exit) {
      	CkStartQD(CkCallback(CkIndex_Main::p_refresh_exit(), proxy_main));
      } else if (             phase == phase_output_enter) {
	CkStartQD(CkCallback(CkIndex_Main::p_output_enter(), proxy_main));
      } else if (             phase == phase_output_exit) {
	CkStartQD(CkCallback(CkIndex_Main::p_output_exit(), proxy_main));
      } else if (             phase == phase_compute_enter) {
	CkStartQD(CkCallback(CkIndex_Main::p_compute_enter(), proxy_main));
      } else if (             phase == phase_compute_continue) {
	CkStartQD(CkCallback(CkIndex_Main::p_compute_continue(), proxy_main));
      } else if (             phase == phase_compute_exit) {
	CkStartQD(CkCallback(CkIndex_Main::p_compute_exit(), proxy_main));
      } else if (             phase == phase_stopping_enter) {
	CkStartQD(CkCallback(CkIndex_Main::p_stopping_enter(), proxy_main));
      } else if (             phase == phase_stopping_exit) {
	CkStartQD(CkCallback(CkIndex_Main::p_stopping_exit(), proxy_main));
      } else if (             phase == phase_exit) {
	CkStartQD(CkCallback(CkIndex_Main::p_exit(), proxy_main));
      } else {
	ERROR2 ("Block::control_sync()",  
		"Unknown phase: phase %s sync type %s", 
		phase_name[phase],sync.c_str());    
      }

    }

  } else if (sync == "neighbor") {

    control_sync_neighbor_(phase);

  } else {
    ERROR2 ("Block::control_sync()",  
	    "Unknown sync type: phase %s sync type %s", 
	    phase_name[phase],sync.c_str());    
  }
}

//----------------------------------------------------------------------

void Block::control_sync_count_ (int phase, int count)
{
  if (count != 0)  max_sync_[phase] = count;

  ++count_sync_[phase];

  // max_sync reached: continue and reset counter
  if (max_sync_[phase] > 0 && count_sync_[phase] >= max_sync_[phase]) {
    max_sync_[phase] = 0;
    count_sync_[phase] = 0;
    control_call_phase_(phase);
  }
}

//----------------------------------------------------------------------

void Block::control_sync_neighbor_(int phase)
{
  if (!is_leaf()) {

#ifdef CELLO_DEBUG
    fprintf (simulation()->fp_debug(),"%d %s called %s\n",__LINE__,name_.c_str(),phase_name[phase]);
#endif
    control_call_phase_ (phase);

  } else {

    const int level        = this->level();
    const int rank = this->rank();

    const int min_face_rank = 0;
    ItFace it_face = this->it_face(min_face_rank,index_);

    int of3[3];

    int num_neighbors = 0;

    while (it_face.next(of3)) {

      int ic3[3] = {0,0,0};

      Index index_neighbor = neighbor_(of3);

      const int level_face = face_level (of3);

      if (level_face == level) {

	// SEND-SAME: Face and level are sent to unique
	// neighboring block in the same level

#ifdef CELLO_DEBUG
	fprintf (simulation()->fp_debug(),
		 "%d %s calling p_control_sync phase %s leaf %d block %s\n",
		 __LINE__,name_.c_str(), phase_name[phase], is_leaf(),
		 index_neighbor.bit_string(-1,2).c_str());
#endif
	thisProxy[index_neighbor].p_control_sync_count(phase);

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

#ifdef CELLO_DEBUG
	  fprintf (simulation()->fp_debug(),
		   "%d %s calling p_control_sync phase %s leaf %d block %s\n",
		   __LINE__,name_.c_str(), phase_name[phase], is_leaf(),
		   index_uncle.bit_string(-1,2).c_str());
#endif
	  thisProxy[index_uncle].p_control_sync_count(phase);

	  ++num_neighbors;

	}

      } else if (level_face == level + 1) {

	// SEND-FINE: Face and level are sent to all nibling
	// blocks in the next-finer level along the face.

	const int if3[3] = {-of3[0],-of3[1],-of3[2]};
	ItChild it_child(rank,if3);
	while (it_child.next(ic3)) {
	  Index index_nibling = index_neighbor.index_child(ic3);

#ifdef CELLO_DEBUG
	  fprintf (simulation()->fp_debug(),
		   "%d %s calling p_control_sync phase %s leaf %d block %s\n",
		   __LINE__,name_.c_str(), phase_name[phase], is_leaf(),
		   index_nibling.bit_string(-1,2).c_str());
#endif
	  thisProxy[index_nibling].p_control_sync_count(phase);

	  ++num_neighbors;
	}

      } else {
	std::string bit_str = index_.bit_string(-1,2);
	WARNING7 ("Block::control_sync_neighbor_()",
		  "phase %d name %s level %d and face (%d %d %d) "
		  "level %d differ by more than 1",
		  phase, name().c_str(), level, of3[0],of3[1],of3[2], level_face);
      }

    }
    control_sync_count_(phase,num_neighbors + 1);
  }
}

//----------------------------------------------------------------------

void Block::control_call_phase_ (int phase)
{
#ifdef CELLO_DEBUG
  fprintf (simulation()->fp_debug(),"%d %s called %s\n",
	   __LINE__,name_.c_str(),phase_name[phase]);
#endif
  CkCallback callback;
  if      (phase == phase_initial_exit)     initial_exit_();
  else if (phase == phase_adapt_called)     adapt_called_();
  else if (phase == phase_adapt_enter)      adapt_enter_();
  else if (phase == phase_adapt_next)       adapt_next_();
  else if (phase == phase_adapt_end)        adapt_end_();
  else if (phase == phase_adapt_exit)       adapt_exit_();
  else if (phase == phase_refresh_enter)    refresh_enter_();
  else if (phase == phase_refresh_exit)     refresh_exit_();
  else if (phase == phase_output_enter)     output_enter_();
  else if (phase == phase_output_exit)      output_exit_();
  else if (phase == phase_compute_enter)    compute_enter_();
  else if (phase == phase_compute_continue) compute_continue_();
  else if (phase == phase_compute_exit)     compute_exit_();
  else if (phase == phase_stopping_enter)   stopping_enter_();
  else if (phase == phase_stopping_exit)    stopping_exit_();
  else if (phase == phase_exit)             exit_();
  else 
    ERROR1 ("Block::control_call_phase_()","Unknown phase: phase %s", 
	    phase_name[phase]);    
}
//======================================================================

