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

// #define TRACE_CONTROL

#ifdef CELLO_VERBOSE
#   define VERBOSE(A) if (index_.is_root()){PARALLEL_PRINTF (A "\n"); fflush(stdout);}
#else
#   define VERBOSE(A) ;
#endif

const char * phase_string [] = {
  "unknown",
  "adapt_called",
  "adapt_enter",
  "adapt_next",
  "adapt_end",
  "refresh_enter",
  "refresh_exit"
  "output_enter",
  "output_exit",
  "compute_enter",
  "compute_exit",
  "stopping_enter",
  "stopping_exit",
  "exit"
};

//----------------------------------------------------------------------

void CommBlock::adapt_enter_()
{
  VERBOSE("ENTER PHASE ADAPT") ;

  performance_switch_ (perf_adapt,__FILE__,__LINE__);

  if ( do_adapt_()) {

    adapt_begin_();

  } else {

    adapt_exit_();

  }
}

//----------------------------------------------------------------------

void CommBlock::adapt_exit_()
{

  VERBOSE("EXIT   PHASE ADAPT") ;

  next_phase_ = phase_stopping;

  const int initial_cycle = simulation()->config()->initial_cycle;
  const bool is_first_cycle = (initial_cycle == cycle());
  const int level_maximum = simulation()->config()->initial_max_level;

  bool adapt_again = (is_first_cycle && adapt_step_++ < level_maximum);

  if (adapt_again) {

    control_sync (phase_sync_adapt_enter);

  } else {

    control_sync (phase_sync_refresh_enter);
  }

}

//----------------------------------------------------------------------

void CommBlock::compute_enter_ ()
{
  VERBOSE("ENTER PHASE COMPUTE");

  performance_switch_(perf_compute,__FILE__,__LINE__);

  compute_begin_();
}

//----------------------------------------------------------------------

void CommBlock::compute_exit_ ()
{

  VERBOSE("EXIT   PHASE COMPUTE") ;

  next_phase_ = phase_adapt;

  control_sync(phase_sync_refresh_enter);
}

//----------------------------------------------------------------------

void CommBlock::refresh_enter_() 
{
  VERBOSE ("ENTER PHASE REFRESH");

  performance_switch_(perf_refresh,__FILE__,__LINE__);

  refresh_begin_();
}

//----------------------------------------------------------------------

void CommBlock::refresh_exit_()
{

  VERBOSE("EXIT   PHASE REFRESH") ;

  if (next_phase_ == phase_stopping) {

    control_sync(phase_sync_stopping_enter);

  }  else if (next_phase_ == phase_adapt) {

    control_sync (phase_sync_adapt_enter);

  } else {

    ERROR1 ("CommBlock::q_refresh_exit()", 
	    "Unknown next_phase %d", next_phase_);
  }
}

//----------------------------------------------------------------------

void CommBlock::stopping_enter_()
{

  VERBOSE ("ENTER PHASE STOPPING");

  performance_switch_(perf_stopping,__FILE__,__LINE__);

  stopping_begin_();

}

//----------------------------------------------------------------------

void CommBlock::stopping_exit_()
{

  VERBOSE ("EXIT  PHASE STOPPING");
  if (stop_) {

    control_sync(phase_sync_exit);

  } else {

    control_sync(phase_sync_output_enter);

  }

}

//----------------------------------------------------------------------

void CommBlock::output_enter_ ()
{

  VERBOSE ("ENTER PHASE OUTPUT");

  performance_switch_ (perf_output,__FILE__,__LINE__);

  output_begin_();

}

//----------------------------------------------------------------------

void CommBlock::output_exit_()
{

  VERBOSE("EXIT   PHASE OUTPUT") ;

  if (index_.is_root()) {

    proxy_simulation.p_performance_output();

    Memory::instance()->reset_high();

    Monitor * monitor = simulation()->monitor();
    monitor-> print("", "-------------------------------------");
    monitor-> print("Simulation", "cycle %04d", cycle_);
    monitor-> print("Simulation", "time-sim %15.12f",time_);
    monitor-> print("Simulation", "dt %15.12g", dt_);
  }

  control_sync(phase_sync_compute_enter);
}

//======================================================================

void CommBlock::control_sync(int phase)
{
  
  std::string sync_type;

  if (phase == phase_sync_adapt_enter) {
    sync_type = simulation()->config()->mesh_sync_adapt_enter;
  } else if (phase == phase_sync_adapt_next) {
    sync_type = simulation()->config()->mesh_sync_adapt_next;
  } else if (phase == phase_sync_adapt_called) {
    sync_type = simulation()->config()->mesh_sync_adapt_called;
  } else if (phase == phase_sync_adapt_end) {
    sync_type = simulation()->config()->mesh_sync_adapt_end;
  } else if (phase == phase_sync_refresh_enter) {
    sync_type = simulation()->config()->mesh_sync_refresh_enter;
  } else if (phase == phase_sync_refresh_exit) {
    sync_type = simulation()->config()->mesh_sync_refresh_exit;
  } else if (phase == phase_sync_output_enter) {
    sync_type = simulation()->config()->mesh_sync_output_enter;
  } else if (phase == phase_sync_output_exit) {
    sync_type = simulation()->config()->mesh_sync_output_exit;
  } else if (phase == phase_sync_compute_enter) {
    sync_type = simulation()->config()->mesh_sync_compute_enter;
  } else if (phase == phase_sync_compute_exit) {
    sync_type = simulation()->config()->mesh_sync_compute_exit;
  } else if (phase == phase_sync_stopping_enter) {
    sync_type = simulation()->config()->mesh_sync_stopping_enter;
  } else if (phase == phase_sync_stopping_exit) {
    sync_type = simulation()->config()->mesh_sync_stopping_exit;
  } else if (phase == phase_sync_exit) {
    sync_type = simulation()->config()->mesh_sync_exit;
  } else {
    ERROR2 ("CommBlock::control_sync()",
	    "Unknown phase: phase %s (%d)",
	    phase_string[phase],phase);    
  }

  if (sync_type == "contribute") {

    CkCallback cb;

    if (phase == phase_sync_adapt_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_enter(NULL), thisProxy);
    } else if (phase == phase_sync_adapt_next) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_next(NULL), thisProxy);
    } else if (phase == phase_sync_adapt_called) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_called(NULL), thisProxy);
    } else if (phase == phase_sync_adapt_end) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_end(NULL), thisProxy);
    } else if (phase == phase_sync_refresh_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_refresh_enter(NULL), thisProxy);
    } else if (phase == phase_sync_refresh_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_refresh_exit(NULL), thisProxy);
    } else if (phase == phase_sync_output_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_output_enter(NULL), thisProxy);
    } else if (phase == phase_sync_output_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_output_exit(NULL), thisProxy);
    } else if (phase == phase_sync_compute_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_compute_enter(NULL), thisProxy);
    } else if (phase == phase_sync_compute_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_compute_exit(NULL), thisProxy);
    } else if (phase == phase_sync_stopping_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_stopping_enter(NULL), thisProxy);
    } else if (phase == phase_sync_stopping_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_stopping_exit(NULL), thisProxy);
    } else if (phase == phase_sync_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_exit(NULL), thisProxy);
    } else {
      ERROR2 ("CommBlock::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_string[phase],sync_type.c_str());    
    }

    contribute(cb);

  } else if (sync_type == "quiescence") {

    CkCallback cb;

    if (phase == phase_sync_adapt_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_adapt_next) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_next(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_adapt_called) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_called(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_adapt_end) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_end(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_refresh_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_refresh_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_refresh_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_refresh_exit(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_output_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_output_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_output_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_output_exit(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_compute_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_compute_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_compute_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_compute_exit(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_stopping_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_stopping_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_stopping_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_stopping_exit(),
		       thisProxy[thisIndex]);
    } else if (phase == phase_sync_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_exit(),
		       thisProxy[thisIndex]);
    } else {
      ERROR2 ("CommBlock::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_string[phase],sync_type.c_str());    
    }

    CkStartQD (cb);

   } else if (sync_type == "main-qd") {

     // Quiescence detection through the Main chare

    if (index_.is_root()) {
      if (phase == phase_sync_adapt_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_adapt_enter(), proxy_main));
      } else if (phase == phase_sync_adapt_next) {
	CkStartQD(CkCallback(CkIndex_Main::q_adapt_next(), proxy_main));
      } else if (phase == phase_sync_adapt_called) {
	CkStartQD(CkCallback(CkIndex_Main::q_adapt_called(), proxy_main));
      } else if (phase == phase_sync_adapt_end) {
	CkStartQD(CkCallback(CkIndex_Main::q_adapt_end(), proxy_main));
      } else if (phase == phase_sync_refresh_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_refresh_enter(), proxy_main));
      } else if (phase == phase_sync_refresh_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_refresh_exit(), proxy_main));
      } else if (phase == phase_sync_output_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_output_enter(), proxy_main));
      } else if (phase == phase_sync_output_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_output_exit(), proxy_main));
      } else if (phase == phase_sync_compute_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_compute_enter(), proxy_main));
      } else if (phase == phase_sync_compute_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_compute_exit(), proxy_main));
      } else if (phase == phase_sync_stopping_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_stopping_enter(), proxy_main));
      } else if (phase == phase_sync_stopping_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_stopping_exit(), proxy_main));
      } else if (phase == phase_sync_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_exit(), proxy_main));
      } else {
	ERROR2 ("CommBlock::control_sync()",  
		"Unknown phase: phase %s sync type %s", 
		phase_string[phase],sync_type.c_str());    
      }

    }

  } else if (sync_type == "neighbor") {

    control_sync_neighbor_(phase);

  } else if (sync_type == "array") {

    if (phase == phase_sync_adapt_enter) {
      if (index().is_root()) thisProxy.p_adapt_enter();
    } else if (phase == phase_sync_adapt_next) {
      if (index().is_root()) thisProxy.p_adapt_next();
    } else if (phase == phase_sync_adapt_called) {
      if (index().is_root()) thisProxy.p_adapt_called();
    } else if (phase == phase_sync_adapt_end) {
      if (index().is_root()) thisProxy.p_adapt_end();
    } else if (phase == phase_sync_refresh_enter) {
      if (index().is_root()) thisProxy.p_refresh_enter();
    } else if (phase == phase_sync_refresh_exit) {
      if (index().is_root()) thisProxy.p_refresh_exit();
    } else if (phase == phase_sync_output_enter) {
      if (index().is_root()) thisProxy.p_output_enter();
    } else if (phase == phase_sync_output_exit) {
      if (index().is_root()) thisProxy.p_output_exit();
    } else if (phase == phase_sync_compute_enter) {
      if (index().is_root()) thisProxy.p_compute_enter();
    } else if (phase == phase_sync_compute_exit) {
      if (index().is_root()) thisProxy.p_compute_exit();
    } else if (phase == phase_sync_stopping_enter) {
      if (index().is_root()) thisProxy.p_stopping_enter();
    } else if (phase == phase_sync_stopping_exit) {
      if (index().is_root()) thisProxy.p_stopping_exit();
    } else if (phase == phase_sync_exit) {
      if (index().is_root()) thisProxy.p_exit();
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
    } else if (phase == phase_sync_adapt_end) {
      adapt_end_();
    } else if (phase == phase_sync_refresh_enter) {
      refresh_enter_();
    } else if (phase == phase_sync_refresh_exit) {
      refresh_exit_();
    } else if (phase == phase_sync_output_enter) {
      output_enter_();
    } else if (phase == phase_sync_output_exit) {
      output_exit_();
    } else if (phase == phase_sync_compute_enter) {
      compute_enter_();
    } else if (phase == phase_sync_compute_exit) {
      compute_exit_();
    } else if (phase == phase_sync_stopping_enter) {
      stopping_enter_();
    } else if (phase == phase_sync_stopping_exit) {
      stopping_exit_();
    } else if (phase == phase_sync_exit) {
      exit_();
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
#ifdef TRACE_CONTROL
    char buffer[255];
    sprintf (buffer,"calling phase %s",phase_string[phase]);
    index_.print(buffer,-1,3,false,simulation());
#endif
    control_call_phase_(phase);
  }
}

//----------------------------------------------------------------------

void CommBlock::control_sync_neighbor_(int phase)
{
  if (!is_leaf()) {

#ifdef TRACE_CONTROL
    char buffer[255];
    sprintf (buffer,"calling phase %s",phase_string[phase]);
    index_.print(buffer,-1,3,false,simulation());
#endif
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
#ifdef TRACE_CONTROL
  char buffer[255];
  sprintf (buffer,"called phase %s",phase_string[phase]);
  index_.print(buffer,-1,3,false,simulation());
#endif
  if (phase == phase_sync_adapt_called) {
    adapt_called_() ;
  } else if (phase == phase_sync_adapt_enter) {
    adapt_enter_() ;
  } else if (phase == phase_sync_adapt_next) {
    adapt_next_() ;
  } else if (phase == phase_sync_adapt_end) {
    adapt_end_() ;
  } else if (phase == phase_sync_refresh_enter) {
    refresh_enter_() ;
  } else if (phase == phase_sync_refresh_exit) {
    refresh_exit_();
  } else if (phase == phase_sync_output_enter) {
    output_enter_() ;
  } else if (phase == phase_sync_output_exit) {
    output_exit_();
  } else if (phase == phase_sync_compute_enter) {
    compute_enter_() ;
  } else if (phase == phase_sync_compute_exit) {
    compute_exit_();
  } else if (phase == phase_sync_stopping_enter) {
    stopping_enter_() ;
  } else if (phase == phase_sync_stopping_exit) {
    stopping_exit_();
  } else if (phase == phase_sync_exit) {
    exit_();
  } else {
    ERROR1 ("CommBlock::control_call_phase_()",  
	    "Unknown phase: phase %s", 
	    phase_string[phase]);    
  }
}
//======================================================================

