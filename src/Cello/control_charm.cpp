// See LICENSE_CELLO file for license and copyright information

/// @file     control_charm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-02-13
/// @brief    Functions controling control flow of charm entry functions
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#ifdef CELLO_VERBOSE
#   define VERBOSE(A)					\
  fprintf (simulation()->fp_debug(),"[%s] %d TRACE %s\n",name_.c_str(),__LINE__,A); \
  if (index_.is_root()) {				\
    Monitor * monitor = simulation()->monitor();	\
    monitor->print("Control", A);			\
  } 
#else
#   define VERBOSE(A) ;
#endif

const char * phase_string [] = {
  "unknown",
  "adapt_enter",
  "adapt_called",
  "adapt_next",
  "adapt_end",
  "adapt_exit",
  "compute_enter",
  "compute_exit",
  "output_enter",
  "output_exit",
  "refresh_enter",
  "refresh_exit",
  "stopping_enter",
  "stopping_exit",
  "exit"
};

//----------------------------------------------------------------------

void CommBlock::adapt_enter_()
{

  VERBOSE("adapt_enter");

  performance_switch_ (perf_adapt,__FILE__,__LINE__);

  if ( do_adapt_()) {

    VERBOSE("adapt_enter");

    adapt_begin_();
    
  } else {

    control_sync(sync_adapt_exit,"none");

  }
}

//----------------------------------------------------------------------

void CommBlock::adapt_exit_()
{

  VERBOSE("adapt_exit");

  next_phase_ = phase_stopping;

  const int initial_cycle = simulation()->config()->initial_cycle;
  const bool is_first_cycle = (initial_cycle == cycle());
  const int level_maximum = simulation()->config()->mesh_max_level;

  bool adapt_again = (is_first_cycle && adapt_step_++ < level_maximum);


  if (adapt_again) {

    control_sync (sync_adapt_enter,"array");

  } else {

    control_sync (sync_refresh_enter,"array");
  }

}

//----------------------------------------------------------------------

void CommBlock::refresh_enter_() 
{
  VERBOSE("refresh_enter");

  performance_switch_(perf_refresh,__FILE__,__LINE__);

  refresh_begin_();
}

//----------------------------------------------------------------------

void CommBlock::refresh_exit_()
{

  VERBOSE("refresh_exit");

  if (next_phase_ == phase_stopping) {

    control_sync(sync_output_enter,"array");

  }  else if (next_phase_ == phase_adapt) {

    control_sync (sync_adapt_enter,"array");

  } else {

    ERROR1 ("CommBlock::q_refresh_exit()", 
	    "Unknown next_phase %d", next_phase_);
  }
}

//----------------------------------------------------------------------

void CommBlock::compute_enter_ ()
{
  VERBOSE("compute_enter");

  performance_switch_(perf_compute,__FILE__,__LINE__);

  compute_begin_();
}

//----------------------------------------------------------------------

void CommBlock::compute_exit_ ()
{

  VERBOSE("compute_exit");

  next_phase_ = phase_adapt;

  control_sync(sync_refresh_enter,"array");
}

//----------------------------------------------------------------------

void CommBlock::output_enter_ ()
{

  VERBOSE("output_enter");

  performance_switch_ (perf_output,__FILE__,__LINE__);

  output_begin_();

}

//----------------------------------------------------------------------

void CommBlock::output_exit_()
{

  VERBOSE("output_exit");

  
  if (index_.is_root()) {

    proxy_simulation[0].p_monitor();
  }

  control_sync(sync_stopping_enter,"none");
}

//----------------------------------------------------------------------

void CommBlock::stopping_enter_()
{

  VERBOSE("stopping_enter");

  performance_switch_(perf_stopping,__FILE__,__LINE__);

  stopping_begin_();

}

//----------------------------------------------------------------------

void CommBlock::stopping_exit_()
{

  VERBOSE("stopping_exit");

  if (stop_) {

    control_sync(sync_exit,"contribute");

  } else {

    control_sync(sync_compute_enter,"none");

  }

}

//======================================================================

void CommBlock::control_sync(int phase, std::string sync)
{

  if (sync == "contribute") {

    CkCallback cb;

    if (       phase == sync_adapt_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_enter(NULL), thisProxy);
    } else if (phase == sync_adapt_next) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_next(NULL), thisProxy);
    } else if (phase == sync_adapt_called) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_called(NULL), thisProxy);
    } else if (phase == sync_adapt_end) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_end(NULL), thisProxy);
    } else if (phase == sync_adapt_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_exit(NULL), thisProxy);
    } else if (phase == sync_adapt_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_adapt_exit(NULL), thisProxy);
    } else if (phase == sync_refresh_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_refresh_enter(NULL), thisProxy);
    } else if (phase == sync_refresh_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_refresh_exit(NULL), thisProxy);
    } else if (phase == sync_output_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_output_enter(NULL), thisProxy);
    } else if (phase == sync_output_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_output_exit(NULL), thisProxy);
    } else if (phase == sync_compute_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_compute_enter(NULL), thisProxy);
    } else if (phase == sync_compute_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_compute_exit(NULL), thisProxy);
    } else if (phase == sync_stopping_enter) {
      cb = CkCallback (CkIndex_CommBlock::r_stopping_enter(NULL), thisProxy);
    } else if (phase == sync_stopping_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_stopping_exit(NULL), thisProxy);
    } else if (phase == sync_exit) {
      cb = CkCallback (CkIndex_CommBlock::r_exit(NULL), thisProxy);
    } else {
      ERROR2 ("CommBlock::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_string[phase],sync.c_str());    
    }

    contribute(cb);

  } else if (sync == "quiescence") {

    CkCallback cb;

    if (phase == sync_adapt_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_adapt_next) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_next(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_adapt_called) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_called(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_adapt_end) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_end(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_adapt_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_adapt_exit(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_refresh_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_refresh_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_refresh_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_refresh_exit(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_output_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_output_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_output_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_output_exit(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_compute_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_compute_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_compute_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_compute_exit(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_stopping_enter) {
      cb = CkCallback (CkIndex_CommBlock::q_stopping_enter(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_stopping_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_stopping_exit(),
		       thisProxy[thisIndex]);
    } else if (phase == sync_exit) {
      cb = CkCallback (CkIndex_CommBlock::q_exit(),
		       thisProxy[thisIndex]);
    } else {
      ERROR2 ("CommBlock::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_string[phase],sync.c_str());    
    }

    CkStartQD (cb);

   } else if (sync == "main-qd") {

     // Quiescence detection through the Main chare

    if (index_.is_root()) {
      if (phase == sync_adapt_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_adapt_enter(), proxy_main));
      } else if (phase == sync_adapt_next) {
	CkStartQD(CkCallback(CkIndex_Main::q_adapt_next(), proxy_main));
      } else if (phase == sync_adapt_called) {
	CkStartQD(CkCallback(CkIndex_Main::q_adapt_called(), proxy_main));
      } else if (phase == sync_adapt_end) {
	CkStartQD(CkCallback(CkIndex_Main::q_adapt_end(), proxy_main));
      } else if (phase == sync_adapt_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_adapt_exit(), proxy_main));
      } else if (phase == sync_refresh_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_refresh_enter(), proxy_main));
      } else if (phase == sync_refresh_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_refresh_exit(), proxy_main));
      } else if (phase == sync_output_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_output_enter(), proxy_main));
      } else if (phase == sync_output_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_output_exit(), proxy_main));
      } else if (phase == sync_compute_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_compute_enter(), proxy_main));
      } else if (phase == sync_compute_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_compute_exit(), proxy_main));
      } else if (phase == sync_stopping_enter) {
	CkStartQD(CkCallback(CkIndex_Main::q_stopping_enter(), proxy_main));
      } else if (phase == sync_stopping_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_stopping_exit(), proxy_main));
      } else if (phase == sync_exit) {
	CkStartQD(CkCallback(CkIndex_Main::q_exit(), proxy_main));
      } else {
	ERROR2 ("CommBlock::control_sync()",  
		"Unknown phase: phase %s sync type %s", 
		phase_string[phase],sync.c_str());    
      }

    }

  } else if (sync == "neighbor") {

    control_sync_neighbor_(phase);

  } else if (sync == "array") {

    if (phase == sync_adapt_enter) {
      if (index().is_root()) thisProxy.p_adapt_enter();
    } else if (phase == sync_adapt_next) {
      if (index().is_root()) thisProxy.p_adapt_next();
    } else if (phase == sync_adapt_called) {
      if (index().is_root()) thisProxy.p_adapt_called();
    } else if (phase == sync_adapt_end) {
      if (index().is_root()) thisProxy.p_adapt_end();
    } else if (phase == sync_adapt_exit) {
      if (index().is_root()) thisProxy.p_adapt_exit();
    } else if (phase == sync_refresh_enter) {
      if (index().is_root()) thisProxy.p_refresh_enter();
    } else if (phase == sync_refresh_exit) {
      if (index().is_root()) thisProxy.p_refresh_exit();
    } else if (phase == sync_output_enter) {
      if (index().is_root()) thisProxy.p_output_enter();
    } else if (phase == sync_output_exit) {
      if (index().is_root()) thisProxy.p_output_exit();
    } else if (phase == sync_compute_enter) {
      if (index().is_root()) thisProxy.p_compute_enter();
    } else if (phase == sync_compute_exit) {
      if (index().is_root()) thisProxy.p_compute_exit();
    } else if (phase == sync_stopping_enter) {
      if (index().is_root()) thisProxy.p_stopping_enter();
    } else if (phase == sync_stopping_exit) {
      if (index().is_root()) thisProxy.p_stopping_exit();
    } else if (phase == sync_exit) {
      if (index().is_root()) thisProxy.p_exit();
    } else {
      ERROR2 ("CommBlock::control_sync()",  
	      "Unknown phase: phase %s sync type %s", 
	      phase_string[phase],sync.c_str());    
    }

  } else if (sync == "none") {

    control_call_phase_(phase);

  } else {
    ERROR2 ("CommBlock::control_sync()",  
	    "Unknown sync type: phase %s sync type %s", 
	    phase_string[phase],sync.c_str());    
  }
}

//----------------------------------------------------------------------

void CommBlock::control_sync_count_ (int phase, int count)
{
  if (count != 0)  max_sync_[phase] = count;

  ++count_sync_[phase];
#ifdef CELLO_DEBUG
  fprintf (simulation()->fp_debug(),"%d %s counting %d/%d leaf %d %s\n",__LINE__,name_.c_str(),
	   count_sync_[phase],max_sync_[phase],is_leaf(),phase_string[phase]);
#endif

  if (0 < max_sync_[phase] && max_sync_[phase] <= count_sync_[phase]) {
    max_sync_[phase] = 0;
    count_sync_[phase] = 0;
#ifdef CELLO_DEBUG
  fprintf (simulation()->fp_debug(),"%d %s called %s\n",__LINE__,name_.c_str(),phase_string[phase]);
#endif
    control_call_phase_(phase);
  }
}

//----------------------------------------------------------------------

void CommBlock::control_sync_neighbor_(int phase)
{
  if (!is_leaf()) {

#ifdef CELLO_DEBUG
  fprintf (simulation()->fp_debug(),"%d %s called %s\n",__LINE__,name_.c_str(),phase_string[phase]);
#endif
    control_call_phase_ (phase);

  } else {

    const int level        = this->level();
    const int rank = this->rank();
    const int rank_refresh = simulation()->config()->field_refresh_rank;

    ItFace it_face = this->it_face();

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
	  fprintf (simulation()->fp_debug(),"%d %s calling p_control_sync phase %s leaf %d block %s\n",
	    __LINE__,name_.c_str(), phase_string[phase], is_leaf(),index_neighbor.bit_string(-1,2).c_str());
#endif
	thisProxy[index_neighbor].p_control_sync_count(phase,0);

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
	  fprintf (simulation()->fp_debug(),"%d %s calling p_control_sync phase %s leaf %d block %s\n",
	    __LINE__,name_.c_str(), phase_string[phase], is_leaf(),index_uncle.bit_string(-1,2).c_str());
#endif
	  thisProxy[index_uncle].p_control_sync_count(phase,0);

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
	  fprintf (simulation()->fp_debug(),"%d %s calling p_control_sync phase %s leaf %d block %s\n",
	    __LINE__,name_.c_str(), phase_string[phase], is_leaf(),index_nibling.bit_string(-1,2).c_str());
#endif
	  thisProxy[index_nibling].p_control_sync_count(phase,0);

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
#ifdef CELLO_DEBUG
  fprintf (simulation()->fp_debug(),"%d %s called %s\n",__LINE__,name_.c_str(),phase_string[phase]);
#endif
  if (phase == sync_adapt_called) {
    adapt_called_() ;
  } else if (phase == sync_adapt_enter) {
    adapt_enter_() ;
  } else if (phase == sync_adapt_next) {
    adapt_next_() ;
  } else if (phase == sync_adapt_end) {
    adapt_end_() ;
  } else if (phase == sync_adapt_exit) {
    adapt_exit_() ;
  } else if (phase == sync_refresh_enter) {
    refresh_enter_() ;
  } else if (phase == sync_refresh_exit) {
    refresh_exit_();
  } else if (phase == sync_output_enter) {
    output_enter_() ;
  } else if (phase == sync_output_exit) {
    output_exit_();
  } else if (phase == sync_compute_enter) {
    compute_enter_() ;
  } else if (phase == sync_compute_exit) {
    compute_exit_();
  } else if (phase == sync_stopping_enter) {
    stopping_enter_() ;
  } else if (phase == sync_stopping_exit) {
    stopping_exit_();
  } else if (phase == sync_exit) {
    exit_();
  } else {
    ERROR1 ("CommBlock::control_call_phase_()",  
	    "Unknown phase: phase %s", 
	    phase_string[phase]);    
  }
}
//======================================================================

