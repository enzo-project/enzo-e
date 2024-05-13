// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodCloseFiles.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    

#include <chrono>
#include <thread>

#include "problem.hpp"

// #define DEBUG_THROTTLE

//----------------------------------------------------------------------

CmiNodeLock MethodCloseFiles::node_lock;

void method_close_files_mutex_init()
{
  MethodCloseFiles::node_lock = CmiCreateLock();
}
//----------------------------------------------------------------------

MethodCloseFiles::MethodCloseFiles(ParameterGroup p) throw()
  : Method(),
    seconds_stagger_( p.value_float("seconds_stagger",0.0) ),
    seconds_delay_( p.value_float("seconds_delay",0.0) ),
    group_size_(p.value_integer("group_size", std::numeric_limits<int>::max()))
{
  cello::simulation()->refresh_set_name(ir_post_,"close_files");
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();
}

//----------------------------------------------------------------------

void MethodCloseFiles::compute( Block * block) throw()
{
#ifdef CONFIG_SMP_MODE
  const bool is_first_cycle = (block->cycle() == cello::config()->initial_cycle);
  if (is_first_cycle) {
    throttle_stagger_();
    CmiLock(MethodCloseFiles::node_lock);
    for (auto it=FileHdf5::file_list.begin();
         it!=FileHdf5::file_list.end(); ++it) {
#ifdef DEBUG_THROTTLE
      CkPrintf ("%d %g DEBUG_THROTTLE closed %s\n",
                CkMyPe(),cello::simulation()->timer(),it->first.c_str());
      fflush(stdout);
#endif
      it->second->file_close();
      delete it->second;
      throttle_delay_();
      FileHdf5::file_list.erase(it);
    }
    CmiUnlock(MethodCloseFiles::node_lock);
  }

  block->compute_done(); 
#else
  ERROR("MethodCloseFiles::compute()",
        "\"close_files\" method can hang if CONFIG_SMP_MODE is not defined (smp=0)");
#endif

}

//----------------------------------------------------------------------

void MethodCloseFiles::throttle_stagger_()
{
  if (seconds_stagger_ > 0.0) {
    //--------------------------------------------------
    static int count_threads = 0;
    const int node_size = CkNumPes() / CkNumNodes();
    if ((seconds_stagger_ > 0.0) &&
	count_threads < node_size ) {
      ++count_threads;
      int ms = 1000*((CkMyPe() / node_size ) % group_size_) * seconds_stagger_;
#ifdef DEBUG_THROTTLE  
      CkPrintf ("%d %g DEBUG_THROTTLE %d ms stagger start\n",
		CkMyPe(),cello::simulation()->timer(),ms);
      fflush(stdout);
#endif      
      std::this_thread::sleep_for(std::chrono::milliseconds(ms));
#ifdef DEBUG_THROTTLE  
      CkPrintf ("%d %g DEBUG_THROTTLE %d ms stagger stop\n",
		CkMyPe(),cello::simulation()->timer(),ms);
      fflush(stdout);
#endif      
    }
  }
}

//----------------------------------------------------------------------

void MethodCloseFiles::throttle_delay_()
{
  if (seconds_delay_ > 0.0) {
    int ms = 1000*seconds_delay_;
#ifdef DEBUG_THROTTLE  
    CkPrintf ("%d %g DEBUG_THROTTLE %d ms delay\n",
	      CkMyPe(),cello::simulation()->timer(),ms);
    fflush(stdout);
#endif      
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
  }
}

//======================================================================

