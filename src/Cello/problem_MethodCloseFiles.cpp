// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodCloseFiles.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    

#include <chrono>
#include <thread>

#include "problem.hpp"

#define DEBUG_THROTTLE

//----------------------------------------------------------------------

CmiNodeLock MethodCloseFiles::node_lock;

void method_close_files_mutex_init()
{
  MethodCloseFiles::node_lock = CmiCreateLock();
}
//----------------------------------------------------------------------

void MethodCloseFiles::compute( Block * block) throw()
{
#ifdef CONFIG_SMP_MODE
  const bool is_first_cycle = (block->cycle() == cello::config()->initial_cycle);
  if (is_first_cycle) {
    CmiLock(MethodCloseFiles::node_lock);
    for (auto it=FileHdf5::file_list.begin();
         it!=FileHdf5::file_list.end(); ++it) {
#ifdef DEBUG_THROTTLE
      CkPrintf ("%d %g DEBUG_THROTTLE closed %s\n",
                CkMyPe(),cello::simulation()->timer(),it->first.c_str());
      fflush(stdout);
#endif
      it->second->file_close();
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

void MethodCloseFiles::throttle_delay_()
{
  if (seconds_delay_ > 0.0) {
    int ms = 1000*seconds_delay_;
#ifdef DEBUG_THROTTLE  
    CkPrintf ("DEBUG_THROTTLE %d %g %d ms delay\n",
	      CkMyPe(),cello::simulation()->timer(),ms);
    fflush(stdout);
#endif      
    std::this_thread::sleep_for(std::chrono::milliseconds(ms));
  }
}

//======================================================================

