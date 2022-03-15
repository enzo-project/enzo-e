// See LICENSE_CELLO file for license and copyright information

/// @file     control_restart.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-25
/// @brief    Charm-related mesh adaptation control functions.
/// @ingroup  Control
///
/// This file controls adaptive mesh refinement on a distributed
/// array of octrees.

//--------------------------------------------------

#include "simulation.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#define TRACE_RESTART


#ifdef TRACE_RESTART
#   define TRACE_RESTART_BLOCK(MSG,BLOCK)       \
  CkPrintf ("%d TRACE_RESTART %s %s \n",        \
            CkMyPe(),BLOCK->name().c_str(),     \
            std::string(MSG).c_str());          \
  fflush(stdout);
#else
#   define TRACE_RESTART_BLOCK(MSG,BLOCK) /* ... */
#endif

//----------------------------------------------------------------------

void Block::restart_enter_()
{
  TRACE_RESTART_BLOCK("restart_enter_",this);
  const std::string restart_dir  = cello::config()->initial_restart_dir;
  if (index_.is_root()) {
    proxy_simulation[0].p_restart_enter(restart_dir);
  }
}

//----------------------------------------------------------------------

void Simulation::p_restart_enter (std::string restart_dir)
{
  // [ Called on root ip only ]

  std::ifstream file_hierarchy = file_open_hierarchy_(restart_dir);

  int num_files;
  file_hierarchy >> num_files;
  for (int i=0; i<num_files; i++) {
    std::string file_name;
    file_hierarchy >> file_name;
    CkPrintf ("TRACE_RESTART file %d %s\n",i,file_name.c_str());
  }
  // file_hierarchy << std::setfill('0');
  // int max_digits = log(check_num_files_-1)/log(10) + 1;
  // file_hierarchy << check_num_files_ << std::endl;
  // for (int i=0; i<check_num_files_; i++) {
  //   file_hierarchy << "block_data-" << std::setw(max_digits) << i << ".h5" << std::endl;
  // }

  // std::ifstream file_hierarchy = file_open_hierarchy_();
  //  ocheck.file_list 
  // open hierarchy file
  CkPrintf ("%d DEBUG_RESTART restart_dir = %s\n",CkMyPe(),restart_dir.c_str());
  // read hierarchy file
  CkPrintf("%d Main::restart_enter_ read hierarchy file\n",CkMyPe());
  // create block_array
  CkPrintf("%d Main::restart_enter_ create block array\n",CkMyPe());
  // create IoEnzoReader array
  CkPrintf("%d Main::restart_enter_ create IoEnzoReader array\n",CkMyPe());
  // initialize sync_file(num_io_reader)
  //  for (i_f = files in restart) {
  //    io_reader[i_f].insert(file_block);
  //  }
  // close hierarcy file
  CkPrintf("%d restart_enter_ close hierarchy file\n",CkMyPe());
  fflush(stdout);
}

//----------------------------------------------------------------------

std::ifstream Simulation::file_open_hierarchy_(std::string name_dir)
{
  std::string name_file = name_dir + "/check.file_list";

  std::ifstream file_hierarchy (name_file);

  ASSERT1("Simulation::file_copen_hierarchy_",
          "Cannot open hierarchy file %s for writing",
          name_file.c_str(),file_hierarchy);

  return file_hierarchy;
}

