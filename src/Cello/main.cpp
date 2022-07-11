// See LICENSE_CELLO file for license and copyright information

/// @file     main.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-08-10
/// @brief    Implementation of main-level CHARM entry functions in main.ci

#include <boost/filesystem.hpp>
#include "cello.hpp"
#include "test.hpp"
#include "parallel.hpp"
#include "monitor.hpp"
#include "main.hpp"

//----------------------------------------------------------------------

#if defined(CHARM_SIMULATION) || defined(CHARM_MESH)
// crude hack to provide dummy implementations of restart-related functions
// that are normally defined in enzo_control_restart.cpp.
//
// This is currently necessary for defining certain classes of test problems
// without including the entire enzo-layer in the binary

#include "simulation.hpp"
void Block::restart_enter_()
{ ERROR("Block::restart_enter", "invalid call"); }
void Simulation::p_restart_enter (std::string name_dir)
{ ERROR("Simulation::p_restart_enter", "invalid call"); }
void Simulation::r_restart_start (CkReductionMsg * msg)
{ ERROR("Simulation::p_restart_start", "invalid call"); }
#endif

//----------------------------------------------------------------------

CProxy_Main proxy_main;

// #define DEBUG_MAIN

#ifdef CHARM_ENZO
#include "simulation.hpp"
extern CProxy_EnzoSimulation proxy_simulation;
#include "enzo_finalize.hpp"
#endif

#ifdef DEBUG_MAIN
#  define TRACE_MAIN(MSG) CkPrintf ("DEBUG_MAIN %s\n",MSG); fflush(stdout);
#else
#  define TRACE_MAIN(MSG) /* ... */
#endif

//----------------------------------------------------------------------

void Main::p_exit(int count)
{
  DEBUG("Main::p_exit");
  count_exit_++;
  unit_finalize();
  if (count_exit_ >= count) {
    count_exit_ = 0;
    exit_();
  }
}

void Main::exit_()
{

#ifdef CHARM_ENZO

  EnzoSimulation * simulation = enzo::simulation();

  if (simulation) {
    enzo_finalize(simulation);
  }

#endif

  if (Monitor::instance()) {
    Monitor::instance()->print ("","END CELLO");
  }

  PARALLEL_EXIT;
}


//----------------------------------------------------------------------

void Main::p_checkpoint_output(int count, std::string dir_name)
{
  
  TRACE_MAIN("DEBUG MAIN p_checkpoint_output");
  
  count_checkpoint_++;
  if (count_checkpoint_ >= count) {
    count_checkpoint_ = 0;
    // Write parameter file

#ifdef CHARM_ENZO
    strncpy(dir_checkpoint_,dir_name.c_str(),255);
    Simulation * simulation = cello::simulation();
    simulation->set_checkpoint(dir_checkpoint_);
#endif    

#ifdef CHARM_ENZO
    CkPrintf ("Calling CkStartCheckpoint\n");
    CkCallback callback(CkIndex_EnzoSimulation::r_write_checkpoint_output(),proxy_simulation);
    CkStartCheckpoint (dir_checkpoint_,callback,false,1);
    // "OLD" CHARM++ (version < 7.0.0) USE:
    //CkStartCheckpoint (dir_checkpoint_,callback);
#endif
  }
  // --------------------------------------------------
}


//----------------------------------------------------------------------

void Main::p_output_enter()
{
#ifdef CHARM_ENZO

  cello::block_array().p_output_enter();
  
#endif
}

//----------------------------------------------------------------------

void Main::p_output_exit()
{
#ifdef CHARM_ENZO
  cello::block_array().p_output_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_compute_enter()
{
#ifdef CHARM_ENZO
  cello::block_array().p_compute_enter();
#endif
}

//----------------------------------------------------------------------

void Main::p_compute_continue()
{
#ifdef CHARM_ENZO
  cello::block_array().p_compute_continue();
#endif
}

//----------------------------------------------------------------------

void Main::p_compute_exit()
{
#ifdef CHARM_ENZO
  cello::block_array().p_compute_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_stopping_enter()
{
#ifdef CHARM_ENZO
  cello::block_array().p_stopping_enter();
#endif
}

//----------------------------------------------------------------------

void Main::p_stopping_balance()
{
#ifdef CHARM_ENZO
  cello::block_array().p_stopping_load_balance();
#endif
}

//----------------------------------------------------------------------

void Main::p_stopping_exit()
{
#ifdef CHARM_ENZO
  cello::block_array().p_stopping_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_text_file_write
(int nd, char * dir,
 int nf, char * file,
 int nl, char * line, int count)
{
  // Open file if first call. Increment count by number of pe's
  std::string full_file = std::string(dir) + "/" + file;

  FILE * fp_text   = fp_text_[full_file];
  Sync * sync_text = sync_text_[full_file];
  
  if (fp_text_[full_file] == NULL) {

    // Create subdirectory if any
    boost::filesystem::path directory(dir);
    if (! boost::filesystem::is_directory(directory)) {
      ASSERT1 ("Main::p_text_file_write()",
	       "Error creating directory %s",
	       dir,
	       (boost::filesystem::create_directory(directory)));
    }
    
    fp_text   = fp_text_[full_file] = fopen(full_file.c_str(),"w");
    sync_text = sync_text_[full_file] = new Sync;
    
    sync_text->set_stop(CkNumPes());
  }

  // First call from any Block on a processor includes count
  // minus one for CkNumPes() above
  if (count > 0) {
    sync_text->inc_stop(count-1);
  }

  fprintf (fp_text,"%s",line);

  if (sync_text->next()) {
    //    text_file_close_();
    fclose(fp_text);
    fp_text_[full_file] = NULL;
    delete sync_text_[full_file];
    sync_text_[full_file] = NULL;
  }
}

//----------------------------------------------------------------------

void Main::p_exit()
{
#ifdef CHARM_ENZO
  cello::block_array().p_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_adapt_enter()
{
  TRACE_MAIN("p_adapt_enter");
#ifdef CHARM_ENZO
  cello::block_array().p_adapt_enter();
#endif
}

//----------------------------------------------------------------------

void Main::p_initial_exit()
{
#ifdef CHARM_ENZO
  cello::block_array().p_initial_exit();
#endif
}

//----------------------------------------------------------------------

void Main::p_adapt_end()
{
  TRACE_MAIN("p_adapt_end");
#ifdef CHARM_ENZO
  cello::block_array().p_adapt_end();
#endif
}

//----------------------------------------------------------------------

void Main::p_adapt_update()
{
  TRACE_MAIN("p_adapt_update");
#ifdef CHARM_ENZO
  cello::block_array().p_adapt_update();
#endif  
}

//----------------------------------------------------------------------

void Main::p_adapt_called()
{
  TRACE_MAIN("p_adapt_called");
#ifdef CHARM_ENZO
  cello::block_array().p_adapt_called();
#endif
}

//----------------------------------------------------------------------

void Main::p_adapt_exit()
{
  TRACE_MAIN("p_adapt_exit");
#ifdef CHARM_ENZO
  cello::block_array().p_adapt_exit();
#endif
}

//----------------------------------------------------------------------

#if defined(CHARM_ENZO)
#  include "main_enzo.def.h"
#elif defined(CHARM_SIMULATION)
#  include "main_simulation.def.h"
#elif defined(CHARM_MESH)
#  include "main_mesh.def.h"
#else
#  include "main.def.h"
#endif

