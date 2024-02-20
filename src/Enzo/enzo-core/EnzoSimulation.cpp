// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    Implementation of EnzoSimulation user-dependent class member functions

/* #define CHECK_MEMORY */  /* calls mtrace() */

// #define TRACE_PARAMETERS
// #define DEBUG_ENZO_SIMULATION

#include "cello.hpp"

#include "enzo.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#ifdef CHECK_MEMORY
#   include <mcheck.h>
#endif

#include "simulation.hpp"

CProxy_EnzoSimulation proxy_enzo_simulation;
CProxy_IoEnzoWriter   proxy_io_enzo_writer;
CProxy_IoEnzoReader   proxy_io_enzo_reader;

//----------------------------------------------------------------------

EnzoSimulation::EnzoSimulation
(
 const char         parameter_file[],
 int                n)
  : CBase_EnzoSimulation(parameter_file, n),
    check_num_files_(0),
    check_ordering_(""),
    check_directory_(),
    restart_level_(0)
{
#ifdef CHECK_MEMORY
  mtrace();
#endif

#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION EnzoSimulation()\n",CkMyPe());
  fflush(stdout);
#endif

  // Synchronize to ensure all EnzoSimulation objects exist before
  // reading parameters

  CkCallback callback (CkIndex_EnzoSimulation::r_startup_begun(NULL),
                       thisProxy);
#ifdef TRACE_CONTRIBUTE
  CkPrintf ("%s:%d DEBUG_CONTRIBUTE\n",__FILE__,__LINE__); fflush(stdout);
#endif

  contribute(callback);

}

//----------------------------------------------------------------------

EnzoSimulation::~EnzoSimulation()
{
}

//----------------------------------------------------------------------

void EnzoSimulation::pup (PUP::er &p)
{
#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION EnzoSimulation::pup()\n",CkMyPe());
  fflush(stdout);
#endif  
  // NOTE: change this function whenever attributes change

  CBase_EnzoSimulation::pup(p);

  TRACEPUP;

  p | sync_check_writer_created_;
  p | sync_check_done_;
  p | check_num_files_;
  p | check_ordering_;
  p | check_directory_;
  p | restart_level_;
}

//----------------------------------------------------------------------
void EnzoSimulation::p_refine_create_block(MsgRefine * msg)
{ refine_create_block(msg); }

void EnzoSimulation::refine_create_block(MsgRefine * msg)
{
  Index index = msg->index();
  if (msg_refine_map_[index] != NULL) {
    int v3[3];
    index.values(v3);
    ASSERT3 ("EnzoSimulation::p_refine_create_block",
	    "index %08x %08x %08x is already in the msg_refine mapping",
	    v3[0],v3[1],v3[2],
	    (msg == NULL));
  }
  msg_refine_map_[index] = msg;

  int ip = CkMyPe();
  enzo::block_array()[index].insert(ip,MsgType::msg_refine,ip);
}

//----------------------------------------------------------------------

void EnzoSimulation::set_msg_check(Index index, EnzoMsgCheck * msg)
{
  if (msg_check_map_[index] != NULL) {
   
    int v3[3];
    index.values(v3);
    ASSERT3 ("EnzoSimulation::set_msg_check",
	    "index %08x %08x %08x is already in the msg_check mapping",
	    v3[0],v3[1],v3[2],
	    (msg == NULL));
  }
  msg_check_map_[index] = msg;
}

//----------------------------------------------------------------------

EnzoMsgCheck * EnzoSimulation::get_msg_check(Index index)
{
  int v3[3];
  index.values(v3);
  EnzoMsgCheck * msg = msg_check_map_[index];
  if (msg == NULL) {
    int v3[3];
    index.values(v3);
    
    ASSERT3 ("EnzoSimulation::get_msg_check",
	    "index %08x %08x %08x is not in the msg_check mapping",
	    v3[0],v3[1],v3[2],
	    (msg != NULL));
  }
  msg_check_map_.erase(index);
  return msg;
}

//----------------------------------------------------------------------

MsgRefine * EnzoSimulation::get_msg_refine(Index index)
{
  int v3[3];
  index.values(v3);
  MsgRefine * msg = msg_refine_map_[index];
  if (msg == NULL) {
    int v3[3];
    index.values(v3);
    
    ASSERT3 ("EnzoSimulation::get_msg_refine",
	    "index %08x %08x %08x is not in the msg_refine mapping",
	    v3[0],v3[1],v3[2],
	    (msg != NULL));
  }
  msg_refine_map_.erase(index);
  return msg;
}
//----------------------------------------------------------------------

void EnzoSimulation::r_startup_begun (CkReductionMsg *msg)
{

#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION r_startup_begun()\n",CkMyPe());
  fflush(stdout);
#endif  

  delete msg;

#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION r_startup_begun()\n",CkMyPe());
  fflush(stdout);
#endif  

  problem_ = new EnzoProblem;

  initialize_config_();

  initialize();
  
#ifdef TRACE_PARAMETERS
  CkPrintf ("%d END   r_startup_begun()\n",CkMyPe());
  fflush(stdout);
#endif
}

//----------------------------------------------------------------------

void EnzoSimulation::initialize_config_() throw()
{
#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION begin initialize_config()\n",CkMyPe());
  fflush(stdout);
#endif  

  config_ = static_cast<EnzoConfig*>(&g_enzo_config);

#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION end initialize_config()\n",CkMyPe());
  fflush(stdout);
#endif  

}

//----------------------------------------------------------------------

void EnzoSimulation::initialize() throw()
{
#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION initialize()\n",CkMyPe());
  fflush(stdout);
#endif  
  // Call initialize() on base Simulation class
  Simulation::initialize();

}

//----------------------------------------------------------------------

const Factory * EnzoSimulation::factory() const throw()
{ 
#ifdef DEBUG_ENZO_SIMULATION
  CkPrintf ("%d DEBUG_ENZO_SIMULATION factory()\n",CkMyPe());
  fflush(stdout);
#endif  
  if (! factory_) factory_ = new EnzoFactory;
  return factory_;
}

//======================================================================

