// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoControl.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-01-12
/// @brief    Enzo-dependent functions controling control flow of charm entry functions
/// @ingroup  Enzo

#include "enzo.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

//======================================================================

void EnzoBlock::control_sync
(int phase, std::string sync, bool next_phase, const char * file, int line)
{

  
  if (phase < phase_enzo_first) {
    CommBlock::control_sync(phase,sync,next_phase,file,line);
    return;
  } 

  if (phase == phase_enzo_matvec) {

    if (sync == "contribute") {

      CkCallback cb;

      cb = CkCallback (CkIndex_EnzoBlock::r_enzo_matvec(NULL), thisProxy);
  
      contribute(cb);

    } else if (sync == "quiescence") {

      // Quiescence detection through the Main chare

      if (index_.is_root()) {
	CkStartQD(CkCallback(CkIndex_Main::p_enzo_matvec(), proxy_main));

      }

    } else if (sync == "neighbor") {

      control_sync_neighbor_(phase);

    } else if (sync == "array") {

      if (index().is_root()) ((CProxy_EnzoBlock)thisProxy).p_enzo_matvec();

    } else if (sync == "none") {

      control_call_phase_(phase);

    }
  }
}

//----------------------------------------------------------------------

void EnzoBlock::control_call_phase_ (int phase)
{

  if (phase < phase_enzo_first) {
    CommBlock::control_call_phase_(phase);
    return;
  }

#ifdef CELLO_DEBUG
  fprintf (simulation()->fp_debug(),"%d %s called %s\n",__LINE__,name_.c_str(),phase_name[phase]);
#endif
  
  if (phase == phase_enzo_matvec) {
    enzo_matvec_();
  }
}

//======================================================================

