// See LICENSE_CELLO file for license and copyright information

/// @file     main.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-08-10
/// @brief    [\ref Main] Declaration of the Main CHARM++ chare
///
/// Different test programs require different subsets of existing chares.
/// For example, lower-level test programs, e.g. Mesh or Field tests,
/// don't need to (and shouldn't) include Enzo-specific chare declarations.
/// main.hpp is an include file used with CHARM.  Four "levels" of
/// chare's are included:
///
///     CHARM_ENZO:       include Enzo, Simulation, and Mesh chare declartions
///     CHARM_SIMULATION: include Simulation, and Mesh chare declarations
///     CHARM_MESH:       include Mesh component chare declarations
///     <default>         don't include any chare declarations

#include "cello.hpp"

#include "parallel.hpp"
#include "monitor.hpp"

class Factory;
class Simulation;

#if defined(CHARM_ENZO)

#  include "main_enzo.decl.h"

#elif defined(CHARM_SIMULATION)

#  include "main_simulation.decl.h"

#elif defined(CHARM_MESH)

#  include "main_mesh.decl.h"

#else

#  include "main.decl.h"

#endif

extern CProxy_Main proxy_main;

//----------------------------------------------------------------------

class Main : public CBase_Main {

  /// @class    Main
  /// @ingroup  Main
  /// @brief    [\ref Main] CHARM++ main chare

public:

  /// Initialize the Main chare (defined in the calling program)
  Main(CkArgMsg* main);

  /// Migration constructor
  Main(CkMigrateMessage * m)
    : CBase_Main(m),
      count_exit_(0),
      count_checkpoint_(0)
  { 
    TRACE("Main::Main(CkMigrateMessage)");
  }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    CBase_Main::pup(p);
    p|count_exit_;
    p|count_checkpoint_;
    WARNING ("Main::pup","skipping monitor_");
    if (p.isUnpacking()) monitor_ = Monitor::instance();
    //    p|*monitor_;
  }

  /// Exit the program
  void p_exit(int count);

  void p_checkpoint (int count, std::string dir_name);

  void q_adapt_enter();
  void q_adapt_called();
  void q_adapt_end();
  void q_adapt_next();
  void q_adapt_exit();
  void q_compute_enter();
  void q_compute_exit();
  void q_output_enter ();
  void q_output_exit();
  void q_refresh_enter();
  void q_refresh_exit();
  void q_stopping_enter();
  void q_stopping_exit();
  void q_exit();

  /// Finalize the simulation
  void enzo_finalize(Simulation * simulation);

private: // functions

  void exit_ ();

private: // attributes

   int count_exit_; 
   int count_checkpoint_; 
   Monitor * monitor_;

};
