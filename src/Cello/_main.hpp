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

#ifdef CONFIG_USE_CHARM

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
  
  /// Exit the program
  void p_exit(int count);

  /// Finalize the simulation
  void enzo_finalize(Simulation * simulation);

private:

   int count_exit_; 
   Monitor * monitor_;

};
#endif
