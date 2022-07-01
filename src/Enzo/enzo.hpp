// See LICENSE_CELLO file for license and copyright information

/// @file     enzo.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Include file for the \ref Enzo component

#ifndef ENZO_HPP
#define ENZO_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdio.h>

#include <vector>
#include <string>
#include <limits>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "cello.hpp"
#include "charm.hpp"

#include "_control.hpp"
#include "_error.hpp"
#include "_monitor.hpp"
#include "_parallel.hpp"
#include "_memory.hpp"
#include "_parameters.hpp"
#include "_performance.hpp"
#include "_problem.hpp"
#include "_mesh.hpp"
#include "_array.hpp"
#include "_data.hpp"
#include "_simulation.hpp"
#include "_disk.hpp"
#include "_io.hpp"
#include "_compute.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "_enzo.hpp"

class CProxy_EnzoBlock;
class EnzoConfig;
class EnzoFactory;
class EnzoPhysicsCosmology;
class EnzoProblem;
class EnzoSimulation;
class EnzoUnits;

/// Namespace for Enzo global constants and accessor functions
namespace enzo {


  const EnzoConfig *        config();
  const EnzoFactory *       factory();
  EnzoProblem *             problem();
  EnzoSimulation *          simulation();
  EnzoPhysicsCosmology *    cosmology();
  EnzoPhysicsFluidProps *   fluid_props();

  const EnzoMethodGrackle * grackle_method();

  CProxy_EnzoBlock          block_array();
  EnzoBlock *               block ( Block * block);
  EnzoPhysicsCosmology *    cosmology();
  EnzoProblem *             problem();
  EnzoSimulation *          simulation();
  EnzoUnits *               units();

  /// Returns whether the dual energy formalism is in use.
  ///
  /// @param default_ret[in] The value to return if no hydro methods are used.
  ///     The default value is false.
  bool uses_dual_energy_formalism(bool default_ret = false);

}

extern CProxy_EnzoSimulation proxy_enzo_simulation;
extern CProxy_IoEnzoWriter proxy_io_enzo_writer;
extern CProxy_IoEnzoReader proxy_io_enzo_reader;
extern void mutex_init();
extern void mutex_init_bcg_iter();
#endif /* ENZO_HPP */

