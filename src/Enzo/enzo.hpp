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
#include "_view.hpp"
#include "_data.hpp"
#include "_simulation.hpp"
#include "_disk.hpp"
#include "_io.hpp"
#include "_compute.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

class CProxy_EnzoBlock;
class EnzoBlock;
class EnzoConfig;
class EnzoFactory;
class EnzoProblem;
class EnzoSimulation;
class EnzoUnits;
class GrackleChemistryData;

class EnzoPhysicsCosmology;
class EnzoPhysicsFluidProps;
class EnzoMethodGrackle;

/// Namespace for Enzo global accessor functions
namespace enzo {

  CProxy_EnzoBlock          block_array();
  EnzoBlock *               block ( Block * block);
  const EnzoConfig *        config();
  const EnzoFactory *       factory();
  EnzoProblem *             problem();
  EnzoSimulation *          simulation();
  EnzoUnits *               units();

  EnzoPhysicsCosmology *    cosmology();
  EnzoPhysicsFluidProps *   fluid_props();
  const EnzoMethodGrackle * grackle_method();

  /// Returns a pointer of GrackleChemistryData, if grackle is being used by
  /// the simulation, otherwise it returns nullptr.
  ///
  /// For a returnd value, `ret`, it's safe to assume that when
  /// `ret != nullptr` that `ret->get<int>("use_grackle") == 1`.
  const GrackleChemistryData * grackle_chemistry();

  /// returns the gravitational constant in code units
  ///
  /// @note
  /// One might naively assume that the gravitational constant is
  /// time-dependent when written in cosmological code units (since they are
  /// comoving). However, the cosmological code units are explicitly defined so
  /// that the gravitational constant is always fixed in code units
  double grav_constant_codeU() noexcept;

  /// returns the gravitational constant in cgs
  ///
  /// This is generally preferable to enzo_constants::standard_grav_constant
  /// since this will return the user-customizable value of the gravitational
  /// constant (the user can only customize the value in non-cosmological sims)
  double grav_constant_cgs() noexcept;
}

// this include statement must follow the above function declarations, so that
// they can be used within header files
#include "_enzo.hpp"

extern CProxy_EnzoSimulation proxy_enzo_simulation;
extern CProxy_IoEnzoWriter proxy_io_enzo_writer;
extern CProxy_IoEnzoReader proxy_io_enzo_reader;
extern CProxy_EnzoLevelArray proxy_level_array;
extern void mutex_init();
extern void mutex_init_bcg_iter();
#endif /* ENZO_HPP */

