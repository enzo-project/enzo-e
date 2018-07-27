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
#include "_data.hpp"
#include "_problem.hpp"
#include "_mesh.hpp"
#include "_simulation.hpp"
#include "_disk.hpp"
#include "_io.hpp"
#include "_compute.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "_enzo.hpp"

class EnzoConfig;
class EnzoPhysicsCosmology;
class EnzoProblem;
class EnzoSimulation;
class EnzoUnits;

/// Namespace for Enzo global constants and accessor functions
namespace enzo {
  EnzoProblem * problem();
  EnzoSimulation * simulation();
  EnzoPhysicsCosmology * cosmology();
  EnzoUnits * units();
  const EnzoConfig * config();
};

extern CProxy_EnzoSimulation proxy_enzo_simulation;

#endif /* ENZO_HPP */

