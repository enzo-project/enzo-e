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

#include <limits>
#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "cello.hpp"
#include "charm.hpp"

#include "_array.hpp"
#include "_compute.hpp"
#include "_control.hpp"
#include "_data.hpp"
#include "_disk.hpp"
#include "_error.hpp"
#include "_io.hpp"
#include "_memory.hpp"
#include "_mesh.hpp"
#include "_monitor.hpp"
#include "_parallel.hpp"
#include "_parameters.hpp"
#include "_performance.hpp"
#include "_problem.hpp"
#include "_simulation.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "_enzo.hpp"

class CProxy_EnzoBlock;
class EnzoConfig;
class EnzoPhysicsCosmology;
class EnzoProblem;
class EnzoSimulation;
class EnzoUnits;

/// Namespace for Enzo global constants and accessor functions
namespace enzo {
EnzoProblem *problem();
EnzoSimulation *simulation();
EnzoPhysicsCosmology *cosmology();
const EnzoMethodGrackle *grackle_method();
EnzoUnits *units();
const EnzoConfig *config();
CProxy_EnzoBlock block_array();
EnzoBlock *block(Block *block);

// Checks if given particle type exists and has an attribute
// of given name
void check_particle_attribute(std::string type, std::string attribute);
} // namespace enzo

extern CProxy_EnzoSimulation proxy_enzo_simulation;
extern void mutex_init();
extern void mutex_init_bcg_iter();
#endif /* ENZO_HPP */
