// See LICENSE_CELLO file for license and copyright information

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

/// @file     simulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 11 17:20:03 PST 2010
/// @brief    Include file for the \ref Simulation component

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdio.h>
#include <vector>

class Factory;

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "performance.hpp"
#include "monitor.hpp"
#include "disk.hpp"

#include "mesh.hpp"
#include "method.hpp"
#include "parameters.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "simulation_Simulation.hpp"

// Output scheduling class
#include "simulation_Schedule.hpp"

// I/O Readers/Writers
#include "simulation_Io.hpp"
#include "simulation_IoSimulation.hpp"

// Input/Output classes
#include "simulation_Output.hpp"
#include "simulation_OutputImage.hpp"



#endif /* SIMULATION_HPP */

