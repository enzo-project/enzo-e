// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef SCHEDULE_HPP
#define SCHEDULE_HPP

/// @file     schedule.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-10 16:12:58
/// @brief    Interface for the Schedule class

#include <stdio.h>

#include "simulation.hpp"
#include "parameters.hpp"

class Schedule {

  /// @class    Schedule
  /// @ingroup  Schedule
  /// @brief    Class for controlling scheduling of parallel tasks

public: // interface

  /// Initialize a Schedule object
  Schedule();

  /// Delete a Schedule object
  ~Schedule();

  /// Create the simulation
  void create_simulation();

  /// Initialize the simulation
  void initialize_simulation();

  /// Run the simulation
  void execute_simulation();

  /// Terminate the simulation
  void terminate_simulation();

private: // attributes

  /// The simulation we're scheduling
  Simulation * simulation_;

};

#endif /* SCHEDULE_HPP */

