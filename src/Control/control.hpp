#ifndef CONTROL_HPP
#define CONTROL_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      control.hpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @bug       
 * @note      
 *
 * 
 *
 * $Id$
 *
 *********************************************************************
 */

#include <stdio.h>

#include "simulation.hpp"
#include "monitor.hpp"
#include "parameters.hpp"

class Control {

/** 
 *********************************************************************
 *
 * @class     Control
 * @brief     
 * @ingroup   Control
 *
 * 
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// 
  Control(Monitor * monitor);

  /// 
  ~Control();

  /// Read parameter file
  void read_parameters(FILE * fp);

  /// Create the simulation
  void create_simulation();

  /// Initialize the simulation
  void initialize_simulation();

  /// Run the simulation
  void execute_simulation();

  /// Terminate the simulation
  void terminate_simulation();

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// 

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// The simulation we're controlling
  Monitor * monitor_;

  /// Input parameters
  Parameters * parameters_;

  /// The simulation we're controlling
  Simulation * simulation_;


};

#endif /* CONTROL_HPP */

