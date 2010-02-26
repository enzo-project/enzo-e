//
// $Id$
//
// See LICENSE_CELLO file for license and copyright information
//

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

///
/// @brief     
/// @author    
/// @date      
/// @ingroup
/// @note      
///

/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */


/** 
 *********************************************************************
 *
 * @file      simulation.hpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @note      
 *
 * 
 *
 * $Id$
 *
 *********************************************************************
 */

#include <vector>

#include "cello.h"

#include "parameters.hpp"
#include "amr.hpp"
#include "field.hpp"
#include "method.hpp"
#include "parallel.hpp"

class Simulation {

/** 
 *********************************************************************
 *
 * @class     Simulation
 * @brief     
 * @ingroup   Simulation
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
  Simulation();

  /// 
  ~Simulation();

  /// Initialize from parameters
  void initialize (Parameters * parameters);

  /// Return the dimension 1 <= d <= 3, of the simulation
  int dimension() 
  {
    assert (lower_.size() == upper_.size()); 
    return lower_.size(); 
  };

  /// Return extents
  void domain_extents (int lower[], int upper[])
  {
    lower[0] = lower_[0];
    lower[1] = lower_[1];
    lower[2] = lower_[2];
    upper[0] = upper_[0];
    upper[1] = upper_[1];
    upper[2] = upper_[2];
  };

  /// Return the number of fields
  int num_fields () { return fields_.size(); };

  /// Return the ith field
  Field * field (int i) { return fields_[i]; };

  /// Return the number of methods
  int num_methods () { return methods_.size(); };

  /// Return the ith method
  Method * method (int i) { return methods_[i]; };

  /// Return the parallel object
  Parallel * parallel () { return parallel_; };

  /// Return the Amr object
  Amr * amr () { return amr_; };

private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// 

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Domain extents
  std::vector<double> lower_;
  std::vector<double> upper_;

  /// AMR grid
  Amr * amr_;

  /// Parallelization
  Parallel * parallel_;

  /// Fields
  std::vector<Field *> fields_;

  /// Methods
  std::vector<Method *> methods_;

};

#endif /* SIMULATION_HPP */

