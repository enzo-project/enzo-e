// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProblem.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Enzo] Declaration of the EnzoProblem class
///

#ifndef ENZO_ENZO_PROBLEM_HPP
#define ENZO_ENZO_PROBLEM_HPP

class EnzoProblem : public Problem {

  /// @class    EnzoProblem
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoProblem() throw();

  /// Destructor
  ~EnzoProblem() throw();

private: // functions

  /// Create named boundary conditions object
  virtual Boundary * create_boundary_ 
  (std::string name, Parameters * parameters) throw ();

  /// Create named initialization object
  virtual Initial *  create_initial_ 
  (std::string name, Parameters * parameters) throw ();

  /// Create named method object
  virtual Method *   create_method_ 
  (std::string name, Parameters * parameters) throw ();

  /// Create named timestep object
  virtual Timestep * create_timestep_
  (std::string name, Parameters * parameters) throw ();

};

#endif /* ENZO_ENZO_PROBLEM_HPP */

