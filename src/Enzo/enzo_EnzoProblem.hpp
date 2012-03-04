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

  /// Copy constructor
  EnzoProblem(const EnzoProblem & EnzoProblem) throw();

  /// Assignment operator
  EnzoProblem & operator= (const EnzoProblem & EnzoProblem) throw();

private: // functions

  /// Create named boundary conditions object
  virtual Boundary * create_boundary_ (std::string name) throw ();

private: // attributes


};

#endif /* ENZO_ENZO_PROBLEM_HPP */

