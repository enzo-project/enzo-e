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

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoProblem);

  /// CHARM++ migration constructor
  EnzoProblem(CkMigrateMessage *m) : Problem (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

private: // functions

  /// Create named boundary conditions object
  virtual Boundary * create_boundary_ 
  (std::string type,
   int index,
   Config * config,
   Parameters * parameters
   ) throw ();

  /// Create named initialization object
  virtual Initial *  create_initial_ 
  (std::string type, 
   Parameters * parameters,
   Config * config,
   const FieldDescr *,
   const GroupProcess * group_process) throw ();

  /// Create named method object
  virtual Method *   create_method_ (std::string type) throw ();

  /// Create named timestep object
  virtual Timestep * create_timestep_
  (std::string type, Config * config) throw ();

  /// Create named interpolation object
  virtual Prolong * create_prolong_
  (std::string type, Config * config) throw ();

  /// Create named restriction object
  virtual Restrict * create_restrict_
  (std::string type, Config * config) throw ();

};

#endif /* ENZO_ENZO_PROBLEM_HPP */

