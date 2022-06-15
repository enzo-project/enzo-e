// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodCheckpoint.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2021-03-09
/// @brief    [\ref Problem] Declaration for the MethodCheckpoint class

#ifndef PROBLEM_METHOD_CHECKPOINT_HPP
#define PROBLEM_METHOD_CHECKPOINT_HPP

class MethodCheckpoint : public Method
{
  /// @class    MethodCheckpoint
  /// @ingroup  MethodCheckpoint
  /// @brief    [\ref MethodCheckpoint] Declaration of MethodCheckpoint
  ///
  /// Method for writing checkpoints to disk via Charm++.

public: // interface

  /// Create a new MethodCheckpoint
  MethodCheckpoint ( std::vector< std::string > path_name);

  /// Destructor
  virtual ~MethodCheckpoint() throw()
  { };

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodCheckpoint);

  /// Charm++ PUP::able migration constructor
  MethodCheckpoint (CkMigrateMessage *m)
    : Method(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  void compute_continue (Block * block);

public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute ( Block * block) throw();

  /// Return the name of this MethodCheckpoint
  virtual std::string name () throw ()
  { return "checkpoint"; }

protected: // attributes

  /// Path name and format
  std::vector <std::string> path_name_;
};


#endif /* PROBLEM_METHOD_CHECKPOINT_HPP */
