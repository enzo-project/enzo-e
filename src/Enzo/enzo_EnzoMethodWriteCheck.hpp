// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodWriteCheck.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2022-02-12
/// @brief    [\ref Enzo] Implementation of EnzoMethodWriteCheck

#ifndef ENZO_ENZO_METHOD_WRITE_CHECK_HPP
#define ENZO_ENZO_METHOD_WRITE_CHECK_HPP

class EnzoMethodWriteCheck : public Method {

  /// @class    EnzoMethodWriteCheck
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Write Enzo-E Checkpoint files

public: // interface

  /// Create a new EnzoMethodWriteCheck object
  EnzoMethodWriteCheck();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodWriteCheck);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodWriteCheck (CkMigrateMessage *m)
    : Method (m)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual methods

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "check"; }

protected: // interface

};

#endif /* ENZO_ENZO_METHOD_WRITE_CHECK_HPP */
