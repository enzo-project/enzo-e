// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGrackle.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu May 15 14:32:28 EDT 2014
/// @brief    [\ref Enzo] Declaration of EnzoMethodGrackle class
///
/// This class interfaces the Grackle primordial chemistry / cooling
/// library with Cello

#ifndef ENZO_ENZO_METHOD_GRACKLE_HPP
#define ENZO_ENZO_METHOD_GRACKLE_HPP

class EnzoMethodGrackle : public Method {

  /// @class    EnzoMethodGrackle
  /// @ingroup  Enzo
///
/// This class interfaces the Grackle primordial chemistry / cooling
/// library with Cello

public: // interface

  /// Create a new EnzoMethodGrackle object
  EnzoMethodGrackle(EnzoConfig *);

  /// Create a new EnzoMethodGrackle object
  EnzoMethodGrackle() : Method() {};

  /// Destructor
  virtual ~EnzoMethodGrackle() throw() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGrackle);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGrackle (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( CommBlock * comm_block) throw();

  /// Compute maximum timestep for this method
  virtual double timestep ( CommBlock * comm_block) throw();

protected: // methods

protected: // attributes

  /// Grackle struct defining chemistry data
  chemistry_data chemistry_;
};

#endif /* ENZO_ENZO_METHOD_GRACKLE_HPP */
