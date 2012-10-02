// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoTimestep.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Enzo] Declaration for the EnzoTimestep component

#ifndef ENZO_ENZO_TIMESTEP_HPP
#define ENZO_ENZO_TIMESTEP_HPP

class EnzoTimestep : public Timestep {

  /// @class    EnzoTimestep
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate determination of timestep

public: // interface

  /// Create a new EnzoTimestep
  EnzoTimestep() throw();

#ifdef CONFIG_USE_CHARM

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoTimestep);
  
  /// Charm++ PUP::able migration constructor
  EnzoTimestep (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
#endif

public: // virtual functions

  /// Evaluate the timestep for the block

  virtual double evaluate
  ( const FieldDescr * field_descr,  Block * block ) throw(); 

};

#endif /* ENZO_ENZO_TIMESTEP_HPP */
