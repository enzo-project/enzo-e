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

#ifdef CONFIG_USE_GRACKLE

class EnzoMethodGrackle : public Method {

  /// @class    EnzoMethodGrackle
  /// @ingroup  Enzo
///
/// This class interfaces the Grackle primordial chemistry / cooling
/// library with Cello

public: // interface

  /// Create a new EnzoMethodGrackle object
  EnzoMethodGrackle(EnzoConfig *, const FieldDescr * field_descr);

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

  void initialize_chemistry_(EnzoConfig * c);
  void initialize_units_(EnzoConfig * c);
  void initialize_fields_(const FieldDescr * field_descr);

protected: // attributes

  /// Grackle struct defining code units
  code_units units_;

  /// Grackle struct defining chemistry data
  chemistry_data chemistry_;

  /// Field id's
  std::map<std::string,int> field_;
};

#endif /* CONFIG_USE_GRACKLE */

#endif /* ENZO_ENZO_METHOD_GRACKLE_HPP */
