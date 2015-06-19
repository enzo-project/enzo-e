// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpml.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon May 17 14:16:01 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo PPML method

#ifndef ENZO_ENZO_METHOD_PPML_HPP
#define ENZO_ENZO_METHOD_PPML_HPP

class EnzoMethodPpml : public Method {

/// @class    EnzoMethodPpml
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulate Enzo's PPML hydro method

public: // interface

  /// Creae a new EnzoMethodPpml object
  EnzoMethodPpml(const FieldDescr * field_descr,EnzoConfig * enzo_config);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodPpml);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodPpml (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block ) throw(); 

  virtual std::string name () throw () 
  { return "ppml"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // interface

  int comoving_coordinates_;
};

#endif /* ENZO_ENZO_METHOD_PPML_HPP */
