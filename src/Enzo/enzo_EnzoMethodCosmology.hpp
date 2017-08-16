// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodCosmology.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-07-24
/// @brief    [\ref Enzo] Declaration of the EnzoMethodCosmology class

#ifndef ENZO_ENZO_METHOD_COSMOLOGY_HPP
#define ENZO_ENZO_METHOD_COSMOLOGY_HPP

class EnzoMethodCosmology : public Method {

  /// @class    EnzoMethodCosmology
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoMethodCosmology(const FieldDescr *) throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodCosmology);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodCosmology (CkMigrateMessage *m)
    : Method (m)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { }
  
public: // virtual methods
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "cosmology"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw()
  { return std::numeric_limits<double>::max(); }

private: // methods


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* ENZO_ENZO_METHOD_COSMOLOGY_HPP */

