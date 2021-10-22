// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpmlIG.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon May 17 14:16:01 PDT 2010
/// @author   Alexei Kritsuk (akritsuk@ucsd.edu) 
/// @date     Tue Sep 18 14:16:01 PDT 2018
/// @brief    [\ref Enzo] Implementation of Enzo PPML method for ideal gas

#ifndef ENZO_ENZO_METHOD_PPMLIG_HPP
#define ENZO_ENZO_METHOD_PPMLIG_HPP

class EnzoMethodPpmlIG : public Method {

/// @class    EnzoMethodPpmlIG
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulate Enzo's PPML MHD method for ideal gas

public: // interface

  /// Creae a new EnzoMethodPpmlIG object
  EnzoMethodPpmlIG();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodPpmlIG);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodPpmlIG (CkMigrateMessage *m)
    : Method (m),
      comoving_coordinates_(false),
      gamma_(0.0),
      b0_()
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block ) throw(); 

  virtual std::string name () throw () 
  { return "ppml_ig"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) throw();

protected: // interface

  bool comoving_coordinates_;
  enzo_float gamma_;
  enzo_float b0_[3];
  
};

#endif /* ENZO_ENZO_METHOD_PPMLIG_HPP */
