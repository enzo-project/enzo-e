// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo PPM hydro method

#ifndef ENZO_ENZO_METHOD_PPM_HPP
#define ENZO_ENZO_METHOD_PPM_HPP

class EnzoMethodPpm : public Method {

  /// @class    EnzoMethodPpm
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's PPM hydro method

public: // interface

  /// Create a new EnzoMethodPpm object
  EnzoMethodPpm(EnzoConfig * enzo_config);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodPpm);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodPpm (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "ppm"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // interface

  int comoving_coordinates_;
};

#endif /* ENZO_ENZO_METHOD_PPM_HPP */
