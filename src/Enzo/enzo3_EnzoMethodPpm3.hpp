// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm3.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Jul 17 13:20:28 PDT 2014
/// @brief    [\ref Enzo] Implementation of Enzo 3.0 PPM hydro method

#ifndef ENZO_ENZO_METHOD_PPM3_HPP
#define ENZO_ENZO_METHOD_PPM3_HPP

class EnzoMethodPpm3 : public Method {

  /// @class    EnzoMethodPpm3
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo 3.0's PPM hydro method

public: // interface

  /// Create a new EnzoMethodPpm3 object
  EnzoMethodPpm3();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodPpm3);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodPpm3 (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( CommBlock * comm_block) throw();

  /// Compute maximum timestep for this method
  virtual double timestep ( CommBlock * comm_block) throw();

};

#endif /* ENZO_ENZO_METHOD_PPM3_HPP */
