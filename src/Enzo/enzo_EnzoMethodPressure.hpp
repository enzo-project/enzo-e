// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPressure.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Wed Jul 23 15:53:28 PDT 2014
/// @brief    [\ref Enzo] Implementation of Enzo's ComputePressure functions

#ifndef ENZO_ENZO_METHOD_PRESSURE_HPP
#define ENZO_ENZO_METHOD_PRESSURE_HPP

class EnzoMethodPressure : public Method {

  /// @class    EnzoMethodPressure
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ComputePressure functions

public: // interface

  /// Create a new EnzoMethodPressure object
  EnzoMethodPressure(double gamma);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodPressure);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodPressure (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( CommBlock * comm_block) throw();

private: // functions

  template <typename T>
  void compute_(CommBlock * comm_block);

private: // attributes

  double gamma_;

};

#endif /* ENZO_ENZO_METHOD_PRESSURE_HPP */
