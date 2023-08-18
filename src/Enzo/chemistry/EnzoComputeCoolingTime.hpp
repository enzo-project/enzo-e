// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeCoolingTime.hpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     2019-05-07
/// @brief    [\ref Enzo] Implementation of Enzo's ComputeCoolingTime functions

#ifndef ENZO_ENZO_COMPUTE_COOLINGTIME_HPP
#define ENZO_ENZO_COMPUTE_COOLINGTIME_HPP

class EnzoComputeCoolingTime : public Compute {

  /// @class    EnzoComputeCoolingTime
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ComputeCoolingTime functions

public: // interface

  /// Create a new EnzoComputeCoolingTime object
  EnzoComputeCoolingTime
  ();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoComputeCoolingTime);

  /// Charm++ PUP::able migration constructor
  EnzoComputeCoolingTime (CkMigrateMessage *m)
    : Compute(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Perform the computation on the block
  virtual void compute( Block * block) throw();

  virtual void compute( Block * block, enzo_float * ct) throw();

  // name of derived field that this function calculates
  std::string name () throw() {
    return "cooling_time";
  }

  void compute_(Block * block,
   enzo_float * ct,
   grackle_field_data * grackle_fields = NULL
 );

private: // functions


private: // attributes

};

#endif /* ENZO_ENZO_COMPUTE_COOLINGTIME_HPP */
