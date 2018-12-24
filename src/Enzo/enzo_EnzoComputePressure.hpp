// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputePressure.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Enzo] Implementation of Enzo's ComputePressure functions

#ifndef ENZO_ENZO_COMPUTE_PRESSURE_HPP
#define ENZO_ENZO_COMPUTE_PRESSURE_HPP

class EnzoComputePressure : public Compute {

  /// @class    EnzoComputePressure
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ComputePressure functions

public: // interface

  /// Create a new EnzoComputePressure object
  EnzoComputePressure(double gamma,
		      bool comoving_coordinates);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoComputePressure);

  /// Charm++ PUP::able migration constructor
  EnzoComputePressure (CkMigrateMessage *m)
    : Compute(m),
      gamma_(0.0),
      comoving_coordinates_(false)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Perform the computation on the block
  virtual void compute( Block * block) throw();

  void compute_(Block * block
#ifdef CONFIG_USE_GRACKLE
               , code_units * grackle_units = NULL,
                 grackle_field_data * grackle_fields = NULL
#endif
    );

    

protected: // attributes

  double gamma_;
  bool comoving_coordinates_;

};

#endif /* ENZO_ENZO_COMPUTE_PRESSURE_HPP */
