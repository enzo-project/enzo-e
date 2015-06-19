// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeAcceleration.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Enzo] Implementation of Enzo's ComputeAcceleration functions

#ifndef ENZO_ENZO_COMPUTE_ACCELERATION_HPP
#define ENZO_ENZO_COMPUTE_ACCELERATION_HPP

class EnzoComputeAcceleration : public Compute {

  /// @class    EnzoComputeAcceleration
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ComputeAcceleration functions

public: // interface

  /// Create a new EnzoComputeAcceleration object
  EnzoComputeAcceleration(FieldDescr * field_descr,
			  int rank, bool symmetric, int order);

  /// Create a new EnzoComputeAcceleration object
  EnzoComputeAcceleration() : rank_(0), symmetric_(false), order_(0) {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoComputeAcceleration);
  
  /// Charm++ PUP::able migration constructor
  EnzoComputeAcceleration (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Perform the computation on the block
  virtual void compute( Block * block) throw();

protected: // functions

  template <typename T>
  void compute_(Block * block);

protected: // attributes

  /// Dimensionality of the problem: 1, 2, or 3
  int rank_;

  /// Whether to use symmetric differencing or one-sided differencing for 
  /// 2nd order
  bool symmetric_;

  /// Order of the differencing, either 2, 4, or 6
  int order_;

  /// Field ID's
  int i_ax_,i_ay_,i_az_,i_p_;
};

#endif /* ENZO_ENZO_COMPUTE_ACCELERATION_HPP */
