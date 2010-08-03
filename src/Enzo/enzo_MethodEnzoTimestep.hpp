// $Id: user_MethodEnzoTimestep.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     user_MethodEnzoTimestep.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    Declaration for the MethodEnzoTimestep component

#ifndef USER_METHOD_ENZO_TIMESTEP_HPP
#define USER_METHOD_ENZO_TIMESTEP_HPP

class MethodEnzoTimestep : public UserTimestep {

  /// @class    MethodEnzoTimestep
  /// @ingroup  Enzo
  /// @brief    Encapsulate determination of timestep

public: // interface

  /// Create a new MethodEnzoTimestep
  MethodEnzoTimestep(EnzoDescr * enzo) throw();

public: // virtual functions

  /// Perform any method-specific initialization

  void initialize (DataDescr * data_descr) throw();

  /// Perform any timestep-specific finalizations steps, e.g. to
  /// deallocate any dynamically-allocated memory

  void finalize (DataDescr * data_descr) throw();

  /// Initialize PPM variable that may change.  Called once per
  /// block per timestep.

  void initialize_block (DataBlock * data_block) throw();

  /// Finalize PPM after advancing a block a timestep, e.g. to deallocate
  /// any dynamically-allocated variables

  void finalize_block (DataBlock * data_block) throw();

  /// Apply the timestep to advance a block one timestep 

  double compute_block( DataBlock * data_block ) throw(); 

protected: // functions

  float * pressure_field_;
  float afloat_;
  float dtBaryons_;
  float dtViscous_;
  float dtExpansion_;

  EnzoDescr * enzo_;
};

#endif /* USER_METHOD_ENZO_TIMESTEP_HPP */
