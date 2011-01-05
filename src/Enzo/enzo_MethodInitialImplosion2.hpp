// $Id: enzo_MethodInitialImplosion2.hpp 1877 2010-11-30 01:20:27Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_METHODINITIALIMPLOSION2_HPP
#define ENZO_METHODINITIALIMPLOSION2_HPP

/// @file     enzo_MethodInitialImplosion2.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    Initialization routine for 2D implosion problem

class MethodInitialImplosion2 : public MethodInitial {

  /// @class    MethodInitialImplosion2
  /// @ingroup  Enzo
  /// @brief    Initialization routine for 2D implosion problem

public: // interface

  /// Constructor
  MethodInitialImplosion2(Global    * global,
			  EnzoDescr * enzo_descr) throw();

  /// Perform any method-specific initialization

  virtual void initialize (DataDescr * data_descr) throw();

  /// Perform any method-specific finalizations steps, e.g. to
  /// deallocate any dynamically-allocated memory

  virtual void finalize (DataDescr * data_descr) throw();

  /// Initialize PPM variable that may change.  Called once per
  /// block per timestep.

  virtual void initialize_block (DataBlock * data_block) throw();

  /// Finalize PPM after advancing a block a timestep, e.g. to deallocate
  /// any dynamically-allocated variables

  virtual void finalize_block (DataBlock * data_block) throw();

  /// Return the name of the method

  virtual std::string method_name() const throw();

private: // attributes

  EnzoDescr * enzo_descr_;
  
};

#endif /* ENZO_METHODINITIALIMPLOSION2_HPP */

