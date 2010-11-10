// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file     enzo_MethodEnzoPpm.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    Implementation of Enzo PPM hydro method

#ifndef ENZO_METHOD_ENZO_PPM_HPP
#define ENZO_METHOD_ENZO_PPM_HPP

class MethodEnzoPpm : public MethodHyperbolic {

  /// @class    MethodEnzoPpm
  /// @ingroup  Enzo
  /// @brief    Encapsulate Enzo's PPM hydro method

public: // interface

  MethodEnzoPpm(Global * global,
		EnzoDescr * enzo):
    MethodHyperbolic(global),
    enzo_(enzo)
  {};

  /// Perform any method-specific initialization

  void initialize (DataDescr * data_descr) throw();

  /// Perform any method-specific finalizations steps, e.g. to
  /// deallocate any dynamically-allocated memory

  void finalize (DataDescr * data_descr) throw();

  /// Initialize PPM variable that may change.  Called once per
  /// block per timestep.

  void initialize_block (DataBlock * data_block) throw();

  /// Finalize PPM after advancing a block a timestep, e.g. to deallocate
  /// any dynamically-allocated variables

  void finalize_block (DataBlock * data_block) throw();
  /// Apply the method to advance a block one timestep 

  void advance_block(DataBlock * data_block,
		     double t, double dt) throw();

  /// Return the name of the method

  std::string method_name() const throw() 
  { return "ppm"; };

private:

  EnzoDescr * enzo_;

};

#endif /* ENZO_METHOD_ENZO_PPM_HPP */
