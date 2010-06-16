// $Id: user_MethodEnzoPpm.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     user_MethodEnzoPpm.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    Implementation of Enzo PPM hydro method

#ifndef USER_ENZO_PPM_HPP
#define USER_ENZO_PPM_HPP

#include <vector>
#include <string>
#include "method.hpp"
#include "data.hpp"

class MethodEnzoPpm : public MethodHyperbolic {

  /// @class    MethodEnzoPpm
  /// @ingroup  Enzo
  /// @brief    Encapsulate Enzo's PPM hydro method

public: // interface

  MethodEnzoPpm()
    : MethodHyperbolic()
  {};

  /// Perform any method-specific initialization

  void initialize_method (DataDescr * data_descr) throw();

  /// Perform any method-specific finalizations steps, e.g. to
  /// deallocate any dynamically-allocated memory

  void finalize_method(DataDescr * data_descr) throw();

  /// Initialize PPM variable that may change.  Called once per
  /// block per timestep.

  void initialize_block (DataBlock * data_block) throw();

  /// Finalize PPM after advancing a block a timestep, e.g. to deallocate
  /// any dynamically-allocated variables

  void finalize_block (DataBlock * data_block) throw();
  /// Apply the method to advance a block one timestep 

  void advance_block(DataBlock * data_block,
		     double t, double dt) throw();

  /// Refresh a block face's boundary / ghost zones given neighboring
  /// block face(s)

  void refresh_face() throw();

};

#endif /* USER_ENZO_PPM_HPP */
