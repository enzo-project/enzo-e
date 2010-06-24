// $Id: user_MethodEnzoPpm.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     user_MethodEnzoPpm.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    Implementation of Enzo PPM hydro method

#ifndef USER_METHOD_ENZO_PPM_HPP
#define USER_METHOD_ENZO_PPM_HPP

#include <vector>
#include <string>
#include "user.hpp"
#include "data.hpp"

class MethodEnzoPpm : public UserMethod {

  /// @class    MethodEnzoPpm
  /// @ingroup  Enzo
  /// @brief    Encapsulate Enzo's PPM hydro method

public: // interface

  MethodEnzoPpm()
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

  /// Refresh a block face's boundary / ghost zones given neighboring
  /// block face(s)

  void refresh_face() throw();
  
  /// Return the name of the method

  virtual std::string method_name() const throw() 
  { return "ppm"; };
  
};

#endif /* USER_METHOD_ENZO_PPM_HPP */
