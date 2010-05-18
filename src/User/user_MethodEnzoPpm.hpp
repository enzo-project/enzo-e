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

class MethodEnzoPpm : public MethodHyperbolic {

  /// @class    MethodEnzoPpm
  /// @ingroup  Enzo
  /// @brief    Encapsulate Enzo's PPM hydro method

public: // interface

  /// Perform any method-specific initialization

  void initialize() throw();

  /// Apply the method to advance a block one timestep 

  void advance_block() throw();

  /// Refresh a block face's boundary / ghost zones given neighboring block face(s) 

  void refresh_face() throw();

private: // attributes

  bool diffusion_;
  bool flattening_;
  bool steepening_; 

};

#endif /* USER_ENZO_PPM_HPP */
