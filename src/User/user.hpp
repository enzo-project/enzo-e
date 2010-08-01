// $Id: user.hpp 1300 2010-03-13 03:22:50Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef USER_HPP
#define USER_HPP

/// @file     user.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Include file for the User component

#include <vector>
#include <string>

#include "cello.hpp"

#include "data.hpp"
#include "error.hpp"
#include "global.hpp"
#include "parameters.hpp"

#include "enzo.hpp"

#include "user_UserControl.hpp"
#include "user_UserTimestep.hpp"
#include "user_UserMethod.hpp"
#include "user_UserDescr.hpp"

// Enzo

#include "user_MethodEnzoControl.hpp"
#include "user_MethodEnzoTimestep.hpp"
#include "user_MethodEnzoPpm.hpp"
#include "user_MethodEnzoPpml.hpp"
#include "user_EnzoUserDescr.hpp"

#endif /* USER_HPP */

