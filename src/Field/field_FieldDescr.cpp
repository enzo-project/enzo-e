// $Id: field_FieldDescr.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the FieldDescr class

#include "error.hpp"
#include "field.hpp"

//----------------------------------------------------------------------

FieldDescr::FieldDescr
() throw ()
  : field_count_(0),
    group_count_(0),
    alignment_(0),
    padding_(0),
    courant_(1),
    field_name_(),
    group_name_(),
    precision_(),
    centering_(),
    ghosts_(),
    min_value_(),
    max_value_(),
    min_action_(),
    max_action_()
{
}

