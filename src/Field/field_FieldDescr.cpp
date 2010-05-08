// $Id: field_FieldDescr.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the FieldDescr class

#include "error.hpp"
#include "field.hpp"

//----------------------------------------------------------------------

FieldDescr::FieldDescr()
  : name_(),
    id_(0),
    dim_(0),
    block_number_(0),
    block_offset_(0),
    centering_(0),
    min_(0),
    max_(0),
    min_action_ (action_field_unknown),
    max_action_ (action_field_unknown),
    precision_  (precision_unknown)
{
}

