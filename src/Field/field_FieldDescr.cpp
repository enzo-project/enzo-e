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
(
 std::string name,
 int dim
 ) throw ()
  : dim_(dim),
    name_(name),
    centering_(0),
    min_value_(0),
    max_value_(0),
    min_action_ (field_action_none),
    max_action_ (field_action_none),
    precision_  (default_precision_())
{
  centering_ = new bool [dim_];
  for (int i=0; i < dim_; i++) centering_[i] = true;
}

