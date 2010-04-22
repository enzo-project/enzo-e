// $Id: field_field.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     field_field.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file field_field.cpp

#include "error.hpp"
#include "field.hpp"

//----------------------------------------------------------------------

Field::Field()
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

//----------------------------------------------------------------------

Field::~Field()
{
  if (centering_) delete [] centering_;
}

//----------------------------------------------------------------------

Field::Field(const Field & field) throw()
{
  INCOMPLETE_MESSAGE("Field::Field","");
}

//----------------------------------------------------------------------

Field & Field::operator= (const Field & field) throw()
{
  INCOMPLETE_MESSAGE("Field::operator =","");
  return *this;
}
