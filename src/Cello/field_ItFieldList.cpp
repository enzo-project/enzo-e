// See LICENSE_CELLO file for license and copyright information

/// @file     field_ItFieldList.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-16
/// @brief    

#include "cello.hpp"
#include "field.hpp"

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void ItFieldList::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change
  ItField::pup(p);
  p | index_;
  p | values_;

}

#endif
