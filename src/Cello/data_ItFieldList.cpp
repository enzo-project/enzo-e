// See LICENSE_CELLO file for license and copyright information

/// @file     data_ItFieldList.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-16
/// @brief    

#include "cello.hpp"
#include "data.hpp"

//----------------------------------------------------------------------

void ItFieldList::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change
  ItField::pup(p);
  p | index_;
  p | values_;

}
