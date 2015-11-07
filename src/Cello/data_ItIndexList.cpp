// See LICENSE_CELLO file for license and copyright information

/// @file     data_ItIndexList.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-16
/// @brief    

#include "cello.hpp"
#include "data.hpp"

//----------------------------------------------------------------------

void ItIndexList::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change
  ItIndex::pup(p);
  p | index_;
  p | values_;

}
