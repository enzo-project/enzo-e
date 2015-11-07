// See LICENSE_CELLO file for license and copyright information

/// @file     data_ItIndexRange.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-16
/// @brief    [\ref Data] Implementation of the ItIndexRange class

#include "cello.hpp"
#include "data.hpp"

//----------------------------------------------------------------------

void ItIndexRange::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change
  ItIndex::pup(p);
  p | index_;
  p | first_;
  p | last_;

}
