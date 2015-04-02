// See LICENSE_CELLO file for license and copyright information

/// @file     data_ItFieldRange.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-16
/// @brief    [\ref Data] Implementation of the ItFieldRange class

#include "cello.hpp"
#include "data.hpp"

//----------------------------------------------------------------------

void ItFieldRange::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change
  ItField::pup(p);
  p | index_;
  p | first_;
  p | last_;

}
