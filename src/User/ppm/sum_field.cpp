// $Id$
// See LICENSE_ENZO file for license and copyright information

/// @file      sum_field.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Sat Aug 29 23:45:03 PDT 2009
/// @brief     Return the sum of the internal values of the grid array

#include "enzo.hpp"
#include "cello_hydro.h"
 
float EnzoDescr::sum_field (int field)
{
  if (BaryonField[field] == NULL) return -1;
  float sum = 0.0;

  int ndx = GridDimension[0];
  int ndy = GridDimension[1];

  for (int iz = GridStartIndex[2]; iz<=GridEndIndex[2]; iz++) {
    for (int iy = GridStartIndex[1]; iy<=GridEndIndex[1]; iy++) {
      for (int ix = GridStartIndex[0]; ix<=GridEndIndex[0]; ix++) {
	int i = ix + ndx * (iy + ndy * iz);
	sum += BaryonField[field][i];
      }
    }
  }

  return sum;
}
    
