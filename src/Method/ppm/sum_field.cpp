//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2009 James Bordner
 * Copyright (C) 2009 Laboratory for Computational Astrophysics
 * Copyright (C) 2009 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */


/** 
 *********************************************************************
 *
 * @file      sum_grid.cpp
 * @brief     Return the sum of the internal values of the grid array
 * @author    James Bordner
 * @date      Sat Aug 29 23:45:03 PDT 2009
 *
 * DESCRIPTION 
 * 
 *    Return the sum of the internal values of the grid array
 *
 * PACKAGES
 *
 *    NONE
 * 
 * INCLUDES
 *  
 *    cello_hydro.h
 *
 * PUBLIC FUNCTIONS
 *  
 *    NONE
 *
 * PRIVATE FUCTIONS
 *  
 *    NONE
 *
 * $Id$
 *
 *********************************************************************
 */

#include "cello_hydro.h"
 
float sum_field (int field)

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
    
